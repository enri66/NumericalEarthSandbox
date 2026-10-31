# Benthic nitrogen flux of a run with the multi-G sediment (MAB_BGC_SEDIMENT=multig) by water depth, to compare with the bottom flux F of the
# single-column tests (bgc/column_compare.jl, MAB_COLUMN_BOTTOM_FLUX, mmol N/m²/day). For the sediment's daily means it prints, per depth bin:
#   - the sinking flux of sPOM + bPOM into the sediment (mmol N/m²/day),
#   - the return to the water, λs Ns + λf Nf (0.2 and 2.0 per day, as in bgc/nbudget_analysis.jl; the refractory pool Nr does not return),
#   - the pools Ns, Nf, Nr (mmol N/m²), and the number of columns.
#   BGC_TAG=/t0/.../bgc_feb19_woa_lim_mg30_fix julia --project=<env> bgc/benthic_flux_by_depth.jl
# Environment: BGC_TAG, BGC_OUT, BGC_FROM_DAY (first sediment day averaged, 5: the pools are spun up), MAB_NZ (50).
using Printf, Statistics
using Oceananigans
using Oceananigans.Grids: znodes
using NumericalEarth.DataWrangling: jldopen

const TAG  = ENV["BGC_TAG"]
const OUT  = get(ENV, "BGC_OUT", "/t0/workdir/enrique/runs/analysis/benthic_flux_" * basename(TAG))
const FROM = parse(Float64, get(ENV, "BGC_FROM_DAY", "5"))
const H, Nz = 7, parse(Int, get(ENV, "MAB_NZ", "50"))
const BINS = [(0.0, 10.0), (10.0, 20.0), (20.0, 40.0), (40.0, 60.0), (60.0, 100.0)]
mkpath(OUT)
io = open(joinpath(OUT, "summary.txt"), "w")
function say(args...); println(args...); flush(stdout); println(io, args...); flush(io); end

bottom_files = ["$(TAG)_bottom_rank$(r).jld2" for r in 0:99 if isfile("$(TAG)_bottom_rank$(r).jld2")]
parts = map(bottom_files) do f; jldopen(f) do file; (I = file["I_OFF"], J = file["J_OFF"], kb = file["k_bottom"], dzb = file["Δz_bottom"]); end; end
iters(f, name) = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/$name"]))))
function exponential_scale(N, H, Δtop)
    top(h) = H * expm1(H / N / h) / expm1(H / h)
    lo, hi = H / 50, 100H
    for _ in 1:200; mid = (lo + hi) / 2; top(mid) > Δtop ? (hi = mid) : (lo = mid); end
    return (lo + hi) / 2
end
zg = RectilinearGrid(size = (1, 1, Nz), x = (0, 1), y = (0, 1), z = ExponentialDiscretization(Nz, -5500, 0; scale = exponential_scale(Nz, 5500, 2.0)))
Δz = diff(collect(znodes(zg, Face())))
Nx = maximum(p.I + size(p.kb, 1) for p in parts); Ny = maximum(p.J + size(p.kb, 2) for p in parts)
depth = zeros(Nx, Ny)
for p in parts, j in axes(p.kb, 2), i in axes(p.kb, 1)
    k = p.kb[i, j]; depth[p.I + i, p.J + j] = k > 0 ? sum(Δz[k+1:end]) + p.dzb[i, j] : 0.0
end
function read2(name, n)
    A = zeros(Nx, Ny)
    for (p, r) in zip(parts, 0:length(parts)-1)
        jldopen("$(TAG)_bgc_sediment_rank$(r).jld2") do f
            it = iters(f, name)[n]; B = f["timeseries/$name/$it"]
            nx, ny = size(p.kb); hx = (size(B, 1) - nx) ÷ 2; hy = (size(B, 2) - ny) ÷ 2
            A[p.I+1:p.I+nx, p.J+1:p.J+ny] .= B[hx+1:hx+nx, hy+1:hy+ny, 1]
        end
    end
    return A
end
st = jldopen(f -> [f["timeseries/t/$it"] for it in iters(f, "Ns")], "$(TAG)_bgc_sediment_rank0.jld2")
days = st ./ 86400
sel = [n for n in eachindex(st) if days[n] ≥ FROM + days[1] - 1]
isempty(sel) && (sel = collect(eachindex(st)))
say("$TAG: $(length(st)) sediment days, averaging $(length(sel)) (day $(round(days[first(sel)], digits = 1)) to $(round(days[last(sel)], digits = 1)))")
acc = Dict(b => (sink = Float64[], ret = Float64[], Ns = Float64[], Nf = Float64[], Nr = Float64[], cols = Ref(0)) for b in BINS)
for n in sel
    Ns, Nf, Nr = read2("Ns", n), read2("Nf", n), read2("Nr", n)
    fl = read2("tracked_sPOM", n) .+ read2("tracked_bPOM", n)             # mmol N/m²/s
    ret = 0.2 .* Ns .+ 2.0 .* Nf                                          # mmol N/m²/day
    for b in BINS
        m = (depth .> b[1]) .& (depth .≤ b[2])
        acc[b].cols[] = count(m)
        count(m) == 0 && continue
        push!(acc[b].sink, mean(fl[m]) * 86400); push!(acc[b].ret, mean(ret[m]))
        push!(acc[b].Ns, mean(Ns[m])); push!(acc[b].Nf, mean(Nf[m])); push!(acc[b].Nr, mean(Nr[m]))
    end
end
say("\nbenthic nitrogen flux by water depth, mean over the days (mmol N/m²/day; pools mmol N/m²)")
say("  depth (m)   columns   sinking flux   return to water   Ns     Nf     Nr")
for b in BINS
    a = acc[b]
    isempty(a.sink) && (say(@sprintf("  %3.0f-%-4.0f %8d   no columns", b[1], b[2], a.cols[])); continue)
    say(@sprintf("  %3.0f-%-4.0f %8d   %12.3f   %15.3f   %.2f  %.2f  %.2f", b[1], b[2], a.cols[], mean(a.sink), mean(a.ret), mean(a.Ns), mean(a.Nf), mean(a.Nr)))
end
