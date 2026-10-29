# Effect of the sediment model on the BGC fields: compares runs of bgc/mab_bgc.jl that differ in MAB_BGC_SEDIMENT (restarted from
# the same checkpoint). For the de-tided daily frames the runs have in common it prints, for the shelf (< 100 m), slope
# (100-1000 m) and deep (> 1000 m) regions, the mean of: surface chlorophyll (1.31 P), the phytoplankton in the top 50 m (mmol N/m²),
# NH₄, NO₃ and O₂ in the surface level and in the bottom wet cell of each column, and the DIC of the bottom cell; and, for
# multi-G, the sediment pools of the shelf columns and the sinking nitrogen flux they receive (<TAG>_bgc_sediment_rank<r>.jld2).
# Changes are relative to the first run.
#
#   BGC_TAGS="instant=/t0/.../tagA,multig=/t0/.../tagB" julia --project=<env> bgc/sediment_comparison.jl
using Printf, Statistics, Dates
using Oceananigans
using Oceananigans.Grids: znodes
using NumericalEarth.DataWrangling: jldopen

const RUNS = [(split(kv, "=")[1], split(kv, "=")[2]) for kv in split(ENV["BGC_TAGS"], ",")]
const H, Nz = 7, parse(Int, get(ENV, "MAB_NZ", "50"))
const CHL_PER_N = 1.31

function exponential_scale(N, H, Δtop)
    top(h) = H * expm1(H / N / h) / expm1(H / h)
    lo, hi = H / 50, 100H
    for _ in 1:200; mid = (lo + hi) / 2; top(mid) > Δtop ? (hi = mid) : (lo = mid); end
    return (lo + hi) / 2
end
Δz_top, Z_BOTTOM = parse(Float64, get(ENV, "MAB_DZ_TOP", "2")), parse(Float64, get(ENV, "MAB_ZBOTTOM", "5500"))
zgrid = RectilinearGrid(size = (1, 1, Nz), x = (0, 1), y = (0, 1),
                        z = ExponentialDiscretization(Nz, -Z_BOTTOM, 0; scale = exponential_scale(Nz, Z_BOTTOM, Δz_top)))
const zf = collect(znodes(zgrid, Face())); const Δz = diff(zf)
const K50 = findall(k -> -zf[k+1] < 50, 1:Nz)

iters(f, name) = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/$name"]))))
bottom_files = sort(filter(f -> occursin(Regex("^" * basename(last(first(RUNS))) * "_bottom_rank\\d+\\.jld2\$"), basename(f)), readdir(dirname(last(first(RUNS))); join = true)))
const NR = length(bottom_files)
parts = map(bottom_files) do f
    jldopen(f) do file; (I = file["I_OFF"], J = file["J_OFF"], kb = file["k_bottom"], dzb = file["Δz_bottom"]); end
end
const Nx = maximum(p.I + size(p.kb, 1) for p in parts); const Ny = maximum(p.J + size(p.kb, 2) for p in parts)
kb = zeros(Int, Nx, Ny); depth = zeros(Nx, Ny)
for p in parts, j in axes(p.kb, 2), i in axes(p.kb, 1)
    k = p.kb[i, j]; kb[p.I + i, p.J + j] = k
    depth[p.I + i, p.J + j] = k > 0 ? sum(Δz[k+1:end]) + p.dzb[i, j] : 0.0
end
wet = kb .> 0
regions = (shelf = (depth .> 0) .& (depth .< 100), slope = (depth .≥ 100) .& (depth .< 1000), deep = depth .≥ 1000)

function frame(tag, name, n)
    A = zeros(Float32, Nx, Ny, Nz)
    for (p, r) in zip(parts, 0:NR-1)
        jldopen("$(tag)_bgc_volume_daily_rank$(r).jld2") do f
            it = iters(f, name)[n]
            B = f["timeseries/$name/$it"][H+1:end-H, H+1:end-H, H+1:end-H]
            A[p.I+1:p.I+size(B, 1), p.J+1:p.J+size(B, 2), :] .= B
        end
    end
    return A
end
frame_times(tag) = jldopen(f -> [f["timeseries/t/$it"] for it in iters(f, "P")], "$(tag)_bgc_volume_daily_rank0.jld2")
times = Dict(label => round.(Int, frame_times(tag)) for (label, tag) in RUNS)
common = sort(collect(intersect(Set.(values(times))...)))
println("common daily frames: $(length(common)), day $(round(first(common) / 86400, digits = 1)) to $(round(last(common) / 86400, digits = 1)) of the run clock")
isempty(common) && error("no common frames")

bottom(A) = [kb[i, j] > 0 ? A[i, j, kb[i, j]] : NaN32 for i in 1:Nx, j in 1:Ny]
top(A) = A[:, :, Nz]
function regional_means(tag, label)
    out = Dict{Tuple{Symbol, Symbol}, Vector{Float64}}()
    for t in common
        n = argmin(abs.(times[label] .- t))
        P, NH4, NO3, O2, DIC = (frame(tag, v, n) for v in ("P", "NH₄", "NO₃", "O₂", "DIC"))
        fields = (chl = top(P) .* CHL_PER_N,
                  P50 = sum(P[:, :, K50[i]] .* Δz[K50[i]] for i in eachindex(K50)),
                  NH4_top = top(NH4), NH4_bot = bottom(NH4), NO3_top = top(NO3), NO3_bot = bottom(NO3),
                  O2_top = top(O2), O2_bot = bottom(O2), DIC_bot = bottom(DIC))
        for (rname, mask) in pairs(regions), (q, A) in pairs(fields)
            v = A[mask .& wet .& isfinite.(A)]
            push!(get!(out, (rname, q), Float64[]), isempty(v) ? NaN : mean(v))
        end
    end
    return out
end
results = [(label, regional_means(tag, label)) for (label, tag) in RUNS]
ref = first(results)[2]
labels = Dict(:chl => "surface chlorophyll (mg/m³)", :P50 => "phytoplankton, top 50 m (mmol N/m²)", :NH4_top => "NH₄, surface level (mmol/m³)", :NH4_bot => "NH₄, bottom cell (mmol/m³)",
              :NO3_top => "NO₃, surface level (mmol/m³)", :NO3_bot => "NO₃, bottom cell (mmol/m³)", :O2_top => "O₂, surface level (mmol/m³)",
              :O2_bot => "O₂, bottom cell (mmol/m³)", :DIC_bot => "DIC, bottom cell (mmol/m³)")
for q in (:chl, :P50, :NH4_top, :NH4_bot, :NO3_top, :NO3_bot, :O2_top, :O2_bot, :DIC_bot)
    println("\n", labels[q], ": mean over the common frames (change against ", first(results)[1], ")")
    println(@sprintf("%-12s", ""), join([@sprintf("%22s", String(r)) for r in keys(regions)]))
    for (label, r) in results
        row = @sprintf("%-12s", label)
        for rname in keys(regions)
            v = mean(filter(isfinite, r[(rname, q)])); v0 = mean(filter(isfinite, ref[(rname, q)]))
            row *= @sprintf("%12.4f (%+6.1f%%)", v, 100 * (v / v0 - 1))
        end
        println(row)
    end
end

# sediment pools of a run that has them: mean over the shelf columns and over the days
for (label, tag) in RUNS
    isfile("$(tag)_bgc_sediment_rank0.jld2") || continue
    names = jldopen(f -> sort(filter(k -> k != "t", collect(keys(f["timeseries"])))), "$(tag)_bgc_sediment_rank0.jld2")
    println("\nsediment fields of $label, shelf columns (< 100 m), first and last saved day: ", join(names, ", "))
    for name in names
        vals = Float64[]; first_last = Float64[]
        series = Vector{Vector{Float64}}()
        for (p, r) in zip(parts, 0:NR-1)
            jldopen("$(tag)_bgc_sediment_rank$(r).jld2") do f
                its = iters(f, name)
                for (k, it) in enumerate(its)
                    B = f["timeseries/$name/$it"]
                    nx, ny = size(p.kb)
                    size(B, 1) == nx + 2H && (B = B[H+1:end-H, H+1:end-H, :])
                    mask = regions.shelf[p.I+1:p.I+nx, p.J+1:p.J+ny]
                    length(series) < k && push!(series, Float64[])
                    append!(series[k], B[:, :, 1][mask])
                end
            end
        end
        means = [isempty(s) ? NaN : mean(s) for s in series]
        println(@sprintf("  %-14s day 1 %10.4g   day %d %10.4g   (all days: %s)", name, means[1], length(means), means[end], join([@sprintf("%.3g", m) for m in means], " ")))
    end
end
