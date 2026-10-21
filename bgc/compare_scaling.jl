# Effect of MAB_BGC_SCALE_NEGATIVES: RUN_A (not scaled) and RUN_B (scaled, restarted from RUN_A's checkpoint) on the daily frames
# they share (the filtered output of a restart is valid from 3 days after it): negative ammonium, the nitrogen total, and the
# surface chlorophyll against the satellite and the Pioneer sensors.
#   julia -t 4 --project=<env with CairoMakie> bgc/compare_scaling.jl
# Environment: RUN_A, RUN_B, MAB_START_DATE, MAB_CELLS_PER_DEGREE, BGC_SAT, BGC_OOI, BGC_OUT.
using Printf, Statistics, Dates
using CairoMakie
using Oceananigans
using Oceananigans.Grids: znodes
using NumericalEarth.DataWrangling: jldopen
const NCD = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "NCDatasets")]

const RUN_A = get(ENV, "RUN_A", "/t0/workdir/enrique/runs/bgc/bgc_feb19_rc4")
const RUN_B = get(ENV, "RUN_B", "/t0/workdir/enrique/runs/bgc/bgc_feb19_rc4sn")
const START = DateTime(get(ENV, "MAB_START_DATE", "2019-02-01"))
const CPD   = parse(Int, get(ENV, "MAB_CELLS_PER_DEGREE", "12"))
const SAT   = get(ENV, "BGC_SAT", "/t0/workdir/enrique/Data/ocean_colour/chl_l3_multi_4km_mab_2019.nc")
const OOI   = get(ENV, "BGC_OOI", "/t0/workdir/enrique/Data/OOI/pioneer_bgc_2019")
const OUT   = get(ENV, "BGC_OUT", "/t0/workdir/enrique/runs/analysis/scaling_" * basename(RUN_B) * "_vs_" * basename(RUN_A))
const H, Nz = 7, parse(Int, get(ENV, "MAB_NZ", "50"))
const CHL_PER_N = 1.31
const λ_west, φ_south = -76 + 2 / 12, 34 + 2 / 12
mkpath(OUT)
io = open(joinpath(OUT, "summary.txt"), "w")
say(args...) = (println(args...); flush(stdout); println(io, args...); flush(io))

function exponential_scale(N, H, Δtop)
    top(h) = H * expm1(H / N / h) / expm1(H / h)
    lo, hi = H / 50, 100H
    for _ in 1:200; mid = (lo + hi) / 2; top(mid) > Δtop ? (hi = mid) : (lo = mid); end
    return (lo + hi) / 2
end
Δz_top, Z_BOTTOM = parse(Float64, get(ENV, "MAB_DZ_TOP", "2")), parse(Float64, get(ENV, "MAB_ZBOTTOM", "5500"))
zgrid = RectilinearGrid(size = (1, 1, Nz), x = (0, 1), y = (0, 1), z = ExponentialDiscretization(Nz, -Z_BOTTOM, 0; scale = exponential_scale(Nz, Z_BOTTOM, Δz_top)))
const Δz = diff(collect(znodes(zgrid, Face())))
bottom_files(tag) = sort(filter(f -> occursin(Regex("^" * basename(tag) * "_bottom_rank\\d+\\.jld2\$"), basename(f)), readdir(dirname(tag); join = true)))
bf = isempty(bottom_files(RUN_B)) ? bottom_files(RUN_A) : bottom_files(RUN_B)
parts = map(bf) do f
    jldopen(f) do file; (I = file["I_OFF"], J = file["J_OFF"], kb = file["k_bottom"], dzb = file["Δz_bottom"]); end
end
const Nx = maximum(p.I + size(p.kb, 1) for p in parts); const Ny = maximum(p.J + size(p.kb, 2) for p in parts)
kb = zeros(Int, Nx, Ny); depth = zeros(Nx, Ny)
for p in parts, j in axes(p.kb, 2), i in axes(p.kb, 1)
    k = p.kb[i, j]; kb[p.I + i, p.J + j] = k
    depth[p.I + i, p.J + j] = k > 0 ? sum(Δz[k+1:end]) + p.dzb[i, j] : 0.0
end
wet = kb .> 0
λm(i) = λ_west + (i - 0.5) / CPD; φm(j) = φ_south + (j - 0.5) / CPD
regions = (shelf = (depth .> 0) .& (depth .< 100), slope = (depth .≥ 100) .& (depth .< 1000), deep = depth .≥ 1000)
cellweight = [cosd(φm(j)) for i in 1:Nx, j in 1:Ny, k in 1:Nz] .* reshape(Δz, 1, 1, Nz)      # relative cell volumes (static grid)

iters(f, name) = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/$name"]))))
strip3(A) = A[H+1:end-H, H+1:end-H, H+1:end-H]
function frame3(tag, name, it; stem = "bgc_volume_daily")
    A = zeros(Float32, Nx, Ny, Nz)
    for (r, p) in enumerate(parts)
        jldopen("$(tag)_$(stem)_rank$(r - 1).jld2") do f
            B = strip3(f["timeseries/$name/$it"]); A[p.I+1:p.I+size(B, 1), p.J+1:p.J+size(B, 2), :] .= B
        end
    end
    return A
end
function frame_index(tag)       # time (s) → iteration key of the daily volume frames
    jldopen("$(tag)_bgc_volume_daily_rank0.jld2") do f
        Dict(round(Int, f["timeseries/t/$it"]) => it for it in iters(f, "P"))
    end
end
ia, ib = frame_index(RUN_A), frame_index(RUN_B)
common = sort(collect(intersect(keys(ia), keys(ib))))
isempty(common) && error("the two runs share no daily frame")
dates = [START + Second(t) for t in common]
say("A = $(basename(RUN_A)) (not scaled), B = $(basename(RUN_B)) (scaled): $(length(common)) shared daily frames, $(Date(first(dates))) → $(Date(last(dates)))")

# satellite
ds = NCD.NCDataset(SAT)
λs = Float64.(ds["longitude"][:]); φs = Float64.(ds["latitude"][:]); sat_times = ds["time"][:]
raw = ds["CHL"][:, :, :]; chl_sat = Array{Float32}(undef, size(raw)); for k in eachindex(raw); chl_sat[k] = ismissing(raw[k]) ? NaN32 : Float32(raw[k]); end
close(ds)
sat_day = Dict(Date(t) => n for (n, t) in enumerate(sat_times))
lon_idx = [findall(abs.(λs .- λm(i)) .< 0.5 / CPD) for i in 1:Nx]; lat_idx = [findall(abs.(φs .- φm(j)) .< 0.5 / CPD) for j in 1:Ny]
function sat_cell(n, i, j)
    v = [x for x in chl_sat[lon_idx[i], lat_idx[j], n] if isfinite(x) && x > 0]
    return length(v) ≥ 2 ? Float64(mean(v)) : NaN
end

nitrogen = ("NO₃", "NH₄", "P", "Z", "DOM", "sPOM", "bPOM")
rows = Dict(tag => Dict(k => fill(NaN, length(common)) for k in (:negfrac, :negtop, :minNH4, :N, :NO3, :NH4, :chl_shelf, :chl_slope, :chl_deep, :maxNO3, :minO2)) for tag in (:A, :B))
matched = Dict((tag, r) => (Float64[], Float64[]) for tag in (:A, :B), r in keys(regions))
for (n, t) in enumerate(common)
    d = Date(dates[n]); satn = get(sat_day, d, 0)
    s_map = satn > 0 ? [wet[i, j] ? sat_cell(satn, i, j) : NaN for i in 1:Nx, j in 1:Ny] : nothing
    S = frame3(RUN_A, "S", ia[t]; stem = "volume_daily"); w3 = S .> 0
    for (tag, run, idx) in ((:A, RUN_A, ia), (:B, RUN_B, ib))
        F = Dict(name => frame3(run, name, idx[t]) for name in (nitrogen..., "O₂"))
        NH4 = F["NH₄"]
        rows[tag][:negfrac][n] = count(<(-1e-3), NH4[w3]) / count(w3)
        rows[tag][:negtop][n] = count(<(-0.1), NH4[:, :, Nz][wet]) / count(wet)
        rows[tag][:minNH4][n] = minimum(NH4[w3])
        rows[tag][:N][n] = sum(sum(F[x] .* cellweight .* w3) for x in nitrogen)
        rows[tag][:NO3][n] = sum(F["NO₃"] .* cellweight .* w3); rows[tag][:NH4][n] = sum(NH4 .* cellweight .* w3)
        rows[tag][:maxNO3][n] = maximum(F["NO₃"][w3]); rows[tag][:minO2][n] = minimum(F["O₂"][w3])
        chl = CHL_PER_N .* F["P"][:, :, Nz]
        for (r, m) in pairs(regions)
            rows[tag][Symbol(:chl_, r)][n] = mean(chl[m])
            if s_map !== nothing
                ok = m .& isfinite.(s_map) .& (chl .> 0)
                append!(matched[(tag, r)][1], s_map[ok]); append!(matched[(tag, r)][2], chl[ok])
            end
        end
    end
    say(@sprintf("  day %5.1f: NH₄ < 0: %.1f%% → %.1f%% of wet cells, surface columns < -0.1: %.1f%% → %.1f%%, min NH₄ %.2f → %.2f, N total B/A − 1 = %+.2e",
                 t / 86400, 100rows[:A][:negfrac][n], 100rows[:B][:negfrac][n], 100rows[:A][:negtop][n], 100rows[:B][:negtop][n],
                 rows[:A][:minNH4][n], rows[:B][:minNH4][n], rows[:B][:N][n] / rows[:A][:N][n] - 1))
end
last_ = length(common)
say(@sprintf("\nlast shared day: NO₃ total B/A − 1 = %+.2e, NH₄ total A %.3g B %.3g (volume-weighted, arbitrary units), max NO₃ %.1f → %.1f, min O₂ %.1f → %.1f",
             rows[:B][:NO3][last_] / rows[:A][:NO3][last_] - 1, rows[:A][:NH4][last_], rows[:B][:NH4][last_],
             rows[:A][:maxNO3][last_], rows[:B][:maxNO3][last_], rows[:A][:minO2][last_], rows[:B][:minO2][last_]))
say("\nsurface chlorophyll against the satellite over the shared days (mg/m³; statistics of log10): run A / run B")
for r in keys(regions)
    line = @sprintf("  %-6s", r)
    for tag in (:A, :B)
        o, m = matched[(tag, r)]; ok = (o .> 0) .& (m .> 0); o, m = o[ok], m[ok]
        length(o) < 5 && continue
        line *= @sprintf("  %s: n=%6d model %.2f (sat %.2f) bias %+.2f rms %.2f r %+.2f ", tag, length(o), mean(m), mean(o), mean(log10.(m) .- log10.(o)),
                         sqrt(mean((log10.(m) .- log10.(o)) .^ 2)), cor(log10.(o), log10.(m)))
    end
    say(line)
end
say("\nregion-mean surface chlorophyll over the shared days (mg/m³), A → B: " * join([@sprintf("%s %.3f → %.3f", r, mean(rows[:A][Symbol(:chl_, r)]), mean(rows[:B][Symbol(:chl_, r)])) for r in keys(regions)], ", "))

days = [t / 86400 for t in common]
fig = Figure(size = (1500, 450))
ax1 = Axis(fig[1, 1], title = "wet cells with NH₄ < 0 (%)", xlabel = "days from $(Date(START))")
lines!(ax1, days, 100 .* rows[:A][:negfrac]; label = "A, not scaled", color = :black); lines!(ax1, days, 100 .* rows[:B][:negfrac]; label = "B, scaled", color = :tomato); axislegend(ax1; position = :lt)
ax2 = Axis(fig[1, 2], title = "surface columns with NH₄ < −0.1 (%)", xlabel = "days from $(Date(START))")
lines!(ax2, days, 100 .* rows[:A][:negtop]; color = :black); lines!(ax2, days, 100 .* rows[:B][:negtop]; color = :tomato)
ax3 = Axis(fig[1, 3], title = "region-mean surface chlorophyll (mg/m³)", xlabel = "days from $(Date(START))", yscale = log10)
for (r, c) in zip(keys(regions), (:steelblue, :darkgreen, :purple))
    lines!(ax3, days, rows[:A][Symbol(:chl_, r)]; color = c, linestyle = :dash); lines!(ax3, days, rows[:B][Symbol(:chl_, r)]; color = c, label = string(r))
end
axislegend(ax3; position = :lt)
save(joinpath(OUT, "scaling.png"), fig)
say("\nfigures and summary in $OUT")
close(io)
