# What the river loads did: two BGC runs on the same grid and layout, RUN_A without river loads and RUN_B with them, on their
# de-tided daily frames. Surface (top-level) fields, mean difference maps, the Hudson and Delaware plume regions (within
# PLUME_KM of the mouths) and, in those regions, surface chlorophyll against the satellite for both runs.
#   julia -t 4 --project=<env with CairoMakie> bgc/compare_rivers.jl
# Environment: RUN_A, RUN_B (path prefixes), MAB_START_DATE, MAB_CELLS_PER_DEGREE, BGC_SAT, BGC_OUT, PLUME_KM, DAY_FROM (first day used).
using Printf, Statistics, Dates
using CairoMakie
using Oceananigans
using Oceananigans.Grids: znodes
using NumericalEarth.DataWrangling: jldopen
const NCD = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "NCDatasets")]

const RUN_A = get(ENV, "RUN_A", "/t0/workdir/enrique/runs/bgc/bgc_feb19_cmb")
const RUN_B = get(ENV, "RUN_B", "/t0/workdir/enrique/runs/bgc/bgc_feb19_rc4")
const START = DateTime(get(ENV, "MAB_START_DATE", "2019-02-01"))
const CPD   = parse(Int, get(ENV, "MAB_CELLS_PER_DEGREE", "12"))
const SAT   = get(ENV, "BGC_SAT", "/t0/workdir/enrique/Data/ocean_colour/chl_l3_multi_4km_mab_2019.nc")
const OUT   = get(ENV, "BGC_OUT", "/t0/workdir/enrique/runs/analysis/rivers_" * basename(RUN_B) * "_vs_" * basename(RUN_A))
const PLUME_KM = parse(Float64, get(ENV, "PLUME_KM", "100"))
const DAY_FROM = parse(Int, get(ENV, "DAY_FROM", "8"))          # skip the first days (filter start-up, plume spin-up)
const H, Nz = 7, parse(Int, get(ENV, "MAB_NZ", "50"))
const CHL_PER_N = 1.31
const λ_west, φ_south = -76 + 2 / 12, 34 + 2 / 12
const MOUTHS = ("Hudson" => (-73.95, 40.50), "Delaware" => (-75.30, 39.20))
mkpath(OUT)
io = open(joinpath(OUT, "summary.txt"), "w")
say(args...) = (println(args...); flush(stdout); println(io, args...); flush(io))

# ---------------- layout and depth (from run B's bottom files; run A must have the same layout) ----------------
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
parts = map(bottom_files(RUN_B)) do f
    jldopen(f) do file; (I = file["I_OFF"], J = file["J_OFF"], kb = file["k_bottom"], dzb = file["Δz_bottom"]); end
end
NR = length(parts); length(bottom_files(RUN_A)) in (0, NR) || error("the two runs have different numbers of ranks")
const Nx = maximum(p.I + size(p.kb, 1) for p in parts); const Ny = maximum(p.J + size(p.kb, 2) for p in parts)
kb = zeros(Int, Nx, Ny); depth = zeros(Nx, Ny)
for p in parts, j in axes(p.kb, 2), i in axes(p.kb, 1)
    k = p.kb[i, j]; kb[p.I + i, p.J + j] = k
    depth[p.I + i, p.J + j] = k > 0 ? sum(Δz[k+1:end]) + p.dzb[i, j] : 0.0
end
wet = kb .> 0
λm(i) = λ_west + (i - 0.5) / CPD; φm(j) = φ_south + (j - 0.5) / CPD
km(λ₁, φ₁, λ₂, φ₂) = 2 * 6371 * asin(sqrt(sind((φ₂ - φ₁) / 2)^2 + cosd(φ₁) * cosd(φ₂) * sind((λ₂ - λ₁) / 2)^2))

iters(f, name) = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/$name"]))))
strip3(A) = A[H+1:end-H, H+1:end-H, H+1:end-H]
function surface(tag, name, n; stem = "bgc_volume_daily")
    A = zeros(Float32, Nx, Ny)
    for (r, p) in enumerate(parts)
        jldopen("$(tag)_$(stem)_rank$(r - 1).jld2") do f
            B = strip3(f["timeseries/$name/$(iters(f, name)[n])"])[:, :, Nz]
            A[p.I+1:p.I+size(B, 1), p.J+1:p.J+size(B, 2)] .= B
        end
    end
    return A
end
frame_times = jldopen(f -> [f["timeseries/t/$it"] for it in iters(f, "P")], "$(RUN_B)_bgc_volume_daily_rank0.jld2")
frame_dates = [START + Millisecond(round(Int, 1000t)) for t in frame_times]
Ndays = length(frame_dates)
frames = [n for n in 1:Ndays if (n + 2) ≥ DAY_FROM]            # frame n is day n + 2 (the first valid frame is day 3)
say("runs: A = $(basename(RUN_A)) (no river loads), B = $(basename(RUN_B)); $(length(frames)) daily frames from $(Date(frame_dates[first(frames)]))")

masks = Dict(name => [wet[i, j] && km(λm(i), φm(j), λ, φ) < PLUME_KM for i in 1:Nx, j in 1:Ny] for (name, (λ, φ)) in MOUTHS)
masks["Hudson plume"] = masks["Hudson"]; masks["Delaware plume"] = masks["Delaware"]
tracers = (:NO₃, :NH₄, :DOM, :Alk, :DIC, :O₂)

# ---------------- satellite ----------------
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

# ---------------- one pass over the frames ----------------
sums = Dict((tag, t) => zeros(Nx, Ny) for tag in (:A, :B), t in (tracers..., :chl, :S))
series = Dict((tag, reg, t) => fill(NaN, Ndays) for tag in (:A, :B), reg in ("Hudson plume", "Delaware plume"), t in (:NO₃, :chl, :Alk, :S))
pairs_ = Dict((tag, reg) => (Float64[], Float64[]) for tag in (:A, :B), reg in ("Hudson plume", "Delaware plume"))
for n in frames
    d = Date(frame_dates[n]); satn = get(sat_day, d, 0)
    s_map = satn > 0 ? [wet[i, j] ? sat_cell(satn, i, j) : NaN for i in 1:Nx, j in 1:Ny] : nothing
    for (tag, run) in ((:A, RUN_A), (:B, RUN_B))
        f = Dict{Symbol, Matrix{Float32}}()
        for t in tracers; f[t] = surface(run, string(t), n); end
        f[:chl] = CHL_PER_N .* surface(run, "P", n); f[:S] = surface(run, "S", n; stem = "volume_daily")
        for (t, A) in f; sums[(tag, t)] .+= A; end
        for reg in ("Hudson plume", "Delaware plume"), t in (:NO₃, :chl, :Alk, :S)
            series[(tag, reg, t)][n] = mean(f[t][masks[reg]])
        end
        if s_map !== nothing
            for reg in ("Hudson plume", "Delaware plume")
                ok = masks[reg] .& isfinite.(s_map) .& (f[:chl] .> 0)
                append!(pairs_[(tag, reg)][1], s_map[ok]); append!(pairs_[(tag, reg)][2], f[:chl][ok])
            end
        end
    end
    n % 10 == 0 && say("  frame $n of $Ndays")
end
mean_(tag, t) = [wet[i, j] ? sums[(tag, t)][i, j] / length(frames) : NaN for i in 1:Nx, j in 1:Ny]

# ---------------- numbers ----------------
say(@sprintf("\nmean surface fields in the plume regions (within %.0f km of the mouth), days ≥ %d: run A / run B", PLUME_KM, DAY_FROM))
for reg in ("Hudson plume", "Delaware plume")
    cells = masks[reg]; say("$reg: $(count(cells)) cells")
    for t in (tracers..., :chl, :S)
        a, b = mean(mean_(:A, t)[cells]), mean(mean_(:B, t)[cells])
        mx_a, mx_b = maximum(mean_(:A, t)[cells]), maximum(mean_(:B, t)[cells])
        say(@sprintf("   %-4s mean %9.2f / %9.2f   max of the mean map %9.2f / %9.2f", t, a, b, mx_a, mx_b))
    end
end
ΔNO₃ = mean_(:B, :NO₃) .- mean_(:A, :NO₃)
say(@sprintf("\nsurface NO₃ increase by the loads (mean map): cells > 1 mmol/m³: %d, > 5: %d, > 20: %d; largest %.1f", count(>(1), filter(isfinite, ΔNO₃)), count(>(5), filter(isfinite, ΔNO₃)),
             count(>(20), filter(isfinite, ΔNO₃)), maximum(filter(isfinite, ΔNO₃))))
sel = [(i, j) for i in 1:Nx, j in 1:Ny if isfinite(ΔNO₃[i, j]) && ΔNO₃[i, j] > 2]
for (name, (λ, φ)) in MOUTHS
    near = [(i, j) for (i, j) in sel if km(λm(i), φm(j), λ, φ) < 300 && km(λm(i), φm(j), λ, φ) < min(km(λm(i), φm(j), MOUTHS[1].second...), km(λm(i), φm(j), MOUTHS[2].second...)) + 1e-6]
    isempty(near) && continue
    d = [km(λm(i), φm(j), λ, φ) for (i, j) in near]
    say(@sprintf("  %s: cells with ΔNO₃ > 2 nearest to this mouth: %d, farthest %.0f km (at %.2f°E %.2f°N), median depth %.0f m", name, length(near), maximum(d),
                 λm(near[argmax(d)][1]), φm(near[argmax(d)][2]), median([depth[i, j] for (i, j) in near])))
end
say("\nsurface chlorophyll in the plume regions against the satellite (mg/m³; statistics of log10), all matched cell-days")
for reg in ("Hudson plume", "Delaware plume"), tag in (:A, :B)
    o, m = pairs_[(tag, reg)]; ok = (o .> 0) .& (m .> 0); o, m = o[ok], m[ok]
    length(o) < 5 && continue
    say(@sprintf("  %-15s run %s: n=%6d  satellite %.2f  model %.2f  bias(log10) %+.2f  rms(log10) %.2f  r(log10) %+.2f", reg, tag, length(o), mean(o), mean(m),
                 mean(log10.(m) .- log10.(o)), sqrt(mean((log10.(m) .- log10.(o)) .^ 2)), cor(log10.(o), log10.(m))))
end

# ---------------- figures ----------------
λ_axis = [λm(i) for i in 1:Nx]; φ_axis = [φm(j) for j in 1:Ny]
fig = Figure(size = (1700, 520))
for (c, (t, label, range)) in enumerate(((:NO₃, "surface NO₃ (mmol/m³)", (-10, 10)), (:Alk, "surface Alk (meq/m³)", (-300, 300)), (:chl, "surface chlorophyll (mg/m³)", (-1, 1))))
    ax = Axis(fig[1, c], title = "B − A, $label", xlabel = "°E", ylabel = c == 1 ? "°N" : "", aspect = DataAspect())
    h = heatmap!(ax, λ_axis, φ_axis, mean_(:B, t) .- mean_(:A, t); colormap = :balance, colorrange = range)
    Colorbar(fig[2, c], h, vertical = false, flipaxis = false)
    for (name, (λ, φ)) in MOUTHS; scatter!(ax, [λ], [φ]; color = :white, strokecolor = :black, strokewidth = 1.5, markersize = 10); end
    xlims!(ax, -76, -69); ylims!(ax, 36.5, 41.5)
end
save(joinpath(OUT, "river_difference_maps.png"), fig)
days_axis = Dates.value.(Date.(frame_dates) .- Date(START))
fig = Figure(size = (1500, 800))
for (r, reg) in enumerate(("Hudson plume", "Delaware plume")), (c, (t, label)) in enumerate(((:NO₃, "NO₃ (mmol/m³)"), (:chl, "chlorophyll (mg/m³)"), (:S, "salinity")))
    ax = Axis(fig[r, c], title = "$reg: surface $label", xlabel = r == 2 ? "days from $(Date(START))" : "")
    lines!(ax, days_axis, series[(:A, reg, t)]; label = "A, no loads", color = :black); lines!(ax, days_axis, series[(:B, reg, t)]; label = "B, loads", color = :tomato)
    r == 1 && c == 1 && axislegend(ax; position = :lt)
end
save(joinpath(OUT, "river_plume_series.png"), fig)
say("\nfigures and summary in $OUT")
close(io)
