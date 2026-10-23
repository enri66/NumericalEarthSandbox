# Comparison figures and an animation of surface chlorophyll for several BGC runs against the satellite (Copernicus L3 daily
# multi-sensor, 4 km) and the OOI Pioneer sensors. The model's chlorophyll is 1.31 mg Chl/mmol N × the surface P of the run's
# de-tided daily output.
#   julia -t 4 --project=<env with CairoMakie> bgc/chl_comparison.jl
# Environment: RUNS = "label=path,label=path,..." (the last one is the one animated), MAB_START_DATE, MAB_CELLS_PER_DEGREE,
# BGC_SAT, BGC_OOI, BGC_OUT, FRAMERATE. Writes chl_maps_runs.png, chl_series_runs.png, pioneer_runs.png, surface_chl.mp4.
using Printf, Statistics, Dates
using CairoMakie
using Oceananigans
using Oceananigans.Grids: znodes
using NumericalEarth.DataWrangling: jldopen
const NCD = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "NCDatasets")]

const RUN_SPEC = get(ENV, "RUNS", "no river loads=/t0/workdir/enrique/runs/bgc/bgc_feb19_cmb,rivers=/t0/workdir/enrique/runs/bgc/bgc_feb19_rc4,rivers + scaling=/t0/workdir/enrique/runs/bgc/bgc_feb19_sc")
const RUNS = [(String(strip(first(split(p, "=", limit = 2)))), String(strip(last(split(p, "=", limit = 2))))) for p in split(RUN_SPEC, ",")]
const START = DateTime(get(ENV, "MAB_START_DATE", "2019-02-01"))
const CPD   = parse(Int, get(ENV, "MAB_CELLS_PER_DEGREE", "12"))
const SAT   = get(ENV, "BGC_SAT", "/t0/workdir/enrique/Data/ocean_colour/chl_l3_multi_4km_mab_2019.nc")
const OOI   = get(ENV, "BGC_OOI", "/t0/workdir/enrique/Data/OOI/pioneer_bgc_2019")
const OUT   = get(ENV, "BGC_OUT", "/t0/workdir/enrique/runs/analysis/chl_comparison")
const FRAMERATE = parse(Int, get(ENV, "FRAMERATE", "6"))
const H, Nz = 7, parse(Int, get(ENV, "MAB_NZ", "50"))
const CHL_PER_N = 1.31
const λ_west, φ_south = -76 + 2 / 12, 34 + 2 / 12
mkpath(OUT)
colors = [:gray45, :steelblue, :tomato, :darkgreen, :purple][1:length(RUNS)]
say(args...) = (println(args...); flush(stdout))

# ---------------- layout and depth ----------------
function exponential_scale(N, H, Δtop)
    top(h) = H * expm1(H / N / h) / expm1(H / h)
    lo, hi = H / 50, 100H
    for _ in 1:200; mid = (lo + hi) / 2; top(mid) > Δtop ? (hi = mid) : (lo = mid); end
    return (lo + hi) / 2
end
Δz_top, Z_BOTTOM = parse(Float64, get(ENV, "MAB_DZ_TOP", "2")), parse(Float64, get(ENV, "MAB_ZBOTTOM", "5500"))
zgrid = RectilinearGrid(size = (1, 1, Nz), x = (0, 1), y = (0, 1), z = ExponentialDiscretization(Nz, -Z_BOTTOM, 0; scale = exponential_scale(Nz, Z_BOTTOM, Δz_top)))
const zc = collect(znodes(zgrid, Center())); const Δz = diff(collect(znodes(zgrid, Face())))
bottom_files(tag) = sort(filter(f -> occursin(Regex("^" * basename(tag) * "_bottom_rank\\d+\\.jld2\$"), basename(f)), readdir(dirname(tag); join = true)))
bf = first(filter(!isempty, [bottom_files(last(r)) for r in reverse(RUNS)]))
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

iters(f, name) = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/$name"]))))
strip3(A) = A[H+1:end-H, H+1:end-H, H+1:end-H]
function field3(tag, name, it)
    A = zeros(Float32, Nx, Ny, Nz)
    for (r, p) in enumerate(parts)
        jldopen("$(tag)_bgc_volume_daily_rank$(r - 1).jld2") do f
            B = strip3(f["timeseries/$name/$it"]); A[p.I+1:p.I+size(B, 1), p.J+1:p.J+size(B, 2), :] .= B
        end
    end
    return A
end
frame_index(tag) = jldopen("$(tag)_bgc_volume_daily_rank0.jld2") do f
    Dict(round(Int, f["timeseries/t/$it"]) => it for it in iters(f, "P"))
end
index = [frame_index(last(r)) for r in RUNS]
common = sort(collect(reduce(intersect, [Set(keys(ix)) for ix in index])))
dates = [START + Second(t) for t in common]; Nf = length(common)
say("runs: $(join([first(r) for r in RUNS], "; ")); $Nf shared daily frames, $(Date(first(dates))) → $(Date(last(dates)))")

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

# ---------------- moorings and observations ----------------
const MOORINGS = ("CNSM" => ("cp01cnsm", 147.0), "ISSM" => ("cp03issm", 92.0))
function read_ooi(path, column)
    isfile(path) || return nothing
    lines = readlines(path)
    (length(lines) > 2 && startswith(lines[1], "time,")) || return nothing
    col = Dict(h => i for (i, h) in enumerate(split(lines[1], ",")))
    haskey(col, column) || return nothing
    qc = col[column * "_qc_agg"]; out = Tuple{DateTime, Float64}[]; lat = Float64[]; lon = Float64[]
    for line in lines[3:end]
        v = split(line, ","); v[qc] in ("1", "2") || continue
        x = tryparse(Float64, v[col[column]]); isnothing(x) && continue
        push!(out, (DateTime(v[1][1:19]), x)); push!(lat, parse(Float64, v[col["latitude"]])); push!(lon, parse(Float64, v[col["longitude"]]))
    end
    return isempty(out) ? nothing : (; data = out, lat = median(lat), lon = median(lon))
end
obs_at(data, T; half = Hour(36), minimum_samples = 24) = (v = [x for (t, x) in data if abs(t - T) ≤ half]; length(v) ≥ minimum_samples ? mean(v) : NaN)
mooring = Dict{String, NamedTuple}()
for (name, (prefix, nominal)) in MOORINGS
    r = read_ooi(joinpath(OOI, "$(prefix)_flort.csv"), "mass_concentration_of_chlorophyll_a_in_sea_water"); isnothing(r) && continue
    cand = [(i, j) for i in 1:Nx, j in 1:Ny if wet[i, j] && abs(λm(i) - r.lon) * cosd(r.lat) < 0.13 && abs(φm(j) - r.lat) < 0.13]
    best = argmin([abs(depth[i, j] - nominal) / nominal + 0.5 * hypot((λm(i) - r.lon) * cosd(r.lat), φm(j) - r.lat) / 0.13 for (i, j) in cand])
    chl_obs = [obs_at([(t, x) for (t, x) in r.data if 0 ≤ x ≤ 30], T) for T in dates]
    no3 = read_ooi(joinpath(OOI, "$(prefix)_nutnr_all.csv"), "mole_concentration_of_nitrate_in_sea_water_suna")
    no3_obs = isnothing(no3) ? fill(NaN, Nf) : [obs_at([(t, x) for (t, x) in no3.data if 0 ≤ x ≤ 40], T) for T in dates]
    mooring[name] = (; cell = cand[best], lon = r.lon, lat = r.lat, chl_obs, no3_obs)
end
function at_depth(col, k0, z0)
    z0 = -z0
    for k in k0:Nz-1; zc[k] ≤ z0 ≤ zc[k+1] && (w = (z0 - zc[k]) / (zc[k+1] - zc[k]); return (1 - w) * col[k] + w * col[k+1]); end
    return col[Nz]
end

# ---------------- one pass over the frames ----------------
chl = [fill(NaN32, Nx, Ny, Nf) for _ in RUNS]                    # model surface chlorophyll per run
sat = fill(NaN32, Nx, Ny, Nf)                                     # satellite on the model cells
mod_chl7 = [Dict(n => fill(NaN, Nf) for n in keys(mooring)) for _ in RUNS]; mod_no3 = [Dict(n => fill(NaN, Nf) for n in keys(mooring)) for _ in RUNS]
sat_at_mooring = Dict(n => fill(NaN, Nf) for n in keys(mooring))
for (n, t) in enumerate(common)
    satn = get(sat_day, Date(dates[n]), 0)
    if satn > 0
        for j in 1:Ny, i in 1:Nx; wet[i, j] && (sat[i, j, n] = sat_cell(satn, i, j)); end
        for (name, m) in mooring; sat_at_mooring[name][n] = sat_cell(satn, m.cell...); end
    end
    for (r, run) in enumerate(RUNS)
        P = field3(last(run), "P", index[r][t]); NO₃ = isempty(mooring) ? nothing : field3(last(run), "NO₃", index[r][t])
        chl[r][:, :, n] .= ifelse.(wet, CHL_PER_N .* P[:, :, Nz], NaN32)
        for (name, m) in mooring
            i, j = m.cell
            mod_chl7[r][name][n] = CHL_PER_N * at_depth(P[i, j, :], kb[i, j], 7.0); mod_no3[r][name][n] = at_depth(NO₃[i, j, :], kb[i, j], 7.0)
        end
    end
    n % 10 == 0 && say("  frame $n of $Nf")
end
sat[.!isfinite.(sat)] .= NaN32

# ---------------- statistics ----------------
valid = isfinite.(sat)
say("\nsurface chlorophyll against the satellite over $Nf days (matched cell-days; log10): label, region, mean model, mean satellite, bias, rms, r")
for (r, run) in enumerate(RUNS), (reg, m) in pairs(regions)
    sel = [valid[i, j, n] && m[i, j] && chl[r][i, j, n] > 0 for i in 1:Nx, j in 1:Ny, n in 1:Nf]
    o, md = Float64.(sat[sel]), Float64.(chl[r][sel]); lo, lm = log10.(o), log10.(md)
    say(@sprintf("  %-18s %-6s model %.2f  satellite %.2f  bias %+.2f  rms %.2f  r %+.2f", first(run), reg, mean(md), mean(o), mean(lm .- lo), sqrt(mean((lm .- lo) .^ 2)), cor(lo, lm)))
end

# ---------------- figures ----------------
λ_axis = [λm(i) for i in 1:Nx]; φ_axis = [φm(j) for j in 1:Ny]
crange = (log10(0.1), log10(8.0))
function matched_mean(A)
    S = zeros(Nx, Ny); C = zeros(Int, Nx, Ny)
    for n in 1:Nf, j in 1:Ny, i in 1:Nx; (valid[i, j, n] && isfinite(A[i, j, n])) && (S[i, j] += A[i, j, n]; C[i, j] += 1); end
    return [C[i, j] ≥ 5 ? S[i, j] / C[i, j] : NaN for i in 1:Nx, j in 1:Ny]
end
msat = matched_mean(sat); mmod = [matched_mean(c) for c in chl]
fig = Figure(size = (400 * (length(RUNS) + 1) + 120, 900))
ax = Axis(fig[1, 1], title = "satellite (days with data)", ylabel = "°N", aspect = DataAspect())
h = heatmap!(ax, λ_axis, φ_axis, log10.(msat); colormap = :viridis, colorrange = crange)
for (r, run) in enumerate(RUNS)
    ax = Axis(fig[1, r + 1], title = first(run) * " (same days and pixels)", aspect = DataAspect()); heatmap!(ax, λ_axis, φ_axis, log10.(mmod[r]); colormap = :viridis, colorrange = crange)
    ax = Axis(fig[2, r + 1], title = "log10 model / satellite", xlabel = "°E", aspect = DataAspect())
    h2 = heatmap!(ax, λ_axis, φ_axis, log10.(mmod[r]) .- log10.(msat); colormap = :balance, colorrange = (-1, 1))
    r == length(RUNS) && Colorbar(fig[2, length(RUNS) + 2], h2, label = "log10 ratio")
end
Colorbar(fig[1, length(RUNS) + 2], h, label = "log10 chlorophyll (mg/m³)")
save(joinpath(OUT, "chl_maps_runs.png"), fig)

days = [t / 86400 for t in common]
fig = Figure(size = (1500, 420))
for (c, (reg, m)) in enumerate(pairs(regions))
    ax = Axis(fig[1, c], title = "$reg: region-mean surface chlorophyll (satellite days)", xlabel = "days from $(Date(START))", ylabel = c == 1 ? "mg/m³" : "", yscale = log10)
    series(A) = [(s = [A[i, j, n] for i in 1:Nx, j in 1:Ny if m[i, j] && valid[i, j, n] && isfinite(A[i, j, n])]; length(s) ≥ 20 ? mean(s) : NaN) for n in 1:Nf]
    lines!(ax, days, series(sat); color = :black, linewidth = 2.5, label = "satellite")
    for (r, run) in enumerate(RUNS); lines!(ax, days, series(chl[r]); color = colors[r], label = first(run)); end
    c == 1 && axislegend(ax; position = :lt)
end
save(joinpath(OUT, "chl_series_runs.png"), fig)

if !isempty(mooring)
    fig = Figure(size = (1300, 800))
    for (c, name) in enumerate(sort(collect(keys(mooring))))
        m = mooring[name]
        ax = Axis(fig[1, c], title = "$name: chlorophyll at 7 m (mg/m³)", ylabel = c == 1 ? "mg/m³" : "")
        lines!(ax, days, m.chl_obs; color = :black, linewidth = 2.5, label = "fluorometer"); scatter!(ax, days, sat_at_mooring[name]; color = :steelblue, markersize = 6, label = "satellite")
        for (r, run) in enumerate(RUNS); lines!(ax, days, mod_chl7[r][name]; color = colors[r], label = first(run)); end
        c == 1 && axislegend(ax; position = :lt)
        ax = Axis(fig[2, c], title = "$name: nitrate at 7 m (mmol/m³)", xlabel = "days from $(Date(START))", ylabel = c == 1 ? "mmol/m³" : "")
        lines!(ax, days, m.no3_obs; color = :black, linewidth = 2.5, label = "SUNA")
        for (r, run) in enumerate(RUNS); lines!(ax, days, mod_no3[r][name]; color = colors[r]); end
    end
    save(joinpath(OUT, "pioneer_runs.png"), fig)
end

# ---------------- animation of the last run ----------------
n_ = Observable(1)
anim_title = @lift(Dates.format(dates[$n_], "yyyy-mm-dd"))
sat_frame = @lift(log10.(sat[:, :, $n_])); mod_frame = @lift(log10.(chl[end][:, :, $n_]))
ratio_frame = @lift(log10.(chl[end][:, :, $n_]) .- log10.(sat[:, :, $n_]))
fig = Figure(size = (1800, 620))
Label(fig[0, 1:5], @lift("surface chlorophyll, " * first(RUNS[end]) * ", " * $anim_title), fontsize = 22)
a1 = Axis(fig[1, 1], title = "satellite (white = no data)", xlabel = "°E", ylabel = "°N", aspect = DataAspect())
hh = heatmap!(a1, λ_axis, φ_axis, sat_frame; colormap = :viridis, colorrange = crange)
a2 = Axis(fig[1, 2], title = "model 1.31 P, surface", xlabel = "°E", aspect = DataAspect()); heatmap!(a2, λ_axis, φ_axis, mod_frame; colormap = :viridis, colorrange = crange)
Colorbar(fig[1, 3], hh, label = "log10 chlorophyll (mg/m³)")
a3 = Axis(fig[1, 4], title = "log10 model / satellite", xlabel = "°E", aspect = DataAspect())
hr = heatmap!(a3, λ_axis, φ_axis, ratio_frame; colormap = :balance, colorrange = (-1, 1)); Colorbar(fig[1, 5], hr, label = "log10 ratio")
for (name, m) in mooring, ax in (a1, a2, a3); scatter!(ax, [m.lon], [m.lat]; color = :white, strokecolor = :black, strokewidth = 1.5, markersize = 9); end
record(fig, joinpath(OUT, "surface_chl.mp4"), 1:Nf; framerate = FRAMERATE) do n
    n_[] = n
end
say("\nfigures and animation in $OUT")
