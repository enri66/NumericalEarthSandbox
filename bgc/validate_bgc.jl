# Validation of a BGC run (bgc/mab_bgc.jl) against observations:
#   - surface chlorophyll against satellite ocean colour (Copernicus Marine multi-sensor L3, daily, 4 km), and
#   - NO₃, O₂ and chlorophyll at the OOI Pioneer surface moorings (nitrate, oxygen and fluorometer on the near-surface
#     frame at 7 m; oxygen on the seafloor node).
# The model's chlorophyll is its phytoplankton P (mmol N/m³) times OceanBioME's 1.31 mg Chl/mmol N. It reads the run's
# de-tided daily volume output (<TAG>_bgc_volume_daily_rank<r>.jld2) and the bottom-cell file (<TAG>_bottom_rank<r>.jld2)
# for the layout of the ranks and the depth. The observations are compared with the daily model value, taken as the
# mean of the observations within ±36 h of it (the model output is low-passed with a 40 h cutoff).
#
#   julia -t 4 --project=<env with CairoMakie> bgc/validate_bgc.jl
# Environment: BGC_TAG (run path prefix), MAB_START_DATE, MAB_CELLS_PER_DEGREE, BGC_SAT (satellite NetCDF), BGC_OOI
# (directory of the Pioneer CSVs, see bgc/download_validation_data.sh), BGC_OUT (figures and summary.txt).
using Printf, Statistics, Dates
using CairoMakie
using Oceananigans
using Oceananigans.Grids: znodes
using NumericalEarth.DataWrangling: jldopen
const NCD = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "NCDatasets")]

const TAG   = get(ENV, "BGC_TAG", "/t0/workdir/enrique/runs/bgc/bgc_feb19_cmb")
const START = DateTime(get(ENV, "MAB_START_DATE", "2019-02-01"))
const CPD   = parse(Int, get(ENV, "MAB_CELLS_PER_DEGREE", "12"))
const SAT   = get(ENV, "BGC_SAT", "/t0/workdir/enrique/Data/ocean_colour/chl_l3_multi_4km_mab_2019.nc")
const OOI   = get(ENV, "BGC_OOI", "/t0/workdir/enrique/Data/OOI/pioneer_bgc_2019")
const OUT   = get(ENV, "BGC_OUT", "/t0/workdir/enrique/runs/analysis/bgc_" * basename(TAG))
const H, Nz = 7, parse(Int, get(ENV, "MAB_NZ", "50"))
const CHL_PER_N = 1.31                                   # mg Chl per mmol N (OceanBioME's PhytoZoo)
const λ_west, φ_south = -76 + 2 / 12, 34 + 2 / 12        # script 05's box: GLORYS box edge plus two GLORYS cells
mkpath(OUT)
summary_io = open(joinpath(OUT, "summary.txt"), "w")
function say(args...)
    println(args...); flush(stdout)
    println(summary_io, args...); flush(summary_io)
end

# ---------------- the model's layout, depth and vertical grid ----------------
function exponential_scale(N, H, Δtop)
    top(h) = H * expm1(H / N / h) / expm1(H / h)
    lo, hi = H / 50, 100H
    for _ in 1:200; mid = (lo + hi) / 2; top(mid) > Δtop ? (hi = mid) : (lo = mid); end
    return (lo + hi) / 2
end
Δz_top, Z_BOTTOM = parse(Float64, get(ENV, "MAB_DZ_TOP", "2")), parse(Float64, get(ENV, "MAB_ZBOTTOM", "5500"))
zgrid = RectilinearGrid(size = (1, 1, Nz), x = (0, 1), y = (0, 1),
                        z = ExponentialDiscretization(Nz, -Z_BOTTOM, 0; scale = exponential_scale(Nz, Z_BOTTOM, Δz_top)))
const zc = collect(znodes(zgrid, Center())); const Δz = diff(collect(znodes(zgrid, Face())))   # bottom first

bottom_files = sort(filter(f -> occursin(Regex("^" * basename(TAG) * "_bottom_rank\\d+\\.jld2\$"), basename(f)), readdir(dirname(TAG); join = true)))
NR = length(bottom_files)
NR > 0 || error("no $(basename(TAG))_bottom_rank*.jld2 in $(dirname(TAG)): the run was started before mab_bgc.jl wrote them")
parts = map(bottom_files) do f
    jldopen(f) do file
        (I = file["I_OFF"], J = file["J_OFF"], kb = file["k_bottom"], dzb = file["Δz_bottom"])
    end
end
const Nx = maximum(p.I + size(p.kb, 1) for p in parts); const Ny = maximum(p.J + size(p.kb, 2) for p in parts)
kb = zeros(Int, Nx, Ny); depth = zeros(Nx, Ny)
for p in parts, j in axes(p.kb, 2), i in axes(p.kb, 1)
    k = p.kb[i, j]
    kb[p.I + i, p.J + j] = k
    depth[p.I + i, p.J + j] = k > 0 ? sum(Δz[k+1:end]) + p.dzb[i, j] : 0.0
end
wet = kb .> 0
λm(i) = λ_west + (i - 0.5) / CPD
φm(j) = φ_south + (j - 0.5) / CPD
say("model: $(Nx)×$(Ny)×$(Nz) on $NR ranks, $(count(wet)) wet columns")

iters(f, name) = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/$name"]))))
volume_files = ["$(TAG)_bgc_volume_daily_rank$(r).jld2" for r in 0:NR-1]
frame_times = jldopen(f -> [f["timeseries/t/$it"] for it in iters(f, "P")], volume_files[1])
function frame(name, n; stem = "bgc_volume_daily")
    A = zeros(Float32, Nx, Ny, Nz)
    for (p, file) in zip(parts, ["$(TAG)_$(stem)_rank$(r).jld2" for r in 0:NR-1])
        jldopen(file) do f
            it = iters(f, name)[n]
            B = f["timeseries/$name/$it"][H+1:end-H, H+1:end-H, H+1:end-H]
            A[p.I+1:p.I+size(B, 1), p.J+1:p.J+size(B, 2), :] .= B
        end
    end
    return A
end
frame_dates = [START + Millisecond(round(Int, 1000t)) for t in frame_times]
say("model frames: $(length(frame_dates)) daily, $(Dates.format(first(frame_dates), "yyyy-mm-dd")) → $(Dates.format(last(frame_dates), "yyyy-mm-dd"))")

regions = (shelf = (depth .> 0) .& (depth .< 100), slope = (depth .≥ 100) .& (depth .< 1000), deep = depth .≥ 1000)
say(@sprintf("regions: shelf (< 100 m) %d columns, upper slope (100-1000 m) %d, deep (> 1000 m) %d", count(regions.shelf), count(regions.slope), count(regions.deep)))

# oxygen solubility (µmol/kg) at temperature T (°C) and salinity S, Garcia & Gordon (1992)
function o2_saturation(T, S)
    Ts = log((298.15 - T) / (273.15 + T))
    lnC = 5.80871 + 3.20291Ts + 4.17887Ts^2 + 5.10006Ts^3 - 0.0986643Ts^4 + 3.80369Ts^5 +
          S * (-0.00701577 - 0.00770028Ts - 0.0113864Ts^2 - 0.00951519Ts^3) - 2.75915e-7 * S^2
    return exp(lnC)
end
const ρ₀ = 1025.0

# ---------------- the model at the Pioneer moorings ----------------
# name => (CSV prefix, nominal water depth (m), instrument depth of the near-surface frame (m))
# (the offshore mooring, OSSM, has data only until 2019-02-01, before the first daily model frame)
const MOORINGS = ("CNSM" => ("cp01cnsm", 147.0, 7.0), "ISSM" => ("cp03issm", 92.0, 7.0))
function read_ooi(path, column)
    isfile(path) || return nothing
    lines = readlines(path)
    (length(lines) > 2 && startswith(lines[1], "time,")) || return nothing
    col = Dict(h => i for (i, h) in enumerate(split(lines[1], ",")))
    haskey(col, column) || return nothing
    qc = col[column * "_qc_agg"]
    out = Tuple{DateTime, Float64}[]; lat = Float64[]; lon = Float64[]
    for line in lines[3:end]
        v = split(line, ",")
        v[qc] in ("1", "2") || continue
        x = tryparse(Float64, v[col[column]]); isnothing(x) && continue
        push!(out, (DateTime(v[1][1:19]), x)); push!(lat, parse(Float64, v[col["latitude"]])); push!(lon, parse(Float64, v[col["longitude"]]))
    end
    isempty(out) && return nothing
    return (; data = out, lat = median(lat), lon = median(lon))
end
mooring_cell = Dict{String, Tuple{Int, Int}}()
mooring_pos  = Dict{String, Tuple{Float64, Float64}}()
for (name, (prefix, nominal, _)) in MOORINGS
    r = nothing
    for (suffix, column) in (("flort", "mass_concentration_of_chlorophyll_a_in_sea_water"), ("dosta_nsif", "mole_concentration_of_dissolved_molecular_oxygen_in_sea_water"))
        isnothing(r) && (r = read_ooi(joinpath(OOI, "$(prefix)_$suffix.csv"), column))
    end
    isnothing(r) && continue
    # the wet cell within 0.13° whose depth is closest to the mooring's, ties to the nearer one
    cand = [(i, j) for i in 1:Nx, j in 1:Ny if wet[i, j] && abs(λm(i) - r.lon) * cosd(r.lat) < 0.13 && abs(φm(j) - r.lat) < 0.13]
    best = argmin([abs(depth[i, j] - nominal) / nominal + 0.5 * hypot((λm(i) - r.lon) * cosd(r.lat), φm(j) - r.lat) / 0.13 for (i, j) in cand])
    mooring_cell[name] = cand[best]; mooring_pos[name] = (r.lon, r.lat)
    i, j = cand[best]
    say(@sprintf("%s: mooring at %.3f°E %.3f°N, nominal depth %.0f m → model cell (%d, %d) at %.3f°E %.3f°N, model depth %.0f m",
                 name, r.lon, r.lat, nominal, i, j, λm(i), φm(j), depth[i, j]))
end

# value at depth z0 (m, positive down) of a column, linear between cell centres
function at_depth(col, k0, z0)
    z0 = -z0
    k0 > 0 || return NaN
    for k in k0:Nz-1
        if zc[k] ≤ z0 ≤ zc[k+1]
            w = (z0 - zc[k]) / (zc[k+1] - zc[k]); return (1 - w) * col[k] + w * col[k+1]
        end
    end
    return col[Nz]
end

# ---------------- satellite chlorophyll ----------------
ds = NCD.Dataset(SAT)
λs = Float64.(ds["longitude"][:]); φs = Float64.(ds["latitude"][:]); sat_times = ds["time"][:]
raw = ds["CHL"][:, :, :]
chl_sat = Array{Float32}(undef, size(raw)); for k in eachindex(raw); chl_sat[k] = ismissing(raw[k]) ? NaN32 : Float32(raw[k]); end
close(ds)
sat_day = Dict(Date(t) => n for (n, t) in enumerate(sat_times))
lon_idx = [findall(abs.(λs .- λm(i)) .< 0.5 / CPD) for i in 1:Nx]
lat_idx = [findall(abs.(φs .- φm(j)) .< 0.5 / CPD) for j in 1:Ny]
function sat_cell(n, i, j)
    v = [x for x in chl_sat[lon_idx[i], lat_idx[j], n] if isfinite(x) && x > 0]
    return length(v) ≥ 2 ? Float64(mean(v)) : NaN
end
say(@sprintf("satellite: %d daily L3 maps, %.1f%% of pixels valid on average", length(sat_times), 100 * count(isfinite, chl_sat) / length(chl_sat)))

# ---------------- one pass over the model frames ----------------
Ndays = length(frame_dates)
sat_sum = zeros(Nx, Ny); mod_sum = zeros(Nx, Ny); cnt = zeros(Int, Nx, Ny)
matched = Dict(name => (Float64[], Float64[]) for name in keys(regions))
series = Dict(name => (fill(NaN, Ndays), fill(NaN, Ndays)) for name in keys(regions))       # region-mean satellite, model
moor_model = Dict(name => Dict(v => fill(NaN, Ndays) for v in (:chl, :NO₃, :O₂, :O₂_bottom, :O₂_sat)) for name in keys(mooring_cell))
moor_sat = Dict(name => fill(NaN, Ndays) for name in keys(mooring_cell))        # satellite chlorophyll at the mooring's cell
for n in 1:Ndays
    P, NO₃, O₂ = frame("P", n), frame("NO₃", n), frame("O₂", n)
    Tm, Sm = frame("T", n; stem = "volume_daily"), frame("S", n; stem = "volume_daily")
    chl_model = CHL_PER_N .* P[:, :, Nz]
    d = Date(frame_dates[n])
    if haskey(sat_day, d)
        s = [wet[i, j] ? sat_cell(sat_day[d], i, j) : NaN for i in 1:Nx, j in 1:Ny]
        ok = isfinite.(s) .& (chl_model .> 0)
        sat_sum[ok] .+= s[ok]; mod_sum[ok] .+= chl_model[ok]; cnt[ok] .+= 1
        for (name, mask) in pairs(regions)
            sel = ok .& mask
            append!(matched[name][1], s[sel]); append!(matched[name][2], chl_model[sel])
            count(sel) ≥ 20 && (series[name][1][n] = mean(s[sel]); series[name][2][n] = mean(chl_model[sel]))
        end
    end
    for (name, (i, j)) in mooring_cell
        moor_model[name][:chl][n]  = CHL_PER_N * at_depth(P[i, j, :], kb[i, j], 7.0)
        moor_model[name][:NO₃][n]  = at_depth(NO₃[i, j, :], kb[i, j], 7.0)
        moor_model[name][:O₂][n]   = at_depth(O₂[i, j, :], kb[i, j], 7.0)
        moor_model[name][:O₂_bottom][n] = O₂[i, j, kb[i, j]]
        T7, S7 = at_depth(Tm[i, j, :], kb[i, j], 7.0), at_depth(Sm[i, j, :], kb[i, j], 7.0)
        moor_model[name][:O₂_sat][n] = 100 * moor_model[name][:O₂][n] / (o2_saturation(T7, S7) * ρ₀ / 1000)
        haskey(sat_day, d) && (moor_sat[name][n] = sat_cell(sat_day[d], i, j))
    end
    n % 10 == 0 && say("  frame $n of $Ndays")
end

function summarize(label, o, m; logspace = false)
    ok = isfinite.(o) .& isfinite.(m)
    logspace && (ok .&= (o .> 0) .& (m .> 0))
    o, m = o[ok], m[ok]
    length(o) < 5 && return say(@sprintf("  %-28s only %d pairs", label, length(o)))
    if logspace
        lo, lm = log10.(o), log10.(m)
        say(@sprintf("  %-28s n=%6d  mean obs %7.3f  mean model %7.3f  bias(log10) %+.2f  rms(log10) %.2f  r(log10) %+.2f",
                     label, length(o), mean(o), mean(m), mean(lm .- lo), sqrt(mean((lm .- lo) .^ 2)), cor(lo, lm)))
    else
        say(@sprintf("  %-28s n=%6d  mean obs %7.2f  mean model %7.2f  bias %+7.2f  rms %6.2f  rms without bias %6.2f  std model/obs %.2f  r %+.2f",
                     label, length(o), mean(o), mean(m), mean(m .- o), sqrt(mean((m .- o) .^ 2)),
                     sqrt(mean(((m .- mean(m)) .- (o .- mean(o))) .^ 2)), std(m) / std(o), cor(o, m)))
    end
end

say("\nsurface chlorophyll, model against satellite (mg/m³; all matched cell-days; statistics of log10)")
for name in keys(regions)
    summarize(string(name), matched[name]..., logspace = true)
end
summarize("whole domain", reduce(vcat, first.(values(matched))), reduce(vcat, last.(values(matched))), logspace = true)
say("\nregion-mean time series (days with ≥ 20 matched cells)")
for name in keys(regions)
    summarize(string(name, " daily mean"), series[name]..., logspace = true)
end

# by water depth: where the shelf deficit is (mean over the matched cell-days of each depth bin; ratio = mean model / mean satellite)
say("\nsurface chlorophyll by water depth (mg/m³, matched cell-days)")
for (lo, hi) in ((0, 20), (20, 40), (40, 60), (60, 80), (80, 100), (100, 200), (200, 1000), (1000, Inf))
    sel = wet .& (depth .> lo) .& (depth .≤ hi) .& (cnt .> 0)
    count(sel) == 0 && continue
    n = sum(cnt[sel])
    say(@sprintf("  %5.0f-%-5.0f m  columns %5d  cell-days %6d  satellite %6.3f  model %6.3f  model/satellite %.2f",
                 lo, hi, count(sel), n, sum(sat_sum[sel]) / n, sum(mod_sum[sel]) / n, sum(mod_sum[sel]) / sum(sat_sum[sel])))
end

# ---------------- Pioneer comparison ----------------
function obs_at(data, T; half = Hour(36), minimum_samples = 24)
    v = [x for (t, x) in data if abs(t - T) ≤ half]
    return length(v) ≥ minimum_samples ? mean(v) : NaN
end
limits = (chl = (0.0, 30.0), NO₃ = (0.0, 40.0), O₂ = (50.0, 500.0), O₂_bottom = (50.0, 500.0))
variables = ((:chl, "flort", "mass_concentration_of_chlorophyll_a_in_sea_water"),
             (:NO₃, "nutnr_all", "mole_concentration_of_nitrate_in_sea_water_suna"),
             (:O₂, "dosta_nsif", "mole_concentration_of_dissolved_molecular_oxygen_in_sea_water"),
             (:O₂_bottom, "dosta_mfn", "mole_concentration_of_dissolved_molecular_oxygen_in_sea_water"))
moor_obs = Dict(name => Dict{Symbol, Vector{Float64}}() for name in keys(mooring_cell))
say("\nPioneer moorings, model daily value against the observations within ±36 h (µg/L = mg/m³, µmol/L = mmol/m³)")
for (name, (prefix, _, _)) in MOORINGS, (v, suffix, column) in variables
    haskey(mooring_cell, name) || continue
    r = read_ooi(joinpath(OOI, "$(prefix)_$suffix.csv"), column)
    isnothing(r) && (say("  $name $v: no data"); continue)
    lo, hi = limits[v]
    data = [(t, x) for (t, x) in r.data if lo ≤ x ≤ hi]
    o = [obs_at(data, T) for T in frame_dates]
    moor_obs[name][v] = o
    summarize("$name $v ($(count(isfinite, o)) days)", o, moor_model[name][v])
end

# oxygen saturation of the sensors' own readings (their temperature and salinity), and the satellite at the mooring's cell
say("\noxygen saturation at 7 m (%) and satellite chlorophyll against the fluorometers")
for (name, (prefix, _, _)) in MOORINGS
    haskey(mooring_cell, name) || continue
    lines = readlines(joinpath(OOI, "$(prefix)_dosta_nsif.csv"))
    col = Dict(h => i for (i, h) in enumerate(split(lines[1], ",")))
    o2c, tc, sc, qc = col["mole_concentration_of_dissolved_molecular_oxygen_in_sea_water"], col["sea_water_temperature"], col["sea_water_practical_salinity"], col["mole_concentration_of_dissolved_molecular_oxygen_in_sea_water_qc_agg"]
    samples = Tuple{DateTime, Float64}[]
    for line in lines[3:end]
        v = split(line, ",")
        v[qc] in ("1", "2") || continue
        o, T, S = tryparse(Float64, v[o2c]), tryparse(Float64, v[tc]), tryparse(Float64, v[sc])
        (isnothing(o) || isnothing(T) || isnothing(S) || !(50 < o < 500)) && continue
        push!(samples, (DateTime(v[1][1:19]), 100 * o / (o2_saturation(T, S) * ρ₀ / 1000)))
    end
    sat_obs = [obs_at(samples, T) for T in frame_dates]
    moor_obs[name][:O₂_sat] = sat_obs
    summarize("$name O₂ saturation", sat_obs, moor_model[name][:O₂_sat])
    if haskey(moor_obs[name], :chl)
        summarize("$name satellite vs fluorometer", moor_obs[name][:chl], moor_sat[name], logspace = true)
    end
end

# ---------------- figures ----------------
λ_axis = [λm(i) for i in 1:Nx]; φ_axis = [φm(j) for j in 1:Ny]
mean_map(sum_) = [cnt[i, j] ≥ 5 ? sum_[i, j] / cnt[i, j] : NaN for i in 1:Nx, j in 1:Ny]
msat, mmod = mean_map(sat_sum), mean_map(mod_sum)
fig = Figure(size = (1500, 520))
crange = (log10(0.1), log10(8.0))
ax1 = Axis(fig[1, 1], title = "satellite chlorophyll, matched days mean", xlabel = "°E", ylabel = "°N", aspect = DataAspect())
h = heatmap!(ax1, λ_axis, φ_axis, log10.(msat); colormap = :viridis, colorrange = crange)
ax2 = Axis(fig[1, 2], title = "model 1.31 P at the surface, same days and pixels", xlabel = "°E", aspect = DataAspect())
heatmap!(ax2, λ_axis, φ_axis, log10.(mmod); colormap = :viridis, colorrange = crange)
Colorbar(fig[1, 3], h, label = "log10 chlorophyll (mg/m³)")
ax3 = Axis(fig[1, 4], title = "log10 model / satellite", xlabel = "°E", aspect = DataAspect())
h3 = heatmap!(ax3, λ_axis, φ_axis, log10.(mmod) .- log10.(msat); colormap = :balance, colorrange = (-1.0, 1.0))
Colorbar(fig[1, 5], h3, label = "log10 ratio")
for (name, (λ, φ)) in mooring_pos
    for ax in (ax1, ax2, ax3); scatter!(ax, [λ], [φ]; color = :white, strokecolor = :black, strokewidth = 1.5, markersize = 9); end
    text!(ax3, λ, φ; text = " " * name, fontsize = 11, align = (:left, :bottom))
end
save(joinpath(OUT, "chl_maps.png"), fig)

dates_axis = Dates.value.(Date.(frame_dates) .- Date(START))
fig = Figure(size = (1400, 420))
for (c, name) in enumerate(keys(regions))
    ax = Axis(fig[1, c], title = "$(name): region-mean surface chlorophyll", xlabel = "days from $(Date(START))", ylabel = "mg/m³", yscale = log10)
    lines!(ax, dates_axis, series[name][1]; label = "satellite", color = :black)
    lines!(ax, dates_axis, series[name][2]; label = "model", color = :tomato)
    c == 1 && axislegend(ax; position = :lt)
end
save(joinpath(OUT, "chl_timeseries.png"), fig)

names_ = collect(keys(mooring_cell)); rows = ((:chl, "chlorophyll (mg/m³)"), (:NO₃, "nitrate (mmol/m³)"), (:O₂, "oxygen at 7 m (mmol/m³)"), (:O₂_sat, "oxygen saturation at 7 m (%)"), (:O₂_bottom, "oxygen at the seafloor (mmol/m³)"))
fig = Figure(size = (1200, 1250))
for (c, name) in enumerate(names_), (r, (v, label)) in enumerate(rows)
    ax = Axis(fig[r, c], title = r == 1 ? name : "", ylabel = c == 1 ? label : "", xlabel = r == length(rows) ? "days from $(Date(START))" : "")
    haskey(moor_obs[name], v) && lines!(ax, dates_axis, moor_obs[name][v]; label = "Pioneer", color = :black)
    lines!(ax, dates_axis, moor_model[name][v]; label = "model", color = :tomato)
    v == :chl && scatter!(ax, dates_axis, moor_sat[name]; label = "satellite", color = :steelblue, markersize = 6)
    v == :O₂_sat && hlines!(ax, [100.0]; color = :gray, linestyle = :dash)
    r == 1 && c == 1 && axislegend(ax; position = :lt)
end
save(joinpath(OUT, "pioneer_bgc.png"), fig)
say("\nfigures and summary in $OUT")
close(summary_io)
