# Sea ice of a GSL run (script 05 with MAB_REGION=gsl, MAB_SEA_ICE=true) against observations and GLORYS, winter 2013-14,
# on the model's 1/12° grid. Daily ice area (Σ ℵ A) and extent (Σ A where ℵ ≥ 0.15) in the Gulf of St. Lawrence, on the
# Labrador side north of 51.7°N, in Lake Melville (part of the Labrador side) and on the Newfoundland shelf and Grand
# Banks (the rest), leaving out the sponge band (8 cells from the box edges), as in gsl_ice_vs_glorys.jl:
#   - model: daily mean of the 6-hourly <tag>_sea_ice_rank0.jld2 (h, ℵ from script 05's JLD2Writer);
#   - GLORYS daily siconc, sithick (the region cache NumericalEarth reads);
#   - OSI SAF OSI-450-a v3.0 and NOAA/NSIDC CDR G02202 v6 (passive microwave, 25 km): each model cell takes the nearest
#     valid pixel within 25 km, among pixels about 50 km or more from land (land spill-over puts false ice along the
#     coasts), so cells near the coast and in narrow water (Lake Melville, the estuary) have no value; model and GLORYS
#     are summed over the same cells, and the tables give the fraction of each region's water those cells cover;
#   - Canadian Ice Service weekly East Coast charts (SIGRID-3, NSIDC G02171): total concentration CT and a thickness from
#     the stages of development (WMO mid-range values), each model cell taking the polygon its centre falls in (cells on
#     chart land the nearest chart cell within 2 cells); also the ice volume Σ Cₖ hₖ A;
#   - AWI CryoSat-2/SMOS v300 thickness (7-day windows, 12.5 km), the window centred on each chart date, for the volume.
# Observations live under <obs dir>/{osisaf_osi450a_v3,nsidc_g02202_v6,cis_g02171,cs2smos_v300} (see
# ~/Data/sea_ice_obs/README.md on antares). Writes <tag>_ice_obs_vs_model.csv (daily) and _cis.csv (chart dates), and
# <tag>_ice_obs_maps.jld2 with the daily fields on the model grid.
#
#   julia --project=<env> sea_ice_obs_vs_model.jl <run dir> <tag> [obs dir] [GLORYS dir] [bathymetry file]
# ICE_OBS_END=yyyy-mm-dd stops the daily tables early (default 2014-06-01).
using NumericalEarth
using NumericalEarth.DataWrangling: jldopen, BoundingBox
using NumericalEarth.DataWrangling.GLORYS: GLORYSDaily
using Dates, Printf
const NCD = NumericalEarth.DataWrangling.NCDatasets

rundir, tag = ARGS[1], ARGS[2]
obsdir = length(ARGS) ≥ 3 ? ARGS[3] : joinpath(homedir(), "Data", "sea_ice_obs")
datadir = length(ARGS) ≥ 4 ? ARGS[4] : joinpath(homedir(), "Data", "NumericalEarth")
bathyfile = length(ARGS) ≥ 5 ? ARGS[5] : joinpath(rundir, "bathy12_bathymetry.jld2")
region = BoundingBox(longitude = (-71.5, -48.0), latitude = (45.0, 56.0))
start = DateTime(2013, 12, 1)

const Nλ, Nφ = 278, 128
const λ₀, φ₀ = -71.5 + 2 / 12, 45.0 + 2 / 12   # the grid's west and south edges
λc = [λ₀ + (i - 0.5) / 12 for i in 1:Nλ]
φc = [φ₀ + (j - 0.5) / 12 for j in 1:Nφ]
area = [(6.371e6 * deg2rad(1 / 12))^2 * cosd(φc[j]) for i in 1:Nλ, j in 1:Nφ]
wet = jldopen(f -> f["bottom_height"] .< 0, bathyfile)
inner = [8 < i ≤ Nλ - 8 && 8 < j ≤ Nφ - 8 for i in 1:Nλ, j in 1:Nφ]
labrador = [φc[j] ≥ 51.7 for i in 1:Nλ, j in 1:Nφ]
gulf = [!labrador[i, j] && ((λc[i] < -59 && φc[j] > 46) || (-59 ≤ λc[i] < -56.5 && 47.5 < φc[j])) for i in 1:Nλ, j in 1:Nφ]
shelf = .!labrador .& .!gulf
# Lake Melville and the narrows to Rigolet (west of 58.3°W; the Labrador coast is east of 58°W at these latitudes)
melville = [-61 ≤ λc[i] ≤ -58.3 && 53 ≤ φc[j] ≤ 54.35 for i in 1:Nλ, j in 1:Nφ]
regions = (("gulf", gulf), ("labrador", labrador), ("melville", melville), ("shelf", shelf))
masks = [(label, wet .& inner .& m) for (label, m) in regions]

# ---------------- model: daily means of the 6-hourly ice output ----------------
# the interior of an output array (halo widths inferred from its size)
function cut(a)
    a = ndims(a) == 3 ? a[:, :, 1] : a
    hx, hy = (size(a, 1) - Nλ) ÷ 2, (size(a, 2) - Nφ) ÷ 2
    return Float64.(a[hx+1:hx+Nλ, hy+1:hy+Nφ])
end

model_days = jldopen(joinpath(rundir, tag * "_sea_ice_rank0.jld2")) do f
    its = sort(parse.(Int, filter(!=("serialized"), collect(keys(f["timeseries/h"])))))
    days = Dict{Date, Any}()
    for it in its
        t = start + Millisecond(round(Int, 1000 * f["timeseries/t/$it"]))
        ℵ, h = cut(f["timeseries/ℵ/$it"]), cut(f["timeseries/h/$it"])
        n, ℵs, vs = get(days, Date(t), (0, zeros(Nλ, Nφ), zeros(Nλ, Nφ)))
        days[Date(t)] = (n + 1, ℵs .+ ℵ, vs .+ ℵ .* h)
    end
    # the last day only if complete (00, 06, 12, 18 h)
    Dict(d => (; ℵ = ℵs ./ n, v = vs ./ n) for (d, (n, ℵs, vs)) in days if n ≥ 4)
end

# ---------------- GLORYS ----------------
function glorys(var, name, date)
    file = joinpath(datadir, NumericalEarth.DataWrangling.metadata_filename(GLORYSDaily(), name, DateTime(date), region))
    isfile(file) || return fill(NaN, Nλ, Nφ)
    NCD.Dataset(file) do ds
        data = Float64.(coalesce.(ds[var][:, :, 1], 0.0))
        λg, φg = Float64.(ds["longitude"][:]), Float64.(ds["latitude"][:])
        ig = [argmin(abs.(λg .- λ)) for λ in λc]; jg = [argmin(abs.(φg .- φ)) for φ in φc]
        [data[ig[i], jg[j]] for i in 1:Nλ, j in 1:Nφ]
    end
end

# ---------------- satellite grids: the pixels near each model cell ----------------
# For each wet model cell, the pixels of a (lon, lat) grid within `radius` metres of its centre, nearest first, among
# those where `keep` is true.
function neighbours(lon, lat, radius, keep = trues(size(lat)))
    lon = mod.(lon .+ 180, 360) .- 180
    near = findall(k -> keep[k] && 43 ≤ lat[k] ≤ 58 && -76 ≤ lon[k] ≤ -43, eachindex(lat))
    R = 6.371e6
    map(CartesianIndices((Nλ, Nφ))) do I
        i, j = Tuple(I)
        wet[i, j] || return Int[]
        d = [hypot(R * cosd(φc[j]) * deg2rad(lon[k] - λc[i]), R * deg2rad(lat[k] - φc[j])) for k in near]
        o = sortperm(d)
        [near[k] for k in o if d[k] ≤ radius]
    end
end

# A field on the model grid from a pixel field `p` (NaN where invalid): the nearest valid pixel, NaN if none.
function onmodel(p, nb)
    map(nb) do pixels
        k = findfirst(k -> !isnan(p[k]), pixels)
        k === nothing ? NaN : p[pixels[k]]
    end
end

nanify(a, ok) = [ok(x) ? Float64(x) : NaN for x in a]

# Pixels with no land within w pixels in either direction. Passive-microwave concentrations near coasts carry land
# spill-over (in early December 2013 both products show 0.2-0.6 along every Gulf coastline that the CIS charts have
# ice free), so pixels within about 50 km of land are left out: the model cells there have no satellite value.
away_from_land(ocean, w) = [all(ocean[max(1, i - w):min(end, i + w), max(1, j - w):min(end, j + w)])
                            for i in axes(ocean, 1), j in axes(ocean, 2)]

const OSI = joinpath(obsdir, "osisaf_osi450a_v3")
const CDR = joinpath(obsdir, "nsidc_g02202_v6")
osi_file(d) = joinpath(OSI, "ice_conc_nh_ease2-250_cdr-v3p0_$(Dates.format(d, "yyyymmdd"))1200.nc")
function cdr_file(d)
    f = filter(x -> startswith(x, "sic_psn25_$(Dates.format(d, "yyyymmdd"))_"), readdir(CDR))
    isempty(f) ? "" : joinpath(CDR, f[1])
end

# OSI-450-a: ice_conc is missing over land and lakes (two pixels ≈ 50 km); G02202: surface_type 50 is ocean, and its
# coast class already marks the first pixel off land (one more pixel)
nb_osi = NCD.Dataset(osi_file(Date(2013, 12, 1))) do ds
    keep = away_from_land(.!ismissing.(ds["ice_conc"][:, :, 1]), 2)
    neighbours(Float64.(ds["lon"][:, :]), Float64.(ds["lat"][:, :]), 25e3, keep)
end
nb_cdr = NCD.Dataset(joinpath(CDR, "G02202-ancillary-psn25-v06r00.nc")) do ds
    keep = away_from_land(ds["surface_type"][:, :] .== 50, 1)
    neighbours(Float64.(ds["longitude"][:, :]), Float64.(ds["latitude"][:, :]), 25e3, keep)
end

# OSI-450-a ice_conc in %; missing over land and lakes
osi(d) = isfile(osi_file(d)) ? NCD.Dataset(ds -> onmodel(nanify(ds["ice_conc"][:, :, 1] ./ 100, x -> !ismissing(x)), nb_osi),
                                         osi_file(d)) : fill(NaN, Nλ, Nφ)
# G02202 cdr_seaice_conc as a fraction; flagged values (pole hole, lake, coast, land) are above 1
cdr(d) = (f = cdr_file(d); isempty(f) ? fill(NaN, Nλ, Nφ) :
          NCD.Dataset(ds -> onmodel(nanify(ds["cdr_seaice_conc"][:, :, 1], x -> !ismissing(x) && x ≤ 1), nb_cdr), f))

# ---------------- Canadian Ice Service charts (SIGRID-3 shapefiles) ----------------
# The charts' Lambert conformal conic (WGS84; standard parallels 49 and 77°N, central meridian 100°W, origin 40°N)
const LCC = let a = 6378137.0, e = sqrt(2 / 298.257223563 - 1 / 298.257223563^2), φ1 = deg2rad(49.0), φ2 = deg2rad(77.0)
    m(φ) = cos(φ) / sqrt(1 - e^2 * sin(φ)^2)
    t(φ) = tan(π / 4 - φ / 2) / ((1 - e * sin(φ)) / (1 + e * sin(φ)))^(e / 2)
    n = log(m(φ1) / m(φ2)) / log(t(φ1) / t(φ2))
    F = m(φ1) / (n * t(φ1)^n)
    (; a, e, n, F, ρ₀ = a * F * t(deg2rad(40.0))^n, λ₀ = deg2rad(-100.0))
end

# chart (x, y) → (longitude, latitude) in degrees
function lcc_inverse(x, y)
    (; a, e, n, F, ρ₀, λ₀) = LCC
    ρ = hypot(x, ρ₀ - y)
    t = (ρ / (a * F))^(1 / n)
    φ = π / 2 - 2atan(t)
    for _ in 1:8
        φ = π / 2 - 2atan(t * ((1 - e * sin(φ)) / (1 + e * sin(φ)))^(e / 2))
    end
    return rad2deg(atan(x, ρ₀ - y) / n + λ₀), rad2deg(φ)
end

# Polygons (type 5) of a .shp file, each as its rings' vertices
function read_shp(path)
    b = read(path)
    i32(o) = reinterpret(Int32, b[o+1:o+4])[1]
    polygons = Vector{Vector{Matrix{Float64}}}()
    o = 100
    while o + 8 ≤ length(b)
        len = 2 * Int(ntoh(i32(o + 4)))
        c = o + 8
        rings = Matrix{Float64}[]
        if i32(c) == 5
            np, npt = i32(c + 36), i32(c + 40)
            parts = [Int(i32(c + 44 + 4k)) for k in 0:np-1]
            xy = reshape(collect(reinterpret(Float64, b[c+44+4np+1:c+44+4np+16npt])), 2, npt)
            for (k, p) in enumerate(parts)
                q = k < np ? parts[k+1] : npt
                push!(rings, xy[:, p+1:q])
            end
        end
        push!(polygons, rings)
        o = c + len
    end
    return polygons
end

# The attribute table of a .dbf file, as a vector of Dicts of strings
function read_dbf(path)
    b = read(path)
    n = Int(reinterpret(UInt32, b[5:8])[1]); hl = Int(reinterpret(UInt16, b[9:10])[1])
    rl = Int(reinterpret(UInt16, b[11:12])[1])
    fields = Tuple{String, Int}[]
    o = 32
    while b[o+1] != 0x0d
        push!(fields, (strip(String(filter(!=(0x00), b[o+1:o+11]))), Int(b[o+17])))
        o += 32
    end
    map(0:n-1) do r
        p = hl + r * rl + 1
        Dict(map(fields) do (name, L)
            v = strip(String(b[p+1:p+L])); p += L; name => v
        end)
    end
end

# The polygon each model cell centre falls in (0 if none): scan lines along the model rows in longitude-latitude, the
# even-odd rule over all of a polygon's rings (so holes stay out)
function cell_polygons(polygons)
    id = zeros(Int, Nλ, Nφ)
    for (k, rings) in enumerate(polygons)
        crossings = [Float64[] for _ in 1:Nφ]
        for r in rings
            ll = [lcc_inverse(r[1, m], r[2, m]) for m in 1:size(r, 2)]
            for m in 1:length(ll)-1
                (λ1, φ1), (λ2, φ2) = ll[m], ll[m+1]
                lo, hi = minmax(φ1, φ2)
                for j in max(1, ceil(Int, (lo - φ₀) * 12 + 0.5)):min(Nφ, ceil(Int, (hi - φ₀) * 12 + 0.5) - 1)
                    (φ1 ≤ φc[j]) == (φ2 ≤ φc[j]) && continue
                    push!(crossings[j], λ1 + (φc[j] - φ1) * (λ2 - λ1) / (φ2 - φ1))
                end
            end
        end
        for j in 1:Nφ
            c = sort(crossings[j])
            for p in 1:2:length(c)-1
                for i in max(1, ceil(Int, (c[p] - λ₀) * 12 + 0.5)):min(Nλ, ceil(Int, (c[p+1] - λ₀) * 12 + 0.5) - 1)
                    id[i, j] = k
                end
            end
        end
    end
    return id
end

# SIGRID-3 concentration codes: 00 ice free, 01 < 1/10, 02 bergy water, X0 X/10, 91 9+/10, 92 10/10, XY (X < Y) X to Y
# tenths, XY (X > Y) X to 10 tenths
function sigrid_concentration(code)
    (code == "" || code == "-9" || code == "99") && return NaN
    c = parse(Int, code)
    c ≤ 2 && return 0.0
    c == 91 && return 0.95
    c == 92 && return 1.0
    X, Y = divrem(c, 10)
    Y == 0 ? X / 10 : X < Y ? (X + Y) / 20 : (X + 10) / 20
end

# WMO stages of development → thickness (m): new ice and nilas, young, grey, grey-white, first-year (thin, first and
# second stage, medium, thick), old
const STAGE_THICKNESS = Dict(81 => 0.05, 82 => 0.05, 83 => 0.2, 84 => 0.125, 85 => 0.225, 86 => 0.7, 87 => 0.5,
                             88 => 0.4, 89 => 0.6, 91 => 0.95, 93 => 1.5, 95 => 2.5, 96 => 2.5, 97 => 3.0)
stage_thickness(code) = (code == "" || code == "-9") ? NaN : get(STAGE_THICKNESS, parse(Int, code), NaN)

# Concentration and ice volume per area (Σ Cₖ hₖ; NaN if a stage is unknown) of a chart polygon
function polygon_ice(r)
    r["POLY_TYPE"] == "L" && return (NaN, NaN)
    r["POLY_TYPE"] == "N" && return (NaN, NaN)
    r["POLY_TYPE"] == "W" && return (0.0, 0.0)
    C = sigrid_concentration(r["CT"])
    C == 0 && return (0.0, 0.0)
    partial = [(sigrid_concentration(r["C$s"]), stage_thickness(r["S$s"])) for s in ("A", "B", "C")]
    partial = filter(p -> !isnan(p[1]), partial)
    v = isempty(partial) ? C * stage_thickness(r["SA"]) : sum(c * h for (c, h) in partial)
    return C, v
end

const CIS = joinpath(obsdir, "cis_g02171")
# The latest version of each weekly chart of a region (East_Coast), by date
function chart_files(r)
    files = sort(filter(endswith(".tar"), readdir(joinpath(CIS, r))))
    Dict(Date(match(r"_(\d{8})_", f)[1], dateformat"yyyymmdd") => joinpath(CIS, r, f) for f in files)
end

# Concentration and volume per area of one chart on the model grid (NaN where no chart or land)
function chart_on_model(tarfile)
    dir = mktempdir()
    run(`tar -xf $tarfile -C $dir`)
    shp = first(filter(endswith(".shp"), readdir(dir; join = true)))
    polygons, records = read_shp(shp), read_dbf(replace(shp, ".shp" => ".dbf"))
    id = cell_polygons(polygons)
    rm(dir; recursive = true)
    ice = [polygon_ice(r) for r in records]
    C = [id[I] == 0 ? NaN : ice[id[I]][1] for I in CartesianIndices(id)]
    v = [id[I] == 0 ? NaN : ice[id[I]][2] for I in CartesianIndices(id)]
    land = [id[I] != 0 && records[id[I]]["POLY_TYPE"] == "L" for I in CartesianIndices(id)]
    return C, v, land
end

# The model's 1/12° coastline puts ~9% of its water cells on chart land (most of them along the coasts, and much of
# Lake Melville, which the model draws wider): those take the nearest chart water or ice cell within `reach` cells.
function reach_out(C, v, reach = 2)
    C′, v′ = copy(C), copy(v)
    offsets = sort([(a, b) for a in -reach:reach, b in -reach:reach if (a, b) != (0, 0)], by = o -> hypot(o...))
    for I in CartesianIndices(C)
        i, j = Tuple(I)
        (wet[i, j] && isnan(C[i, j])) || continue
        for (a, b) in offsets
            (1 ≤ i + a ≤ Nλ && 1 ≤ j + b ≤ Nφ && !isnan(C[i+a, j+b])) || continue
            C′[i, j], v′[i, j] = C[i+a, j+b], v[i+a, j+b]
            break
        end
    end
    return C′, v′
end

# ---------------- CryoSat-2/SMOS ----------------
const CS2SMOS = joinpath(obsdir, "cs2smos_v300")
function cs2smos_file(d)
    isdir(CS2SMOS) || return ""
    window = Dates.format(d - Day(3), "yyyymmdd") * "_" * Dates.format(d + Day(3), "yyyymmdd")
    f = filter(x -> occursin(window, x) && endswith(x, ".nc"), readdir(CS2SMOS))
    isempty(f) ? "" : joinpath(CS2SMOS, f[1])
end
nb_cs2 = Ref{Any}(nothing)
# Concentration and volume per area on the model grid; open ocean is ice free, land and failed retrievals NaN
function cs2smos(d)
    f = cs2smos_file(d)
    isempty(f) && return fill(NaN, Nλ, Nφ), fill(NaN, Nλ, Nφ)
    NCD.Dataset(f) do ds
        nb_cs2[] === nothing && (nb_cs2[] = neighbours(Float64.(ds["lon"][:, :]), Float64.(ds["lat"][:, :]), 12.5e3))
        status = ds["status_flag"][:, :, 1]
        sic = Float64.(coalesce.(ds["sea_ice_concentration"][:, :, 1], NaN)) ./ 100
        sit = Float64.(coalesce.(ds["sea_ice_thickness"][:, :, 1], NaN))
        open = status .== 1
        C = ifelse.(open, 0.0, sic)
        V = ifelse.(open, 0.0, sic .* sit)
        bad = .!open .& (isnan.(C) .| isnan.(V))
        C[bad] .= NaN; V[bad] .= NaN
        onmodel(C, nb_cs2[]), onmodel(V, nb_cs2[])
    end
end

# ---------------- daily comparison ----------------
# Sums over the cells of mask m; the satellite products are compared on the cells they cover (m .& covered), with the
# model and GLORYS summed over the same cells
kkm2(x) = x / 1e9                             # m² → 10³ km²
areasum(c, m) = all(isnan, c) || !any(m) ? NaN : kkm2(sum(c[I] * area[I] for I in eachindex(c) if m[I] && !isnan(c[I]); init = 0.0))
extent(c, m) = all(isnan, c) || !any(m) ? NaN : kkm2(sum(area[I] for I in eachindex(c) if m[I] && !isnan(c[I]) && c[I] ≥ 0.15; init = 0.0))
covered(c) = .!isnan.(c)
coverage(c, m) = 100 * sum(area[m .& covered(c)]; init = 0.0) / sum(area[m])
fmt(xs) = (@sprintf("%.3f", x) for x in xs)

days = Date(2013, 12, 1):Day(1):Date(get(ENV, "ICE_OBS_END", "2014-06-01"))
nanmap() = fill(NaN32, Nλ, Nφ, length(days))
maps = (; model = nanmap(), model_h = nanmap(), glorys = nanmap(), glorys_h = nanmap(), osi = nanmap(), cdr = nanmap())
csv = open(joinpath(rundir, tag * "_ice_obs_vs_model.csv"), "w")
println(csv, "date,region,area_model,area_glorys," *
             "osi_area_model,osi_area_glorys,osi_area_obs,osi_extent_model,osi_extent_glorys,osi_extent_obs," *
             "cdr_area_model,cdr_area_glorys,cdr_area_obs,cdr_extent_model,cdr_extent_glorys,cdr_extent_obs," *
             "coverage_osi,coverage_cdr")
println("10³ km²: area (Σ ℵ A) of model and GLORYS over the region | on the OSI-450-a cells: area of model, GLORYS, OSI; " *
        "extent (ℵ ≥ 0.15) of the same | on the G02202 cells: the same | % of the region's water those cells cover")
for (n, d) in enumerate(days)
    md = get(model_days, d, nothing)
    ℵm = md === nothing ? fill(NaN, Nλ, Nφ) : md.ℵ
    hm = md === nothing ? fill(NaN, Nλ, Nφ) : md.v ./ max.(md.ℵ, 1e-6)
    ℵg, hg = glorys("siconc", :sea_ice_concentration, d), glorys("sithick", :sea_ice_thickness, d)
    ℵo, ℵc = osi(d), cdr(d)
    for (k, a) in enumerate((ℵm, hm, ℵg, hg, ℵo, ℵc))
        maps[k][:, :, n] .= a
    end
    for (label, m) in masks
        full = [areasum(c, m) for c in (ℵm, ℵg)]
        o, c = m .& covered(ℵo), m .& covered(ℵc)
        on_osi = [[areasum(x, o) for x in (ℵm, ℵg, ℵo)]; [extent(x, o) for x in (ℵm, ℵg, ℵo)]]
        on_cdr = [[areasum(x, c) for x in (ℵm, ℵg, ℵc)]; [extent(x, c) for x in (ℵm, ℵg, ℵc)]]
        cov = [coverage(ℵo, m), coverage(ℵc, m)]
        println(csv, join([string(d), label, fmt((full..., on_osi..., on_cdr..., cov...))...], ","))
        @printf("%s %-8s %5.1f %5.1f | %5.1f %5.1f %5.1f  %5.1f %5.1f %5.1f | %5.1f %5.1f %5.1f  %5.1f %5.1f %5.1f | %3.0f%% %3.0f%%\n",
                d, label, full..., on_osi..., on_cdr..., cov...)
    end
end
close(csv)

# ---------------- weekly charts ----------------
# The East Coast charts cover the whole box; the Eastern Arctic charts (also downloaded) lie north of it.
charts = chart_files("East_Coast")
csv = open(joinpath(rundir, tag * "_ice_obs_vs_model_cis.csv"), "w")
println(csv, "date,region,area_model,area_glorys,area_cis,extent_model,extent_glorys,extent_cis," *
             "volume_model,volume_glorys,volume_cis,cs2smos_volume_model,cs2smos_volume_glorys,cs2smos_volume_obs," *
             "coverage_cis,coverage_cs2smos")
println("\nCIS charts, on the cells they cover: area | extent (10³ km²) | volume (km³) of model, GLORYS, CIS | on the " *
        "CryoSat-2/SMOS cells (window centred on the chart date): volume of model, GLORYS, CS2SMOS | % covered by CIS, CS2SMOS")
chart_dates = sort(filter(d -> first(days) ≤ d ≤ last(days), collect(keys(charts))))
chart_maps = Dict{String, Any}()
for (k, d) in enumerate(chart_dates)
    C, v, land = chart_on_model(charts[d])
    k == 1 && @printf("chart %s: %.0f%% of the model's water cells fall on chart land, %.0f%% of its land cells on chart water or ice\n",
                      d, 100 * count(land .& wet) / count(wet), 100 * count(.!land .& .!isnan.(C) .& .!wet) / count(.!wet))
    C, v = reach_out(C, v)
    n = findfirst(==(d), days)
    md = get(model_days, d, nothing)
    ℵm = md === nothing ? fill(NaN, Nλ, Nφ) : md.ℵ
    vm = md === nothing ? fill(NaN, Nλ, Nφ) : md.v
    ℵg, hg = Float64.(maps.glorys[:, :, n]), Float64.(maps.glorys_h[:, :, n])
    vg = ℵg .* hg
    Cs, vs = cs2smos(d)
    for (label, m) in masks
        c, s = m .& covered(C), m .& covered(vs)
        A = [areasum(x, c) for x in (ℵm, ℵg, C)]
        E = [extent(x, c) for x in (ℵm, ℵg, C)]
        V = [areasum(x, c) for x in (vm, vg, v)]       # 10³ km² × m = km³
        Vs = [areasum(x, s) for x in (vm, vg, vs)]
        cov = [coverage(C, m), coverage(vs, m)]
        println(csv, join([string(d), label, fmt((A..., E..., V..., Vs..., cov...))...], ","))
        @printf("%s %-8s %5.1f %5.1f %5.1f | %5.1f %5.1f %5.1f | %5.1f %5.1f %5.1f | %5.1f %5.1f %5.1f | %3.0f%% %3.0f%%\n",
                d, label, A..., E..., V..., Vs..., cov...)
    end
    chart_maps[string(d)] = (; cis_C = C, cis_v = v, cs2smos_C = Cs, cs2smos_v = vs)
end
close(csv)

jldopen(joinpath(rundir, tag * "_ice_obs_maps.jld2"), "w") do file
    file["λ"] = λc; file["φ"] = φc; file["wet"] = wet; file["dates"] = collect(days)
    for (k, v) in pairs(maps)
        file[string(k)] = v
    end
    for (d, m) in chart_maps, (k, v) in pairs(m)
        file["charts/$d/$k"] = v
    end
end
