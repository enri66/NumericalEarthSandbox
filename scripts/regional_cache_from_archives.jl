# Cut the global forcing archives on antares into NumericalEarth's per-region cache files, so a new domain or period
# needs no download:
#   GLORYS12 daily (A3, Mercator's daily files, all variables in one file)  → <var>_GLORYSDaily_<date>_<date>_<region>.nc
#   ERA5 hourly (A2, yearly or monthly global files)                         → <name>_ERA5HourlySingleLevel_<hour>_<region>.nc
#   GloFAS v4 daily (A3, monthly global files)                               → average_river_discharge_in_the_last_24_hours_GloFASReanalysis_<date>_<region>.nc
# The files have the layout NumericalEarth writes when it downloads them itself (the padded box, the coordinate order and
# the variable names), checked against downloaded MAB files with CACHE_COMPARE (see the end of this file).
#
# Environment:
#   CACHE_LON, CACHE_LAT   the region as NumericalEarth's BoundingBox sees it, e.g. "-70.5,-48" and "45,56"
#   CACHE_EXTRA_REGIONS    more regions cut for the first day only (GLORYS), e.g. the model box script 05 uses for the
#                          initial state: "-70.33333,-48.16667,45.16667,55.83333" (several separated by ';')
#   CACHE_START, CACHE_END first and last day (inclusive)
#   CACHE_OUT              the cache directory (GLORYS files go in it, ERA5 in <out>/era5, GloFAS in <out>/glofas)
#   CACHE_WHAT             any of "glorys,era5,glofas" (default all)
#   CACHE_GLORYS_VARS      default "thetao,so,uo,vo,zos,siconc,sithick,usi,vsi"
#   CACHE_ERA5_VARS        CDS names, default the ones ERA5PrescribedAtmosphere/Radiation read
#   CACHE_COMPARE          a directory of NumericalEarth-downloaded files to compare the cut files with, instead of cutting
using NCDatasets, Dates, Printf

const GLORYS_DIR = "/Volumes/A3/enrique/GLORYS12"
const ERA5_DIR   = "/Volumes/A2/ERA5"
const GLOFAS_DIR = "/Volumes/A3/enrique/GloFAS"

parse_pair(s) = Tuple(parse.(Float64, split(s, ",")))
const LON = parse_pair(get(ENV, "CACHE_LON", "-70.5,-48"))
const LAT = parse_pair(get(ENV, "CACHE_LAT", "45,56"))
const EXTRA = [Tuple(parse.(Float64, split(r, ","))) for r in split(get(ENV, "CACHE_EXTRA_REGIONS", ""), ";") if !isempty(r)]
const DAY1 = Date(get(ENV, "CACHE_START", "2013-11-25"))
const DAY2 = Date(get(ENV, "CACHE_END",   "2014-06-05"))
const OUT  = get(ENV, "CACHE_OUT", joinpath(homedir(), "Data", "NumericalEarth_cache_cut"))
const WHAT = split(get(ENV, "CACHE_WHAT", "glorys,era5,glofas"), ",")
const GLORYS_VARS = split(get(ENV, "CACHE_GLORYS_VARS", "thetao,so,uo,vo,zos,siconc,sithick,usi,vsi"), ",")

# CDS name => short name in the files
const ERA5_NAMES = Dict("10m_u_component_of_wind" => "u10", "10m_v_component_of_wind" => "v10",
                        "2m_temperature" => "t2m", "2m_dewpoint_temperature" => "d2m",
                        "surface_pressure" => "sp", "mean_sea_level_pressure" => "msl",
                        "total_precipitation" => "tp", "snowfall" => "sf",
                        "surface_solar_radiation_downwards" => "ssrd", "surface_thermal_radiation_downwards" => "strd")
const ERA5_VARS = split(get(ENV, "CACHE_ERA5_VARS",
    "10m_u_component_of_wind,10m_v_component_of_wind,2m_temperature,2m_dewpoint_temperature,surface_pressure," *
    "total_precipitation,surface_solar_radiation_downwards,surface_thermal_radiation_downwards,snowfall"), ",")

suffix(lon, lat) = @sprintf("_%.1f_%.1f_%.1f_%.1f", lon[1], lon[2], lat[1], lat[2])

# indices of the coordinate values within [a - pad, b + pad] (a small tolerance for Float32 coordinates)
function window(c, a, b, pad)
    tol = 1e-4
    idx = findall(x -> a - pad - tol <= x <= b + pad + tol, c)
    isempty(idx) && error("no coordinate in [$a, $b] ± $pad")
    return first(idx):last(idx)
end

#####
##### GLORYS: the CopernicusMarine extension requests the box padded by two native cells (1/6°)
#####

function glorys_file(d::Date)
    dir = joinpath(GLORYS_DIR, Dates.format(d, "yyyy"), Dates.format(d, "mm"))
    tag = "_mean_" * Dates.format(d, "yyyymmdd") * "_"
    f = filter(x -> occursin(tag, x), readdir(dir))
    length(f) == 1 || error("expected one GLORYS file for $d in $dir, found $(length(f))")
    return joinpath(dir, f[1])
end

glorys_name(var, d, lon, lat) = string(var, "_GLORYSDaily_", Dates.format(DateTime(d), "yyyy-mm-ddTHH-MM-SS"), "_",
                                       Dates.format(DateTime(d), "yyyy-mm-ddTHH-MM-SS"), suffix(lon, lat), ".nc")

copyattribs(v) = Dict(k => v.attrib[k] for k in keys(v.attrib))

function cut_glorys(d::Date, regions, vars)
    NCDataset(glorys_file(d)) do src
        λ, φ = src["longitude"][:], src["latitude"][:]
        for (lon, lat) in regions
            ii, jj = window(λ, lon..., 2/12), window(φ, lat..., 2/12)
            for var in vars
                path = joinpath(OUT, glorys_name(var, d, lon, lat))
                isfile(path) && continue
                tmp = path * ".part"
                sv = src[var]
                three_d = ndims(sv) == 4
                NCDataset(tmp, "c") do dst
                    three_d && (defDim(dst, "depth", src.dim["depth"]))
                    defDim(dst, "latitude", length(jj)); defDim(dst, "longitude", length(ii)); defDim(dst, "time", 1)
                    defVar(dst, "longitude", λ[ii], ("longitude",); attrib = copyattribs(src["longitude"]))
                    defVar(dst, "latitude",  φ[jj], ("latitude",);  attrib = copyattribs(src["latitude"]))
                    three_d && defVar(dst, "depth", src["depth"][:], ("depth",); attrib = copyattribs(src["depth"]))
                    tv = defVar(dst, "time", Float32, ("time",); attrib = copyattribs(src["time"]))
                    # Mercator's files stamp the daily mean at noon, the CopernicusMarine subsets at 00:00 of the day
                    tv.var[:] = [Float32(Dates.value(DateTime(d) - DateTime(1950)) ÷ 3_600_000)]
                    dims = three_d ? ("longitude", "latitude", "depth", "time") : ("longitude", "latitude", "time")
                    v = defVar(dst, var, Int16, dims; attrib = copyattribs(sv), deflatelevel = 1)
                    raw = three_d ? sv.var[ii, jj, :, 1:1] : sv.var[ii, jj, 1:1]
                    three_d ? (v.var[:, :, :, :] = raw) : (v.var[:, :, :] = raw)
                end
                mv(tmp, path; force = true)
            end
        end
    end
end

#####
##### ERA5: NumericalEarth pads the box by 1/2° (two native cells); longitudes −180…180, latitudes decreasing
#####

function era5_source(name, t::DateTime)
    y = Dates.format(t, "yyyy"); ym = Dates.format(t, "yyyy-mm")
    yearly, monthly = joinpath(ERA5_DIR, "ERA5_$(name)_$(y).nc"), joinpath(ERA5_DIR, "ERA5_$(name)_$(ym).nc")
    isfile(monthly) && return monthly, :monthly
    isfile(yearly) && filesize(yearly) > 18e9 && return yearly, :yearly
    error("no complete ERA5 file for $name at $t")
end

era5_name(name, t, lon, lat) = string(name, "_ERA5HourlySingleLevel_", Dates.format(t, "yyyy-mm-ddTHH"), suffix(lon, lat), ".nc")

# one day of hours from one source file
function cut_era5_day(name, d::Date, lon, lat)
    short = ERA5_NAMES[name]
    hours = DateTime(d):Hour(1):DateTime(d) + Hour(23)
    todo = [t for t in hours if !isfile(joinpath(OUT, "era5", era5_name(name, t, lon, lat)))]
    isempty(todo) && return
    path, kind = era5_source(name, first(todo))
    NCDataset(path) do src
        tname = haskey(src, "valid_time") ? "valid_time" : "time"
        times = src[tname][:]
        λ = Float64.(src["longitude"][:]); φ = Float64.(src["latitude"][:])
        λ180 = ifelse.(λ .> 180, λ .- 360, λ)
        ii = window(λ180, lon..., 0.5)
        issorted(λ180[ii]) || error("the ERA5 longitude window crosses the source's seam")
        jj = window(φ, lat..., 0.5)
        for t in todo
            n = findfirst(==(t), times)
            n === nothing && error("$t not in $path")
            data = Float32.(coalesce.(src[short][ii, jj, n], NaN32))
            out = joinpath(OUT, "era5", era5_name(name, t, lon, lat))
            tmp = out * ".part"
            NCDataset(tmp, "c") do dst
                defDim(dst, "valid_time", 1); defDim(dst, "latitude", length(jj)); defDim(dst, "longitude", length(ii))
                defVar(dst, "valid_time", [Dates.value(t - DateTime(1970)) ÷ 1000], ("valid_time",);
                       attrib = Dict("long_name" => "time", "standard_name" => "time",
                                     "units" => "seconds since 1970-01-01", "calendar" => "proleptic_gregorian"))
                defVar(dst, "latitude", φ[jj], ("latitude",);
                       attrib = Dict("units" => "degrees_north", "standard_name" => "latitude", "stored_direction" => "decreasing"))
                defVar(dst, "longitude", λ180[ii], ("longitude",); attrib = Dict("units" => "degrees_east", "standard_name" => "longitude"))
                defVar(dst, short, reshape(data, size(data)..., 1), ("longitude", "latitude", "valid_time");
                       fillvalue = NaN32, attrib = Dict("units" => get(src[short].attrib, "units", ""),
                                                         "long_name" => get(src[short].attrib, "long_name", "")))
            end
            mv(tmp, out; force = true)
        end
    end
end

#####
##### GloFAS: daily mean discharge; the archive stamps each daily mean at the end of its window (D+1 00:00)
#####

const GLOFAS_PAD = parse(Float64, get(ENV, "CACHE_GLOFAS_PAD", "0.2"))
const GLOFAS_SHIFT = Day(parse(Int, get(ENV, "CACHE_GLOFAS_SHIFT", "1")))   # archive valid_time = file date + shift
glofas_name(d, lon, lat) = string("average_river_discharge_in_the_last_24_hours_GloFASReanalysis_",
                                  Dates.format(d, "yyyy-mm-dd"), suffix(lon, lat), ".nc")

function cut_glofas(d::Date, lon, lat)
    out = joinpath(OUT, "glofas", glofas_name(d, lon, lat))
    isfile(out) && return
    t = DateTime(d) + GLOFAS_SHIFT
    # the window for valid_time t lies in the file of the day it averages
    path = joinpath(GLOFAS_DIR, "GloFAS_v4_river_discharge_$(Dates.format(t - Day(1), "yyyy-mm")).nc")
    NCDataset(path) do src
        n = findfirst(==(t), src["valid_time"][:])
        n === nothing && error("$t not in $path")
        λ = src["longitude"][:]; φ = src["latitude"][:]
        ii, jj = window(λ, lon..., GLOFAS_PAD), window(φ, lat..., GLOFAS_PAD)
        data = Float32.(coalesce.(src["avg_dis"][ii, jj, n], NaN32))
        tmp = out * ".part"
        NCDataset(tmp, "c") do dst
            defDim(dst, "valid_time", 1); defDim(dst, "latitude", length(jj)); defDim(dst, "longitude", length(ii))
            defVar(dst, "valid_time", [Dates.value(t - DateTime(1970)) ÷ 1000], ("valid_time",);
                   attrib = Dict("units" => "seconds since 1970-01-01", "calendar" => "proleptic_gregorian"))
            defVar(dst, "latitude", φ[jj], ("latitude",); attrib = Dict("units" => "degrees_north", "stored_direction" => "decreasing"))
            # the EWDS subsets come with longitudes in 0…360
            defVar(dst, "longitude", mod.(λ[ii], 360), ("longitude",); attrib = Dict("units" => "degrees_east"))
            defVar(dst, "avg_dis", reshape(data, size(data)..., 1), ("longitude", "latitude", "valid_time");
                   fillvalue = NaN32, attrib = Dict("units" => "m**3 s**-1"))
        end
        mv(tmp, out; force = true)
    end
end

#####
##### Comparison with NumericalEarth-downloaded files
#####

function compare(reference_dir)
    for f in sort(filter(x -> endswith(x, ".nc"), readdir(OUT; join = true)) ∪
                  filter(x -> endswith(x, ".nc"), vcat([isdir(joinpath(OUT, s)) ? readdir(joinpath(OUT, s); join = true) : String[]
                                                        for s in ("era5", "glofas")]...)))
        ref = joinpath(reference_dir, relpath(f, OUT))
        isfile(ref) || (println("no reference for ", relpath(f, OUT)); continue)
        NCDataset(f) do a
            NCDataset(ref) do b
                name = first(k for k in keys(a) if k ∉ ("longitude", "latitude", "depth", "time", "valid_time"))
                sa, sb = size(a[name]), size(b[name])
                if sa != sb
                    println(relpath(f, OUT), ": SIZE ", sa, " vs reference ", sb)
                    return
                end
                dλ = maximum(abs, Float64.(a["longitude"][:]) .- Float64.(b["longitude"][:]))
                dφ = maximum(abs, Float64.(a["latitude"][:]) .- Float64.(b["latitude"][:]))
                x, y = Float64.(coalesce.(a[name][:], NaN)), Float64.(coalesce.(b[name][:], NaN))
                nan_mismatch = count(isnan.(x) .!= isnan.(y))
                ok = .!isnan.(x) .& .!isnan.(y)
                dmax = any(ok) ? maximum(abs, x[ok] .- y[ok]) : NaN
                scale = any(ok) ? maximum(abs, y[ok]) : NaN
                tname = haskey(a, "valid_time") ? "valid_time" : "time"
                @printf("%-95s size %-16s |Δλ| %.1e |Δφ| %.1e  max|Δ| %.3e (max|ref| %.3e)  NaN mismatches %d  time %s vs %s\n",
                        relpath(f, OUT), string(sa), dλ, dφ, dmax, scale, nan_mismatch, a[tname].var[1], b[tname].var[1])
            end
        end
    end
end

#####
##### Main
#####

mkpath(OUT); mkpath(joinpath(OUT, "era5")); mkpath(joinpath(OUT, "glofas"))

if haskey(ENV, "CACHE_COMPARE")
    compare(ENV["CACHE_COMPARE"])
else
    regions = [(LON, LAT)]
    for d in DAY1:Day(1):DAY2
        t0 = time()
        if "glorys" in WHAT
            cut_glorys(d, regions, GLORYS_VARS)
            d == DAY1 && cut_glorys(d, [((r[1], r[2]), (r[3], r[4])) for r in EXTRA], GLORYS_VARS)
        end
        "era5" in WHAT && foreach(name -> cut_era5_day(name, d, LON, LAT), ERA5_VARS)
        "glofas" in WHAT && cut_glofas(d, LON, LAT)
        @printf("%s done in %.1f s\n", d, time() - t0); flush(stdout)
    end
end
