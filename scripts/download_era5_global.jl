# Completes the global ERA5 archive (hourly single levels, 0.25°) from 1993 to the last complete month.
# The archive already holds yearly files `ERA5_<variable>_<year>.nc` for 1993-2018 (16-bit, downloaded in 2020);
# this adds monthly files `ERA5_<variable>_<YYYY-MM>.nc` (32-bit floats, as CDS delivers them now) for every month
# that no complete yearly file covers. Requests go to CDS one at a time, since CDS rejects more queued requests per
# dataset. It can be rerun at any time: finished months are skipped, and a month is written to `.part` and renamed
# only once it holds every hour. Usage (on antares):
#   julia --project=. scripts/download_era5_global.jl
# ERA5_DIR:       the archive (default /Volumes/A2/ERA5)
# ERA5_START:     first month, YYYY-MM (default 1993-01)
# ERA5_END:       last month, YYYY-MM (default: the last complete month)
# ERA5_VARIABLES: comma-separated CDS variable names (default: the forcing variables, then the archive's others)
using CopernicusClimateDataStore, Dates, Printf
using NumericalEarth
const NCD = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "NCDatasets")]   # loaded by NumericalEarth

const DIR = get(ENV, "ERA5_DIR", "/Volumes/A2/ERA5")

# The ocean model's atmospheric forcing, plus mean sea-level pressure; downloaded month by month first
const FORCING = ["10m_u_component_of_wind", "10m_v_component_of_wind", "2m_temperature", "2m_dewpoint_temperature",
                 "surface_pressure", "mean_sea_level_pressure", "surface_solar_radiation_downwards",
                 "surface_thermal_radiation_downwards", "total_precipitation"]
# The archive's other variables
const OTHERS = ["sea_surface_temperature", "sea_ice_cover", "total_cloud_cover", "snowfall", "convective_rain_rate",
                "large_scale_rain_rate", "convective_snowfall_rate_water_equivalent",
                "large_scale_snowfall_rate_water_equivalent"]
const GROUPS = haskey(ENV, "ERA5_VARIABLES") ? [String.(split(ENV["ERA5_VARIABLES"], ","))] : [FORCING, OTHERS]

const FIRST = Date(get(ENV, "ERA5_START", "1993-01") * "-01")
const LAST  = haskey(ENV, "ERA5_END") ? Date(ENV["ERA5_END"] * "-01") : firstdayofmonth(today()) - Month(1)

# A complete yearly file is about 18.2 GB (18.24 GB in leap years); the 2008 u-wind file is only 0.25 GB.
const COMPLETE_YEAR_BYTES = 18.0e9

# One yearly file in the archive carries a stray suffix in its name
yearly_names(var, year) = [joinpath(DIR, "ERA5_$(var)_$(year).nc"), joinpath(DIR, "ERA5_$(var)runoff_$(year).nc")]
monthly_name(var, date) = joinpath(DIR, @sprintf("ERA5_%s_%04d-%02d.nc", var, year(date), month(date)))

covered(var, date) = any(f -> isfile(f) && filesize(f) >= COMPLETE_YEAR_BYTES, yearly_names(var, year(date))) ||
                     isfile(monthly_name(var, date))

function hours_in(path)
    ds = NCD.Dataset(path)
    n = length(ds["valid_time"])
    close(ds)
    return n
end

function download_month(var, date; attempts = 3)
    out = monthly_name(var, date)
    part = out * ".part"
    ndays = Dates.daysinmonth(date)
    params = Dict("product_type" => ["reanalysis"], "variable" => [var],
                  "year" => [string(year(date))], "month" => [lpad(month(date), 2, "0")],
                  "day" => [lpad(d, 2, "0") for d in 1:ndays], "time" => [lpad(h, 2, "0") * ":00" for h in 0:23],
                  "data_format" => "netcdf", "download_format" => "unarchived")
    for attempt in 1:attempts
        try
            t0 = time()
            retrieve("reanalysis-era5-single-levels", params, part; max_wait = 12 * 3600, poll_interval = 30, verbose = false)
            n = hours_in(part)
            n == 24ndays || error("$(n) hours in the file, expected $(24ndays)")
            mv(part, out; force = true)
            @printf("%s  %-45s %s  %.2f GB  %.1f min\n", Dates.format(now(), "yyyy-mm-dd HH:MM"), var,
                    Dates.format(date, "yyyy-mm"), filesize(out) / 1e9, (time() - t0) / 60)
            return true
        catch e
            @warn "$(var) $(Dates.format(date, "yyyy-mm")), attempt $(attempt) of $(attempts): $(sprint(showerror, e))"
            isfile(part) && rm(part)
            attempt < attempts && sleep(600)
        end
    end
    return false
end

months = FIRST:Month(1):LAST
@info "ERA5 archive $(DIR): $(Dates.format(FIRST, "yyyy-mm")) to $(Dates.format(LAST, "yyyy-mm"))"
failed = Tuple{String, Date}[]
for group in GROUPS
    todo = [(var, date) for date in months for var in group if !covered(var, date)]
    @info "$(length(todo)) months to download for $(join(group, ", "))"
    for (var, date) in todo
        covered(var, date) && continue          # done meanwhile by another run
        download_month(var, date) || push!(failed, (var, date))
    end
end
isempty(failed) || @warn "failed, rerun to retry: " * join(["$(v) $(Dates.format(d, "yyyy-mm"))" for (v, d) in failed], ", ")
@info "done; the most recent three months or so are preliminary (ERA5T): to refresh them later, delete their files and rerun"
