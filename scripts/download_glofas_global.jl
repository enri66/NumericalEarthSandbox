# Downloads the global GloFAS v4 river discharge reanalysis (daily means, 0.05°), one file per month, from 1993 on,
# alongside the global ERA5 and GLORYS12 archives. GloFAS lives on the Early Warning Data Store (EWDS), not on CDS:
# the requests go to https://ewds.climate.copernicus.eu/api with the key of ~/.cdsapirc, and the `cems-glofas-historical`
# licence must have been accepted on the EWDS site. Requests go one at a time. It can be rerun at any time: finished
# months are skipped, and a month is written to `.part` and renamed only once it holds every day. Usage (on antares):
#   julia --project=. scripts/download_glofas_global.jl
# GLOFAS_DIR:   the archive (default /Volumes/A3/enrique/GloFAS)
# GLOFAS_START: first month, YYYY-MM (default 1993-01)
# GLOFAS_END:   last month, YYYY-MM (default: three months before the current one)
using CopernicusClimateDataStore, Dates, Printf
using NumericalEarth
const NCD = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "NCDatasets")]   # loaded by NumericalEarth
const CCDS = CopernicusClimateDataStore

const DIR   = get(ENV, "GLOFAS_DIR", "/Volumes/A3/enrique/GloFAS")
const FIRST = Date(get(ENV, "GLOFAS_START", "1993-01") * "-01")
const LAST  = haskey(ENV, "GLOFAS_END") ? Date(ENV["GLOFAS_END"] * "-01") : firstdayofmonth(today()) - Month(3)

const EWDS = CCDS.CDSCredentials("https://ewds.climate.copernicus.eu/api", CCDS.read_cds_credentials().key)

monthly_name(date) = joinpath(DIR, @sprintf("GloFAS_v4_river_discharge_%04d-%02d.nc", year(date), month(date)))

function days_in(path)
    ds = NCD.Dataset(path)
    n = length(ds["valid_time"])
    close(ds)
    return n
end

function download_month(date; attempts = 3)
    out = monthly_name(date)
    part = out * ".part"
    ndays = Dates.daysinmonth(date)
    params = Dict("system_version" => ["version_4_0"], "hydrological_model" => ["lisflood"],
                  "product_type" => ["consolidated"], "timespan" => ["time_mean"],
                  "variable" => ["average_river_discharge_in_the_last_24_hours"],
                  "year" => [string(year(date))], "month" => [lpad(month(date), 2, "0")],
                  "day" => [lpad(d, 2, "0") for d in 1:ndays],
                  "data_format" => "netcdf", "download_format" => "unarchived")
    for attempt in 1:attempts
        try
            t0 = time()
            retrieve("cems-glofas-historical", params, part; credentials = EWDS, max_wait = 12 * 3600, poll_interval = 30,
                     verbose = false)
            n = days_in(part)
            n == ndays || error("$(n) days in the file, expected $(ndays)")
            mv(part, out; force = true)
            @printf("%s  %s  %.2f GB  %.1f min\n", Dates.format(now(), "yyyy-mm-dd HH:MM"), Dates.format(date, "yyyy-mm"),
                    filesize(out) / 1e9, (time() - t0) / 60)
            return true
        catch e
            @warn "$(Dates.format(date, "yyyy-mm")), attempt $(attempt) of $(attempts): $(sprint(showerror, e))"
            isfile(part) && rm(part)
            attempt < attempts && sleep(600)
        end
    end
    return false
end

mkpath(DIR)
todo = [date for date in FIRST:Month(1):LAST if !isfile(monthly_name(date))]
@info "GloFAS archive $(DIR): $(Dates.format(FIRST, "yyyy-mm")) to $(Dates.format(LAST, "yyyy-mm")), $(length(todo)) months to download"
failed = Date[]
for date in todo
    isfile(monthly_name(date)) && continue          # done meanwhile by another run
    download_month(date) || push!(failed, date)
end
isempty(failed) || @warn "failed, rerun to retry: " * join(Dates.format.(failed, "yyyy-mm"), ", ")
@info "done"
