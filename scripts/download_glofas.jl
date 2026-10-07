# Download GloFAS river discharge (daily, 0.05°) for the MAB box, one month per EWDS request.
#
# GloFAS lives on the Early Warning Data Store, not the ERA5 CDS: NumericalEarth points CDSAPI at
# https://ewds.climate.copernicus.eu/api for these requests (the ECMWF token in ~/.cdsapirc is shared), and the
# `cems-glofas-historical` licence has to be accepted once in a browser on the EWDS site. The CDSAPI package must be
# in the active environment (it activates NumericalEarth's GloFAS download).
#
#   julia --project=~/glofas_env scripts/download_glofas.jl
#
# GLOFAS_START / GLOFAS_END (default 2019-08-29 / 2019-11-01) select the days; MAB_DATA_DIR the cache. Files already on
# disk are skipped. Needs JULIA_SSL_CA_ROOTS_PATH on antares (see the ERA5 download notes).
using NumericalEarth, Dates, Printf
using CDSAPI                       # activates the GloFAS download backend
using NumericalEarth.DataWrangling: Metadata, BoundingBox
import Downloads

const DATA_DIR = get(ENV, "MAB_DATA_DIR", joinpath(homedir(), "Data", "NumericalEarth"))
const dir      = joinpath(DATA_DIR, "glofas")
const region   = BoundingBox(longitude = (-76.0, -64.0), latitude = (34.0, 42.0))   # same box as script 05
const first_date = DateTime(get(ENV, "GLOFAS_START", "2019-08-29"))
const last_date  = DateTime(get(ENV, "GLOFAS_END",   "2019-11-01"))

mkpath(dir)
month_start = first_date
while month_start <= last_date
    month_end = min(DateTime(Dates.lastdayofmonth(month_start)), last_date)
    t = time()
    metadata = Metadata(:river_discharge; dataset = GloFASReanalysis(), start_date = month_start, end_date = month_end, dir, region)
    paths = Downloads.download(metadata)
    @printf("%s to %s: %d files in %.0f s\n", Dates.format(month_start, "yyyy-mm-dd"), Dates.format(month_end, "yyyy-mm-dd"),
            length(paths), time() - t)
    flush(stdout)
    global month_start = DateTime(Dates.firstdayofmonth(month_start)) + Month(1)
end
println("done: ", dir)
