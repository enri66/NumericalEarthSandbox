# Download the ERA5 forcing scripts 04/05 need for a long MAB run, one CDS request at a time.
#
# NumericalEarth's ERA5 download submits one request per variable concurrently (era5cli threads). CDS now caps
# the number of queued requests per user for ERA5 ("Number queued requests for this dataset is temporarily
# limited"), rejects the whole batch, and the client just polls the rejected jobs until its 1-hour timeout.
# `threads = 1` makes era5cli submit the variables one after another instead. Files already on disk are skipped,
# so this can be re-run after an interruption.
#
#   julia --project=. scripts/download_era5_year.jl
#
# ERA5_START / ERA5_END (default 2019-04-01T00 / 2020-04-01T01) and MAB_DATA_DIR select the range and cache;
# ERA5_NAMES (comma-separated NumericalEarth variable names) replaces the forcing variables.

using NumericalEarth, Dates, Printf
using CopernicusClimateDataStore    # activates the ERA5 download backend
using NumericalEarth.DataWrangling: MetadataSet, BoundingBox
using NumericalEarth.DataWrangling.ERA5: ERA5HourlySingleLevel
import Downloads

const DATA_DIR = get(ENV, "MAB_DATA_DIR", joinpath(homedir(), "Data", "NumericalEarth"))
const dir      = joinpath(DATA_DIR, "era5")
const region   = BoundingBox(longitude = (-76.0, -64.0), latitude = (34.0, 42.0))   # same box as script 04
const dataset  = ERA5HourlySingleLevel()

# ERA5PrescribedAtmosphere + ERA5PrescribedRadiation, or the comma-separated NumericalEarth names in ERA5_NAMES
const names = haskey(ENV, "ERA5_NAMES") ? Tuple(Symbol.(split(ENV["ERA5_NAMES"], ","))) :
              (:eastward_velocity, :northward_velocity, :temperature, :dewpoint_temperature,
               :surface_pressure, :total_precipitation,
               :downwelling_shortwave_radiation, :downwelling_longwave_radiation)

first_date = DateTime(get(ENV, "ERA5_START", "2019-04-01T00"))
last_date  = DateTime(get(ENV, "ERA5_END",   "2020-04-01T01"))

month_start = first_date
while month_start <= last_date
    month_end = min(firstdayofmonth(month_start) + Month(1) - Hour(1), last_date)
    for name in names
        t = time()
        mset = MetadataSet(name; dataset, start_date = month_start, end_date = month_end, dir, region)
        Downloads.download(mset; threads = 1)
        @printf("%s  %-32s %s → %s  (%.0f s)\n", Dates.format(now(), "yyyy-mm-dd HH:MM"), name,
                month_start, month_end, time() - t)
        flush(stdout)
    end
    global month_start = month_end + Hour(1)
end

println("✅ ERA5 complete: $first_date → $last_date in $dir")
