# Download the Copernicus Marine global BGC hindcast (CopernicusBGCMonthly) for the MAB with MAB_BGC_DATA=copernicus:
# NO₃, O₂, phytoplankton carbon and pH (3-D, total scale; for deriving DIC), one file per month from MAB_BGC_FIRST_MONTH to MAB_BGC_LAST_MONTH (default
# 2018-12 → 2020-01, which covers the mid-month frames of any run in 2019), into the directory and with the region
# that bgc_data.jl uses, so mab_bgc.jl finds them and downloads nothing. Serial; run it on a compute node:
#   srun -p triton24 -N 1 --exclude=node[38-41] julia --project=/t0/workdir/enrique/mpi05_bgc_cmems bgc/bgc_copernicus_download.jl

using NumericalEarth
using NumericalEarth.DataWrangling: Metadatum, BoundingBox, metadata_path
using CopernicusMarine
using Dates
import Downloads

const DATA_DIR = get(ENV, "MAB_DATA_DIR", joinpath(homedir(), "Data", "NumericalEarth"))
const data_λ = (-76.0, -64.0)   # mab_bgc.jl's GLORYS box
const data_φ = ( 34.0,  42.0)
const region = BoundingBox(longitude = (data_λ[1] - 1, data_λ[2] + 1), latitude = (data_φ[1] - 1, data_φ[2] + 1))
const dir    = joinpath(DATA_DIR, "copernicus_bgc")

first_month = DateTime(get(ENV, "MAB_BGC_FIRST_MONTH", "2018-12") * "-01")
last_month  = DateTime(get(ENV, "MAB_BGC_LAST_MONTH",  "2020-01") * "-01")

for name in (:nitrate, :dissolved_oxygen, :phytoplankton_carbon, :ph), date in first_month:Month(1):last_month
    md = Metadatum(name; dataset = CopernicusBGCMonthly(), region, date, dir)
    isfile(metadata_path(md)) || Downloads.download(md)
    println(name, " ", Dates.format(date, "yyyy-mm"), " ", isfile(metadata_path(md)) ? "ok" : "MISSING", "  ", metadata_path(md))
end
