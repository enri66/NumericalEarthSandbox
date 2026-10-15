# Download the time-independent data script 05 reads for a region, with NumericalEarth's own downloaders: GLORYS's
# static depth (`deptho`, for match_boundary_bathymetry!) and ERA5's land-sea mask (for era5_land_flooding.jl).
# The time-dependent fields come from the global archives (regional_cache_from_archives.jl).
#
#   STAGE_LON="-71.5,-48" STAGE_LAT="45,56" STAGE_DIR=~/Data/NumericalEarth_GSL julia --project=. stage_static_region_data.jl
using NumericalEarth, Dates
using CopernicusMarine              # activates the GLORYS download backend
using CopernicusClimateDataStore    # activates the ERA5 download backend
using NumericalEarth.DataWrangling: Metadata, Metadatum, BoundingBox, metadata_path
using Downloads: download
using NumericalEarth.DataWrangling.GLORYS: GLORYSStatic
using NumericalEarth.DataWrangling.ERA5: ERA5HourlySingleLevel
include(joinpath(@__DIR__, "era5_land_flooding.jl"))

parse_pair(s) = Tuple(parse.(Float64, split(s, ",")))
region = BoundingBox(longitude = parse_pair(ENV["STAGE_LON"]), latitude = parse_pair(ENV["STAGE_LAT"]))
dir = ENV["STAGE_DIR"]

depth = Metadata(:depth; dataset = GLORYSStatic(), dir, region)
download(depth)
println("GLORYS static depth: ", metadata_path(depth))

lsm = Metadatum(:land_sea_mask; dataset = ERA5HourlySingleLevel(), region, dir = joinpath(dir, "era5"),
                date = DateTime(2019, 1, 1))
download(lsm)
println("ERA5 land-sea mask: ", metadata_path(lsm))
