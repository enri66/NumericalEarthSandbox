# Fetch ERA5's land-sea mask for the GSL + Newfoundland-Labrador region once, so parallel runs find it on disk.
using NumericalEarth, CopernicusClimateDataStore, Oceananigans, Dates
using NumericalEarth.DataWrangling: Metadatum, BoundingBox
include(joinpath(pwd(), "scripts", "era5_land_flooding.jl"))
md = Metadatum(:land_sea_mask; dataset = ERA5HourlySingleLevel(), date = DateTime(2019, 1, 1), dir = joinpath(homedir(), "Data", "NumericalEarth", "era5"),
               region = BoundingBox(longitude = (-71.5, -48.0), latitude = (45.0, 56.0)))
lsm = Field(md; inpainting = nothing)
println("land-sea mask: ", size(interior(lsm)), ", land points ", count(>(0), interior(lsm)))
