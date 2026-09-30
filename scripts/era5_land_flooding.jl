# Flood ERA5 from the ocean over land, so that coastal ocean cells never interpolate land values (surface pressure
# over the Appalachians, land air temperatures, winds and radiation).
#
# NumericalEarth fills the masked cells of a dataset field when it reads a file (`inpaint_mask!`), on the dataset's
# own grid and before the model interpolates it, and caches the result next to the file (`*_inpainted.jld2`). For
# ERA5 that is switched off, and its default mask (missing values) would not mark land anyway, because ERA5 has
# values there. Including this file marks every ERA5 point whose land-sea mask exceeds ERA5_LAND_THRESHOLD (default 0,
# i.e. any land) and fills it from the surrounding ocean points, whenever ERA5 is read with an `inpainting` method,
# e.g. `ERA5PrescribedAtmosphere(...; inpainting = ERA5_FLOODING, cache_inpainted_data = false)`.
#
# If the filled fields are cached (`cache_inpainted_data = true`), delete the ERA5 `*_inpainted.jld2` files after
# changing the threshold: the cache check compares only the inpainting method.
using NumericalEarth.DataWrangling: DataWrangling, Metadatum, NearestNeighborInpainting
using NumericalEarth.DataWrangling.ERA5: ERA5Metadatum, ERA5_dataset_variable_names, ERA5_netcdf_variable_names
using Oceananigans.Architectures: architecture
using Oceananigans.Fields: location
using Dates

# ERA5's land-sea mask (land fraction, 0 to 1) is not among NumericalEarth's ERA5 variables yet
ERA5_dataset_variable_names[:land_sea_mask] = "land_sea_mask"
ERA5_netcdf_variable_names[:land_sea_mask]  = "lsm"

const ERA5_LAND_THRESHOLD = parse(Float64, get(ENV, "ERA5_LAND_THRESHOLD", "0"))
const ERA5_FLOODING = NearestNeighborInpainting(Inf)

# The mask does not change in time, so one hour of it is read once per region and directory
const era5_land_masks = Dict{Any, Any}()

function era5_land(md::ERA5Metadatum, grid)
    get!(era5_land_masks, (md.dataset, md.region, md.dir)) do
        lsm_md = Metadatum(:land_sea_mask; dataset = md.dataset, region = md.region, dir = md.dir,
                           date = DateTime(2019, 1, 1))
        lsm = Field(lsm_md, architecture(grid); inpainting = nothing)
        Array(interior(lsm)) .> ERA5_LAND_THRESHOLD
    end
end

function DataWrangling.compute_mask(md::ERA5Metadatum, field, args...)
    md.name === :land_sea_mask && return nothing
    LX, LY, LZ = location(field)
    mask = Field{LX, LY, LZ}(field.grid, Bool)
    interior(mask) .= era5_land(md, field.grid)
    return mask
end

# NumericalEarth reports every filled field ("Inpainting ERA5... data from <date>..." and " ... (<time>)"), which is two
# log lines per variable per hour of ERA5; leave those out of the log and pass everything else through.
using Logging

struct SkipERA5InpaintingMessages{L} <: Logging.AbstractLogger
    logger :: L
end

Logging.min_enabled_level(l::SkipERA5InpaintingMessages) = Logging.min_enabled_level(l.logger)
Logging.shouldlog(l::SkipERA5InpaintingMessages, args...) = Logging.shouldlog(l.logger, args...)
Logging.catch_exceptions(l::SkipERA5InpaintingMessages) = Logging.catch_exceptions(l.logger)

function Logging.handle_message(l::SkipERA5InpaintingMessages, level, message, args...; kw...)
    skip = level == Logging.Info && message isa AbstractString &&
           (startswith(message, "Inpainting ERA5") || startswith(message, " ... ("))
    skip || Logging.handle_message(l.logger, level, message, args...; kw...)
    return nothing
end

global_logger(SkipERA5InpaintingMessages(global_logger()))
