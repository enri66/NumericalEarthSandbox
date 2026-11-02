# Bathymetry as the mean of a dataset's points in each cell of a LatitudeLongitudeGrid no finer than the dataset.
#
# `regrid_bathymetry(...; interpolation_passes)` interpolates through intermediate grids that span the dataset's whole
# (global) domain, with cell counts stepping linearly from the dataset's to the target's; for a regional target the last
# of them is a coarse global grid, so the 1/60° MAB bathymetry came out with about 1/8° resolution and the 1/12° one with
# about 1/6°. Here every cell is a plain average of all the dataset points inside it: at the dataset's own resolution that
# is the point itself, on coarser grids an area mean, and nothing finer than the grid is lost.
using NumericalEarth.DataWrangling: Metadatum, metadata_path, dataset_variable_name

const BATHYMETRY_NC = NumericalEarth.DataWrangling.NCDatasets

"""
    dataset_cell_statistics(grid, dataset; height_above_water = 1)

Mean height of `dataset`'s points in each cell of `grid`, with the points above sea level counted at `height_above_water`
(as `regrid_bathymetry` does), and the fraction of the points above sea level: two (Nx, Ny) arrays. The dataset file must
already be on disk (download it once, serially, before an MPI run).
"""
function dataset_cell_statistics(grid, dataset; height_above_water = 1)
    metadatum = Metadatum(:bottom_height; dataset)
    path = metadata_path(metadatum)
    isfile(path) || error("$path is missing: download $(nameof(typeof(dataset))) first")
    λf = Float64.(collect(Oceananigans.Grids.λnodes(grid, Oceananigans.Face())))
    φf = Float64.(collect(Oceananigans.Grids.φnodes(grid, Oceananigans.Face())))
    Nx, Ny = size(grid, 1), size(grid, 2)
    ds = BATHYMETRY_NC.Dataset(path)
    lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:])
    ilon = findall(x -> λf[1] <= x < λf[end], lon); ilat = findall(y -> φf[1] <= y < φf[end], lat)
    z = ds[dataset_variable_name(metadatum)][ilon[1]:ilon[end], ilat[1]:ilat[end]]
    close(ds)
    total = zeros(Int, Nx, Ny); land = zeros(Int, Nx, Ny); height = zeros(Nx, Ny)
    for (b, y) in enumerate(lat[ilat]), (a, x) in enumerate(lon[ilon])
        v = z[a, b]
        ismissing(v) && continue
        i = searchsortedlast(λf, x); j = searchsortedlast(φf, y)
        total[i, j] += 1
        if v > 0
            land[i, j] += 1; height[i, j] += height_above_water
        else
            height[i, j] += v
        end
    end
    any(==(0), total) && error("some cells of the grid hold no point of $(nameof(typeof(dataset))): the grid is finer than the dataset")
    return height ./ total, land ./ total
end

"""
    cell_mean_bathymetry(grid, dataset; height_above_water = 1, minimum_depth = 10, major_basins = 1)

A `Field{Center, Center, Nothing}` on `grid` with the mean height of `dataset`'s points in each cell, treated otherwise as
`regrid_bathymetry` treats its result: wet cells are at least `minimum_depth` deep and only the `major_basins` largest
basins stay wet.
"""
function cell_mean_bathymetry(grid, dataset; height_above_water = 1, minimum_depth = 10, major_basins = 1)
    h, _ = dataset_cell_statistics(grid, dataset; height_above_water)
    h = @. ifelse(h < 0, min(h, -minimum_depth), h)
    bottom_height = Field{Center, Center, Nothing}(grid)
    set!(bottom_height, h)
    NumericalEarth.Bathymetry.remove_minor_basins!(bottom_height, major_basins)
    Oceananigans.fill_halo_regions!(bottom_height)
    return bottom_height
end
