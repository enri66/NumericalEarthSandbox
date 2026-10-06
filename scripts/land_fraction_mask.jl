# Land mask from the ETOPO land fraction, for regridded bathymetry on a LatitudeLongitudeGrid.
#
# `regrid_bathymetry(...; height_above_water = 1, interpolation_passes = 10)` sets land to +1 m, averages and smooths,
# and calls a cell land only if the mean height is >= 0. For 20 m of water that needs a land fraction of d/(1 + d) =
# 95%, so narrow land disappears: on the 1/12° MAB grid almost all of Long Island east of 73.2°W comes out as water
# 10-15 m deep (see the land-mask print in the 2026-10-06 notes). This file recomputes the land mask from the fraction
# of ETOPO2022's 1/60° points above sea level in each model cell (exactly 5 × 5 points per 1/12° cell, edges aligned):
#   a cell is land when its land fraction exceeds `threshold` (default use 0.5);
#   cells that regrid_bathymetry made land but whose fraction is at most `threshold` become wet at the minimum depth;
#   cells that become land are set to +1 m (the height_above_water value), and `remove_minor_basins!` keeps the
#   largest basin as before.
# Depths of the cells that stay wet are untouched (smoothed, with land mixed in at +1 m, as before).
using NumericalEarth.DataWrangling: Metadatum, metadata_path
using Printf

const ETOPO_NC = NumericalEarth.DataWrangling.NCDatasets

"Fraction (0-1) of ETOPO2022's points above sea level in each cell of `grid`; (Nx, Ny)."
function etopo_land_fraction(grid)
    path = metadata_path(Metadatum(:bottom_height; dataset = ETOPO2022()))
    ds = ETOPO_NC.Dataset(path)
    lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:])
    λf = Float64.(collect(Oceananigans.Grids.λnodes(grid, Oceananigans.Face())))
    φf = Float64.(collect(Oceananigans.Grids.φnodes(grid, Oceananigans.Face())))
    Nx, Ny = size(grid, 1), size(grid, 2)
    ilon = findall(x -> λf[1] <= x <= λf[end], lon); ilat = findall(y -> φf[1] <= y <= φf[end], lat)
    z = ds["z"][ilon, ilat]
    close(ds)
    Δλ = (λf[end] - λf[1]) / Nx; Δφ = (φf[end] - φf[1]) / Ny
    land = zeros(Int, Nx, Ny); total = zeros(Int, Nx, Ny)
    for (b, y) in enumerate(lat[ilat]), (a, x) in enumerate(lon[ilon])
        i = clamp(floor(Int, (x - λf[1]) / Δλ) + 1, 1, Nx); j = clamp(floor(Int, (y - φf[1]) / Δφ) + 1, 1, Ny)
        total[i, j] += 1
        (!ismissing(z[a, b]) && z[a, b] > 0) && (land[i, j] += 1)
    end
    return land ./ max.(total, 1)
end

"""
    apply_land_fraction!(bottom_height, grid, threshold; minimum_depth = 10, major_basins = 1, say = println)

Replace the land/wet decision of `bottom_height` (a `Field{Center, Center, Nothing}` on `grid`, wet where negative) by
the ETOPO land fraction > `threshold`, as described at the top of this file. Reports what changed.
"""
function apply_land_fraction!(bottom_height, grid, threshold; minimum_depth = 10, major_basins = 1, say = println)
    f = etopo_land_fraction(grid)
    h = Array(interior(bottom_height))[:, :, 1]
    was_land = h .>= 0
    now_land = f .> threshold
    to_land = now_land .& .!was_land
    to_wet  = .!now_land .& was_land
    h[to_land] .= 1.0
    h[to_wet]  .= -minimum_depth
    set!(bottom_height, h)
    NumericalEarth.Bathymetry.remove_minor_basins!(bottom_height, major_basins)
    h2 = Array(interior(bottom_height))[:, :, 1]
    removed = count(@. (h2 >= 0) & (h < 0))
    Oceananigans.fill_halo_regions!(bottom_height)
    say(@sprintf("land fraction mask (> %.2f): %d cells wet → land, %d land → wet (at %g m), %d more removed as minor basins; wet cells %d → %d",
                 threshold, count(to_land), count(to_wet), minimum_depth, removed, count(.!was_land), count(h2 .< 0)))
    return bottom_height
end
