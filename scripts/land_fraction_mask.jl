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
# Depths of the cells that stay wet are untouched. The land fraction can come from another dataset (`dataset`), e.g. GEBCO.
using Printf
@isdefined(dataset_cell_statistics) || include(joinpath(@__DIR__, "cell_mean_bathymetry.jl"))

"Fraction (0-1) of `dataset`'s points above sea level in each cell of `grid`; (Nx, Ny)."
etopo_land_fraction(grid; dataset = ETOPO2022()) = dataset_cell_statistics(grid, dataset)[2]

"""
    apply_land_fraction!(bottom_height, grid, threshold; dataset = ETOPO2022(), minimum_depth = 10, major_basins = 1, say = println)

Replace the land/wet decision of `bottom_height` (a `Field{Center, Center, Nothing}` on `grid`, wet where negative) by
the land fraction of `dataset` > `threshold`, as described at the top of this file. Reports what changed.
"""
function apply_land_fraction!(bottom_height, grid, threshold; dataset = ETOPO2022(), minimum_depth = 10, major_basins = 1, say = println)
    f = etopo_land_fraction(grid; dataset)
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
