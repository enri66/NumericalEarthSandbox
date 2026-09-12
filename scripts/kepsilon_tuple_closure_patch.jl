# ==================================================================
# WORKAROUND for a gap in Oceananigans (as of the RPT5V pin): `top_tke_flux`
# has explicit methods unwrapping a Tuple closure (`tke_top_boundary_condition.jl`,
# a generic zero fallback plus a sum over 1/2/3-element tuples), but the matching
# `top_dissipation_flux` for TKEDissipationVerticalDiffusivity's `:ϵ` equation
# (`tke_dissipation_equations.jl`) has neither — only a single method for a bare
# `TKEDissipationVerticalDiffusivity`. Since this script always pairs the vertical
# closure with `HorizontalScalarDiffusivity` in a Tuple, k-ϵ hits a MethodError the
# moment it is combined with anything else, on any grid.
#
# These four methods mirror `top_tke_flux`'s pattern exactly, so `top_dissipation_flux`
# gains the same Tuple support. Isolated in this file rather than the package: it is a
# workaround for testing, not a claim that this is the right home for the fix.
# ==================================================================

using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities

const TDVD = Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities

@inline TDVD.top_dissipation_flux(i, j, grid, clock, fields, parameters, closure, buoyancy) = zero(grid)

@inline TDVD.top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple::Tuple{<:Any}, buoyancy) =
    TDVD.top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[1], buoyancy)

@inline TDVD.top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple::Tuple{<:Any, <:Any}, buoyancy) =
    TDVD.top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[1], buoyancy) +
    TDVD.top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[2], buoyancy)

@inline TDVD.top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple::Tuple{<:Any, <:Any, <:Any}, buoyancy) =
    TDVD.top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[1], buoyancy) +
    TDVD.top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[2], buoyancy) +
    TDVD.top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[3], buoyancy)
