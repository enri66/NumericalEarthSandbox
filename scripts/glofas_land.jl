# GloFAS river forcing with hand-placed mouths, for rivers whose mouth NumericalEarth's automatic detection misses.
#
# `GloFASPrescribedLand` finds river mouths as finite GloFAS cells that touch an ocean (NaN) cell. For a river that enters a
# long estuary the GloFAS cells below the last big value are small or NaN, so the discharge is never routed: on the MAB box
# the Hudson (about 290 m³/s in September 2019, at 73.93°W 40.88°N) and the Delaware (about 180 m³/s, at 75.23°W 39.88°N) are
# missing, more than the Connecticut River's 308 m³/s. `glofas_land_with_mouths` adds, for each entry of `extra_mouths`, the
# discharge of the GloFAS cell nearest (glofas_λ, glofas_φ), delivered at (mouth_λ, mouth_φ) on the model grid.
#
#   extra_mouths = [(name = "Hudson", glofas_λ = -73.925, glofas_φ = 40.875, mouth_λ = -73.95, mouth_φ = 40.50), ...]
using Oceananigans
using Oceananigans.Architectures: architecture
using Oceananigans.Fields: Field
using Oceananigans.Grids: λnodes, φnodes, Center
using NumericalEarth.DataWrangling: Metadata
using NumericalEarth.Lands: PrescribedLand, build_river_routing, coastal_outlet_indices, routable_grid
using Printf

const MAB_EXTRA_MOUTHS = [(name = "Hudson",   glofas_λ = -73.925, glofas_φ = 40.875, mouth_λ = -73.95, mouth_φ = 40.50),
                          (name = "Delaware", glofas_λ = -75.225, glofas_φ = 39.875, mouth_λ = -75.30, mouth_φ = 39.20)]

function glofas_land_with_mouths(grid; extra_mouths = MAB_EXTRA_MOUTHS, dataset = GloFASReanalysis(), start_date, end_date, dir, region,
                                 time_indices_in_memory = 10, time_indexing = Oceananigans.OutputReaders.Cyclical(),
                                 freshwater_density = 1000, maximum_search_radius = 5, spread_radius = 1.2, maximum_spread_cells = 8, say = println)
    arch = architecture(grid)
    discharge_meta = Metadata(:river_discharge; dataset, start_date, end_date, dir, region)
    discharge = FieldTimeSeries(discharge_meta, arch; time_indexing, time_indices_in_memory)

    snapshot = Field(first(discharge_meta), arch)
    outlet_i, outlet_j, outlet_λ, outlet_φ = coastal_outlet_indices(snapshot)
    outlet_i = collect(outlet_i); outlet_j = collect(outlet_j); outlet_λ = collect(outlet_λ); outlet_φ = collect(outlet_φ)
    say("rivers: $(length(outlet_i)) mouths found from the GloFAS coastline")

    gλ = collect(λnodes(snapshot.grid, Center())); gφ = collect(φnodes(snapshot.grid, Center()))
    for m in extra_mouths
        i = argmin(abs.(gλ .- m.glofas_λ)); j = argmin(abs.(gφ .- m.glofas_φ))
        q = snapshot[i, j, 1]
        isfinite(q) || error("river $(m.name): the GloFAS cell at $(gλ[i]), $(gφ[j]) has no discharge")
        push!(outlet_i, i); push!(outlet_j, j); push!(outlet_λ, m.mouth_λ); push!(outlet_φ, m.mouth_φ)
        say(@sprintf("rivers: added %s, GloFAS cell %.3f°E %.3f°N (%.0f m³/s on the first day) delivered at %.3f°E %.3f°N",
                     m.name, gλ[i], gφ[j], q, m.mouth_λ, m.mouth_φ))
    end

    outlet_weight = fill(convert(eltype(grid), freshwater_density), length(outlet_i))
    river_routing = routable_grid(grid) ?
        (; rivers = build_river_routing(grid, outlet_i, outlet_j, outlet_λ, outlet_φ, outlet_weight;
                                        maximum_search_radius, spread_radius, maximum_spread_cells)) : nothing
    return PrescribedLand((; rivers = discharge); river_routing)
end
