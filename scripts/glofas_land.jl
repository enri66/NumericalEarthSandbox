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
using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.Fields: Field
using Oceananigans.Grids: λnodes, φnodes, Center
using NumericalEarth.DataWrangling: Metadata
using NumericalEarth.Lands: PrescribedLand, RiverRouting, build_river_routing, coastal_outlet_indices, routable_grid
using Printf

const MAB_EXTRA_MOUTHS = [(name = "Hudson",   glofas_λ = -73.925, glofas_φ = 40.875, mouth_λ = -73.95, mouth_φ = 40.50),
                          (name = "Delaware", glofas_λ = -75.225, glofas_φ = 39.875, mouth_λ = -75.30, mouth_φ = 39.20)]

"""
    localize_routing(routing, i_offset, j_offset, nx, ny, arch)

Keep the destination cells of a whole-domain `RiverRouting` that fall in a rank's block (whole-domain indices
`i_offset+1 : i_offset+nx`, `j_offset+1 : j_offset+ny`) and shift them to the rank's local indices. NumericalEarth's own
distributed path rebuilds the immersed boundary on the CPU and fails with halo-fill errors on our grids, so the routing is
built once on the whole-domain grid and cut up here.
"""
function localize_routing(routing, i_offset, j_offset, nx, ny, arch)
    ti = Array(routing.target_i); tj = Array(routing.target_j); offs = Array(routing.offsets)
    coi = Array(routing.contribution_outlet_i); coj = Array(routing.contribution_outlet_j); w = Array(routing.contribution_weight)
    target_i = Int[]; target_j = Int[]; offsets = Int[1]; o_i = Int[]; o_j = Int[]; weight = eltype(w)[]
    for c in eachindex(ti)
        (i_offset < ti[c] <= i_offset + nx && j_offset < tj[c] <= j_offset + ny) || continue
        push!(target_i, ti[c] - i_offset); push!(target_j, tj[c] - j_offset)
        for k in offs[c]:offs[c+1]-1
            push!(o_i, coi[k]); push!(o_j, coj[k]); push!(weight, w[k])
        end
        push!(offsets, length(o_i) + 1)
    end
    return RiverRouting(on_architecture(arch, o_i), on_architecture(arch, o_j), on_architecture(arch, weight),
                        on_architecture(arch, target_i), on_architecture(arch, target_j), on_architecture(arch, offsets))
end

function glofas_land_with_mouths(grid; extra_mouths = MAB_EXTRA_MOUTHS, dataset = GloFASReanalysis(), start_date, end_date, dir, region,
                                 time_indices_in_memory = 10, time_indexing = Oceananigans.OutputReaders.Cyclical(),
                                 freshwater_density = 1000, maximum_search_radius = 5, spread_radius = 1.2, maximum_spread_cells = 8, say = println,
                                 routing_grid = grid, block = nothing)
    # The GloFAS data are global-window fields, not partitioned over ranks (as for ERA5): build them on the CPU, and only the
    # routing onto `grid` (which may be distributed) uses the model grid.
    arch = CPU()
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
    # `routing_grid` is the grid the mouths are routed onto: `grid` itself, or on a distributed run the whole-domain grid on
    # the CPU with `block = (i_offset, j_offset, nx, ny)` of this rank's part, which is then cut out of the whole-domain routing
    rivers = build_river_routing(routing_grid, outlet_i, outlet_j, outlet_λ, outlet_φ, outlet_weight;
                                 maximum_search_radius, spread_radius, maximum_spread_cells)
    isnothing(block) || (rivers = localize_routing(rivers, block..., architecture(grid)))
    river_routing = routable_grid(routing_grid) ? (; rivers) : nothing
    return PrescribedLand((; rivers = discharge); river_routing)
end
