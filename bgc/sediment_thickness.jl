# Cell thickness in the sediment-water coupling, included by mab_bgc.jl with MAB_BGC_SEDIMENT_THICKNESS=true.
#
# OceanBioME's update_tendencies!(bgc, ::BiogeochemicalSediment, model) (Sediments/tracer_coupling.jl) adds the sediment's flux of
# each coupled tracer (NH₄, NO₃, O₂ of the multi-G sediment; the receiver of the instant remineralisation) to the bottom wet cell k as
#     G[i, j, k] += flux / Δzᶜᶜᶠ(i, j, k, grid)
# Δzᶜᶜᶠ at k is the centre-to-centre spacing at the bottom FACE of the cell. With a partial bottom cell (MAB_BOTTOM=partial) the
# cell is thinner than the full cell (hb between 0.1 and 1 of it), but Δzᶜᶜᶠ(i, j, k) is the spacing of the FULL cells k-1 and k,
# so the tracer received is flux / Δz_full instead of flux / hb: on average about half of the nitrogen the sediment returns never
# reaches the water (and the oxygen it takes is halved too), while the sinking flux that leaves the water column is the exact
# advective flux through the face. The nitrogen of the shelf then loses ~ (1 - hb / Δz_full) of the return every day.
# This modifier runs after the sediment's own update_tendencies! (OceanBioME calls the modifiers last) and adds the difference,
#     G[i, j, k] += flux * (1 / Δzᶜᶜᶜ - 1 / Δzᶜᶜᶠ),
# so the net is flux / Δzᶜᶜᶜ, the thickness of the cell. Columns of land are skipped.
using Oceananigans.Grids: inactive_cell
using Oceananigans.Operators: Δzᶜᶜᶜ, Δzᶜᶜᶠ
import Oceananigans.Biogeochemistry: update_tendencies!, update_biogeochemical_state!
import OceanBioME.Sediments as OBSediments

struct SedimentThickness end

update_biogeochemical_state!(model, ::SedimentThickness) = nothing

function update_tendencies!(bgc, ::SedimentThickness, model)
    sediment = bgc.sediment
    sediment isa OBSediments.BiogeochemicalSediment || return nothing

    grid = model.grid
    Nx, Ny, _ = size(grid)
    Gⁿ = model.timestepper.Gⁿ
    kb = sediment.bottom_indices

    for name in OBSediments.coupled_tracers(sediment)
        G = Gⁿ[name]
        args = (Val(name), sediment.biogeochemistry, OBSediments.fields(sediment), sediment.tracked_fields, sediment.clock)
        for j in 1:Ny, i in 1:Nx
            k = kb[i, j, 1]
            inactive_cell(i, j, k, grid) && continue
            flux = OBSediments.sediment_tendencies(i, j, grid, args...)
            G[i, j, k] += flux * (1 / Δzᶜᶜᶜ(i, j, k, grid) - 1 / Δzᶜᶜᶠ(i, j, k, grid))
        end
    end

    return nothing
end
