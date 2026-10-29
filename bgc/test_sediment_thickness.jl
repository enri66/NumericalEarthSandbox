# How much of the nitrogen the multi-G sediment returns reaches the water on a grid with partial bottom cells? A 4×4×20 column grid
# with PartialCellBottom (bottom heights between 21 and 38 m, the 2 m cells so the bottom cell is 0.1-1 of a full cell), LOBSTER +
# SimpleMultiGSediment, a few steps so the pools fill; then the tendencies are zeroed and only the sediment's coupling is applied:
# the nitrogen the bottom cells receive (Σ (G_NH₄ + G_NO₃) Δzᶜᶜᶜ Az) against what the sediment decays (Σ (λs Ns + λf Nf) Az),
# without and with bgc/sediment_thickness.jl. Usage: julia --project=<env> bgc/test_sediment_thickness.jl
using Oceananigans, OceanBioME, Printf
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, PartialCellBottom
using Oceananigans.Operators: Δzᶜᶜᶜ, Azᶜᶜᶜ
using Oceananigans.Biogeochemistry: update_tendencies!
include(joinpath(@__DIR__, "sediment_thickness.jl"))

underlying = RectilinearGrid(size = (4, 4, 20), x = (0, 4000), y = (0, 4000), z = (-40, 0), topology = (Periodic, Periodic, Bounded))
bottom(x, y) = -(21 + 17 * (x + 1.7y) / 7000)
grid = ImmersedBoundaryGrid(underlying, PartialCellBottom(bottom))
sediment = SimpleMultiGSediment(grid; sinking_nitrogen = (:sPOM, :bPOM))
biogeochemistry = LOBSTER(grid; oxygen = Oxygen(), inorganic_carbon = CarbonateSystem(), sediment)
model = HydrostaticFreeSurfaceModel(grid; biogeochemistry, timestepper = :SplitRungeKutta3, buoyancy = nothing, tracers = (),
                                    momentum_advection = nothing)
set!(model, P = 0.1, Z = 0.1, NO₃ = 5.0, NH₄ = 0.5, sPOM = 1.0, bPOM = 1.0, DOM = 0.1, O₂ = 300.0, DIC = 2000.0, Alk = 2300.0)
for _ in 1:20; time_step!(model, 60.0); end

function delivered(model, sediment, grid)
    Gⁿ = model.timestepper.Gⁿ
    for name in keys(Gⁿ); parent(Gⁿ[name]) .= 0; end
    update_tendencies!(model.biogeochemistry, sediment, model)
    return Gⁿ
end

λs, λf = sediment.biogeochemistry.slow_decay_rate, sediment.biogeochemistry.fast_decay_rate
kb = sediment.bottom_indices
Ns, Nf = sediment.fields.Ns, sediment.fields.Nf
function totals(Gⁿ)
    got = 0.0; decayed = 0.0; thin = Float64[]
    for j in 1:4, i in 1:4
        k = kb[i, j, 1]
        got += (Gⁿ.NH₄[i, j, k] + Gⁿ.NO₃[i, j, k]) * Δzᶜᶜᶜ(i, j, k, grid) * Azᶜᶜᶜ(i, j, k, grid)
        decayed += (λs * Ns[i, j, 1] + λf * Nf[i, j, 1]) * Azᶜᶜᶜ(i, j, k, grid)
        push!(thin, Δzᶜᶜᶜ(i, j, k, grid) / 2)
    end
    return got, decayed, thin
end

Gⁿ = delivered(model, sediment, grid)
got, decayed, thin = totals(Gⁿ)
@printf("OceanBioME as it is:   N received by the bottom cells / N decayed in the sediment = %.4f   (bottom cell / full cell: %.2f - %.2f)\n", got / decayed, minimum(thin), maximum(thin))

update_tendencies!(model.biogeochemistry, SedimentThickness(), model)
got, decayed, _ = totals(Gⁿ)
@printf("with SedimentThickness: N received by the bottom cells / N decayed in the sediment = %.4f\n", got / decayed)
