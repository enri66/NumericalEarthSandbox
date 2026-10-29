# Does the sediment get one time step of input per step? A small column with LOBSTER + SimpleMultiGSediment under the split
# Runge-Kutta (β = 3, 2, 1), a few steps, with MAB_BGC_SEDIMENT_STEPPING=stages (OceanBioME) and =once (bgc/sediment_stepping.jl).
# The refractory pool Nr has no loss, so Nr = (refractory fraction) × (the sinking flux of nitrogen to the bottom) × Δt summed over
# the steps; the ratio of the two is printed (expected 11/6 = 1.833 for stages, 1 for once). Usage:
#   MAB_BGC_SEDIMENT_STEPPING=once julia --project=<env> bgc/test_sediment_stepping.jl
using Oceananigans, OceanBioME, Printf
using Oceananigans.Fields: interior
using Oceananigans.Models: prognostic_fields

const MODE = get(ENV, "MAB_BGC_SEDIMENT_STEPPING", "stages")
MODE == "once" && include(joinpath(@__DIR__, "sediment_stepping.jl"))

grid = RectilinearGrid(size = (2, 2, 20), x = (0, 2000), y = (0, 2000), z = (-40, 0), topology = (Periodic, Periodic, Bounded))
sediment = SimpleMultiGSediment(grid; sinking_nitrogen = (:sPOM, :bPOM))
biogeochemistry = LOBSTER(grid; oxygen = Oxygen(), inorganic_carbon = CarbonateSystem(), sediment)
model = HydrostaticFreeSurfaceModel(grid; biogeochemistry, timestepper = :SplitRungeKutta3, buoyancy = nothing, tracers = (),
                                    momentum_advection = nothing)
set!(model, P = 0.1, Z = 0.1, NO₃ = 5.0, NH₄ = 0.5, sPOM = 1.0, bPOM = 1.0, DOM = 0.1, O₂ = 300.0, DIC = 2000.0, Alk = 2300.0)

Δt = 60.0
nsteps = 30
Nr = sediment.fields.Nr
flux = zeros(nsteps)
for n in 1:nsteps
    time_step!(model, Δt)
    # the nitrogen flux the sediment receives: sum of the tracked sinking fluxes (mmol N/m²/s) of the first column, at this state
    tf = sediment.tracked_fields
    flux[n] = sum(Array(interior(tf[name]))[1, 1, 1] for name in keys(tf) if occursin("sPOM", String(name)) || occursin("bPOM", String(name)))
end
refractory_fraction = sediment.biogeochemistry.refactory_fraction
Nr1 = Array(interior(Nr))[1, 1, 1]
expected = refractory_fraction * sum(flux) * Δt
@printf("stepping %-6s  Nr = %.6g mmol N/m²   refractory fraction × Σ flux Δt = %.6g   ratio %.4f   (tracked fluxes: %s)\n",
        MODE, Nr1, expected, Nr1 / expected, join(String.(keys(sediment.tracked_fields)), ", "))
