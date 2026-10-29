# Stepping of OceanBioME's sediment (SimpleMultiGSediment) under Oceananigans' SplitRungeKuttaTimeStepper, included by mab_bgc.jl
# with MAB_BGC_SEDIMENT_STEPPING=once (default "stages", the stepping OceanBioME does today).
#
# OceanBioME steps the sediment inside update_biogeochemical_state!, once per stage. It supports AB2 and RungeKutta3TimeStepper;
# for any other stepper (so for the SplitRungeKutta that Oceananigans 0.114 uses) it falls back to time_step!(sediment, last_stage_Δt)
# at every stage, which ADDS the stage's Δt to a state that carries on from the previous stage. The split scheme restarts every
# stage from the state at the start of the step Uⁿ: stage m is Uⁿ + Δt/βᵐ G, with β = (3, 2, 1), so the water column advances Δt per
# step and the sediment Δt/3 + Δt/2 + Δt = 11/6 Δt: its pools rise and decay 1.83 times too fast (the equilibria are not changed,
# the spin-up is). This method does what the water column does: it keeps the sediment at the start of the step and sets, at
# stage m, S = Sⁿ + Δt/βᵐ G, with G the sediment tendency computed at the end of the previous stage (from the sediment and the
# tracked water-column fields of the state that stage m steps from).
using Oceananigans.Fields: Field, Center
using Oceananigans.TimeSteppers: SplitRungeKuttaTimeStepper, rk3_substep!
using Oceananigans.Models: prognostic_fields
import OceanBioME.Sediments: step_sediment!

mutable struct SedimentStepStart{F}
    start :: F                  # the sediment fields at the start of the step
    iteration :: Int            # the iteration they were taken at (update_state! can be called twice in one stage)
end
const SEDIMENT_STEP_START = IdDict{Any, SedimentStepStart}()

function step_sediment!(sediment, model, ::SplitRungeKuttaTimeStepper)
    fields = prognostic_fields(sediment)
    s = get!(SEDIMENT_STEP_START, sediment) do
        SedimentStepStart(map(f -> Field{Center, Center, Nothing}(sediment.grid), fields), -1)
    end

    # the start of the step: taken at the first call of the iteration (the stage-1 update at iteration 0 is skipped by OceanBioME,
    # which leaves the sediment as it was)
    if s.iteration != model.clock.iteration
        for name in keys(fields)
            parent(s.start[name]) .= parent(fields[name])
        end
        s.iteration = model.clock.iteration
    end

    # S = Sⁿ + Δt/β G: back to the start of the step, then OceanBioME's own increment U += Δτ G (columns without sea floor untouched)
    for name in keys(fields)
        parent(fields[name]) .= parent(s.start[name])
    end
    rk3_substep!(sediment, model.clock.last_stage_Δt, 1.0, nothing, nothing)
    sediment.clock.time = model.clock.time

    return nothing
end
