# Records every anchoring decision of the open-boundary radiation schemes during the first steps of a run:
#   julia -t 1 --project=<Oceananigans worktree> anchored_fill_probe.jl
# Each row: (where, iteration, stage, finite stage Δt?, anchored?) and how many boundary points made that decision.
using Oceananigans
using Oceananigans.BoundaryConditions: NormalRadiation, GravityWaveRadiationBoundaryCondition
using Oceananigans.TimeSteppers: time_step!

const decisions = Dict{Any, Int}()
const phase = Ref("")

@eval Oceananigans.BoundaryConditions begin
    @inline function anchored_fill(clock)
        anchored = clock.stage ≤ 1
        key = ($(phase)[], clock.iteration, clock.stage, isfinite(clock.last_stage_Δt), anchored)
        $(decisions)[key] = get($(decisions), key, 0) + 1
        return anchored
    end
end

function build_model(timestepper)
    grid = RectilinearGrid(CPU(), size=(16, 8, 4), x=(0, 16000), y=(0, 8000), z=(-100, 0),
                           topology=(Bounded, Periodic, Bounded))
    scheme = NormalRadiation(inflow_timescale=100.0, outflow_timescale=1000.0)
    u_bcs = FieldBoundaryConditions(east = NormalFlowBoundaryCondition(0; scheme), west = NormalFlowBoundaryCondition(0; scheme))
    b_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(0.01; scheme), west = ValueBoundaryCondition(0.01; scheme))
    U_bcs = FieldBoundaryConditions(grid, (Face(), Center(), nothing),
                                    east = GravityWaveRadiationBoundaryCondition((0.0, 0.0)),
                                    west = GravityWaveRadiationBoundaryCondition((0.0, 0.0)))
    model = HydrostaticFreeSurfaceModel(grid; buoyancy = BuoyancyTracer(), tracers = :b, timestepper,
                                        free_surface = SplitExplicitFreeSurface(grid; substeps=4),
                                        boundary_conditions = (u=u_bcs, b=b_bcs, U=U_bcs))
    phase[] = "set!"
    set!(model, u = (x, y, z) -> 0.1 * cos(2π * x / 16000), b = (x, y, z) -> 0.005 * (1 + sin(2π * x / 16000)))
    return model
end

for timestepper in (:SplitRungeKutta3, :QuasiAdamsBashforth2)
    empty!(decisions)
    phase[] = "build"
    model = build_model(timestepper)
    simulation = Simulation(model, Δt=10.0, stop_iteration=2, verbose=false)
    phase[] = "run! (initialize + 2 steps)"
    run!(simulation)
    println("\n==== $timestepper   (phase, iteration, stage, finite stage Δt, anchored) => boundary points")
    for key in sort(collect(keys(decisions)), by = k -> (k[2], k[1], k[3]))
        println("  ", key, " => ", decisions[key])
    end
end
