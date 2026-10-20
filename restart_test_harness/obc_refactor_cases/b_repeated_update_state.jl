using Oceananigans, Test
XT = get(ENV, "OPEN", "") == "true" ? Open : Bounded
YT = XT
using Oceananigans.TimeSteppers: update_state!
using Oceananigans.BoundaryConditions: NormalRadiation, ObliqueRadiation, TracerReservoir, GravityWaveRadiationBoundaryCondition

# A second update_state! in the same iteration, as after restoring a checkpoint, leaves the tracer and its boundary unchanged
function repeated_update_leaves_unchanged(velocity_scheme, tracer_scheme, timestepper)
    grid = RectilinearGrid(size=(16, 8, 4), x=(0, 16000), y=(0, 8000), z=(-100, 0), topology=(XT, Periodic, Bounded))
    u_bcs = FieldBoundaryConditions(east = NormalFlowBoundaryCondition(0; scheme=velocity_scheme),
                                    west = NormalFlowBoundaryCondition(0; scheme=velocity_scheme))
    b_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(0.01; scheme=tracer_scheme),
                                    west = ValueBoundaryCondition(0.01; scheme=tracer_scheme))
    U_bcs = FieldBoundaryConditions(grid, (Face(), Center(), nothing);
                                    east = GravityWaveRadiationBoundaryCondition((0.0, 0.0)),
                                    west = GravityWaveRadiationBoundaryCondition((0.0, 0.0)))
    model = HydrostaticFreeSurfaceModel(grid; timestepper, buoyancy=BuoyancyTracer(), tracers=:b,
                                        free_surface=SplitExplicitFreeSurface(grid; substeps=10),
                                        boundary_conditions=(u=u_bcs, U=U_bcs, b=b_bcs))
    set!(model, u = (x, y, z) -> 0.1 * cos(2π * x / 16000), b = (x, y, z) -> 0.005 * (1 + sin(2π * x / 16000)))
    for _ in 1:3
        time_step!(model, 10)
    end
    b = Array(parent(model.tracers.b)); u = Array(parent(model.velocities.u))
    update_state!(model)
    return Array(parent(model.tracers.b)) == b && Array(parent(model.velocities.u)) == u
end

@testset "A repeated update_state! leaves the open boundary unchanged [$timestepper]" for timestepper in (:QuasiAdamsBashforth2, :SplitRungeKutta3)
    normal = NormalRadiation(inflow_timescale=100, outflow_timescale=1000)
    @test repeated_update_leaves_unchanged(normal, normal, timestepper)
    @test repeated_update_leaves_unchanged(ObliqueRadiation(inflow_timescale=100, outflow_timescale=1000),
                                           TracerReservoir(inflow_length_scale=5000, outflow_length_scale=2000), timestepper)
end
