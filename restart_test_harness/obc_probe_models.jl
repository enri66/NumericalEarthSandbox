using Oceananigans
using Oceananigans.BoundaryConditions: NormalRadiation, ObliqueRadiation, TracerReservoir, GravityWaveRadiationBoundaryCondition

configurations = (("NormalRadiation/NormalRadiation", NormalRadiation(inflow_timescale=100.0, outflow_timescale=1000.0),
                                                      NormalRadiation(inflow_timescale=100.0, outflow_timescale=1000.0)),
                  ("ObliqueRadiation/TracerReservoir", ObliqueRadiation(inflow_timescale=100.0, outflow_timescale=1000.0),
                                                       TracerReservoir(inflow_length_scale=5000.0, outflow_length_scale=2000.0)))

function probe_model(timestepper, velocity_scheme, tracer_scheme)
    grid = RectilinearGrid(CPU(), size=(16, 8, 4), x=(0, 16000), y=(0, 8000), z=(-100, 0), topology=(Bounded, Periodic, Bounded))
    u_bcs = FieldBoundaryConditions(east = NormalFlowBoundaryCondition(0; scheme=velocity_scheme),
                                    west = NormalFlowBoundaryCondition(0; scheme=velocity_scheme))
    b_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(0.01; scheme=tracer_scheme),
                                    west = ValueBoundaryCondition(0.01; scheme=tracer_scheme))
    U_bcs = FieldBoundaryConditions(grid, (Face(), Center(), nothing),
                                    east = GravityWaveRadiationBoundaryCondition((0.0, 0.0)),
                                    west = GravityWaveRadiationBoundaryCondition((0.0, 0.0)))
    model = HydrostaticFreeSurfaceModel(grid; buoyancy = BuoyancyTracer(), tracers = :b, timestepper,
                                        free_surface = SplitExplicitFreeSurface(grid; substeps=10),
                                        boundary_conditions = (u=u_bcs, b=b_bcs, U=U_bcs))
    set!(model, u = (x, y, z) -> 0.1 * cos(2π * x / 16000), b = (x, y, z) -> 0.005 * (1 + sin(2π * x / 16000)))
    return model
end
