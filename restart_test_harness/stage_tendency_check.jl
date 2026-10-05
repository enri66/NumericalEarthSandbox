# Runs four configurations for 10 steps and serializes every prognostic and closure array, so a code change
# can be checked for bit-identical results:  STC_TAG=before julia --project=fixtest_env stage_tendency_check.jl
using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.BoundaryConditions: ObliqueRadiation, TracerReservoir, GravityWaveRadiationBoundaryCondition
using NumericalEarth
using Serialization, Random

const TAG = get(ENV, "STC_TAG", "run")

function regional_model(timestepper)
    underlying_grid = LatitudeLongitudeGrid(size = (12, 10, 8), halo = (4, 4, 4),
                                            longitude = (-71, -70), latitude = (37, 38), z = (-100, 0))
    bump(λ, φ) = -100 + 60 * exp(-((λ + 70.5)^2 + (φ - 37.5)^2) / 0.15^2)
    grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bump))

    velocity_scheme = ObliqueRadiation(inflow_timescale = 600.0, outflow_timescale = 6000.0)
    tracer_scheme   = TracerReservoir(inflow_length_scale = 5000.0, outflow_length_scale = 2000.0)
    u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(-1e-4),
                                    east = NormalFlowBoundaryCondition(0; scheme = velocity_scheme),
                                    west = NormalFlowBoundaryCondition(0; scheme = velocity_scheme))
    v_bcs = FieldBoundaryConditions(east = GradientBoundaryCondition(0), west = GradientBoundaryCondition(0))
    T_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(18.0; scheme = tracer_scheme),
                                    west = ValueBoundaryCondition(18.0; scheme = tracer_scheme))
    U_bcs = FieldBoundaryConditions(grid, (Face(), Center(), nothing),
                                    east = GravityWaveRadiationBoundaryCondition((0.0, 0.0)),
                                    west = GravityWaveRadiationBoundaryCondition((0.0, 0.0)))

    model = HydrostaticFreeSurfaceModel(grid; timestepper,
                                        free_surface = SplitExplicitFreeSurface(grid; substeps = 30),
                                        closure = CATKEVerticalDiffusivity(),
                                        coriolis = HydrostaticSphericalCoriolis(),
                                        tracers = (:T, :S), buoyancy = SeawaterBuoyancy(),
                                        momentum_advection = WENO(), tracer_advection = WENO(),
                                        boundary_conditions = (u = u_bcs, v = v_bcs, T = T_bcs, U = U_bcs))
    Random.seed!(1234)
    set!(model, T = (x, y, z) -> 20 + 0.01z + 1e-3 * randn(), S = 35,
                u = (x, y, z) -> 1e-3 * randn(), v = (x, y, z) -> 1e-3 * randn())
    return model
end

function column_model()
    grid = RectilinearGrid(size = 32, z = (-128, 0), topology = (Flat, Flat, Bounded))
    model = HydrostaticFreeSurfaceModel(grid; timestepper = :SplitRungeKutta3,
                                        closure = CATKEVerticalDiffusivity(),
                                        coriolis = FPlane(f = 1e-4), tracers = (:T, :S),
                                        buoyancy = SeawaterBuoyancy(),
                                        boundary_conditions = (u = FieldBoundaryConditions(top = FluxBoundaryCondition(-1e-4)),
                                                               T = FieldBoundaryConditions(top = FluxBoundaryCondition(1e-4))))
    set!(model, T = z -> 20 + 0.01z, S = 35)
    return model
end

function coupled_model()
    underlying_grid = LatitudeLongitudeGrid(size = (12, 10, 8), halo = (7, 7, 7),
                                            longitude = (-71, -70), latitude = (37, 38), z = (-100, 0))
    bump(λ, φ) = -100 + 60 * exp(-((λ + 70.5)^2 + (φ - 37.5)^2) / 0.15^2)
    grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bump))
    ocean = ocean_simulation(grid; Δt = 60)
    Random.seed!(1234)
    set!(ocean.model, T = (x, y, z) -> 20 + 0.01z + 1e-3 * randn(), S = 35,
                      u = (x, y, z) -> 1e-3 * randn(), v = (x, y, z) -> 1e-3 * randn())
    atmosphere = PrescribedAtmosphere(grid)
    set!(atmosphere; u = -3.0, v = 2.0, T = 288.0, q = 0.008, p = 101325.0)
    radiation = PrescribedRadiation(grid)
    set!(radiation; downwelling_shortwave = 100.0, downwelling_longwave = 300.0)
    return OceanOnlyModel(ocean; atmosphere, radiation)
end

ocean_of(model) = model isa HydrostaticFreeSurfaceModel ? model : model.ocean.model

function state(model)
    om = ocean_of(model)
    fields = merge(Oceananigans.prognostic_fields(om), (; κu = om.closure_fields.κu, κc = om.closure_fields.κc))
    return Dict(string(name) => copy(Array(parent(f))) for (name, f) in pairs(fields))
end

results = Dict{String, Any}()
for (name, build) in ("regional SplitRK3" => () -> regional_model(:SplitRungeKutta3),
                      "regional QAB2"     => () -> regional_model(:QuasiAdamsBashforth2),
                      "column SplitRK3"   => column_model,
                      "coupled OceanOnlyModel" => coupled_model)
    model = build()
    simulation = Simulation(model; Δt = 60, stop_iteration = 10, verbose = false)
    run!(simulation)
    results[name] = state(model)
end

serialize("/tmp/stage_tendency_$(TAG).jls", results)
println("wrote /tmp/stage_tendency_$(TAG).jls")
