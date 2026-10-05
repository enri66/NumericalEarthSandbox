# Standalone version of test/simulation/checkpointer.jl's test_checkpoint_open_boundary_schemes:
# uninterrupted run to iteration 10 vs checkpoint at 5 + restart to 10; prints max |Δ| per prognostic field.
#   julia --project=<Oceananigans worktree> obc_checkpoint_probe.jl [QuasiAdamsBashforth2|SplitRungeKutta3]
using Oceananigans
using Oceananigans.BoundaryConditions: NormalRadiation, ObliqueRadiation, TracerReservoir, GravityWaveRadiationBoundaryCondition
using Oceananigans: prognostic_fields

const TS = Symbol(get(ARGS, 1, "QuasiAdamsBashforth2"))
const DIR = mktempdir()

function build_model(velocity_scheme, tracer_scheme; open = true)
    grid = RectilinearGrid(CPU(), size=(16, 8, 4), x=(0, 16000), y=(0, 8000), z=(-100, 0),
                           topology=(Bounded, Periodic, Bounded))
    bcs = if open
        u_bcs = FieldBoundaryConditions(east = NormalFlowBoundaryCondition(0; scheme=velocity_scheme),
                                        west = NormalFlowBoundaryCondition(0; scheme=velocity_scheme))
        b_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(0.01; scheme=tracer_scheme),
                                        west = ValueBoundaryCondition(0.01; scheme=tracer_scheme))
        U_bcs = FieldBoundaryConditions(grid, (Face(), Center(), nothing),
                                        east = GravityWaveRadiationBoundaryCondition((0.0, 0.0)),
                                        west = GravityWaveRadiationBoundaryCondition((0.0, 0.0)))
        (u=u_bcs, b=b_bcs, U=U_bcs)
    else
        NamedTuple()
    end
    model = HydrostaticFreeSurfaceModel(grid; buoyancy = BuoyancyTracer(), tracers = :b, timestepper = TS,
                                        free_surface = SplitExplicitFreeSurface(grid; substeps=10),
                                        boundary_conditions = bcs)
    set!(model, u = (x, y, z) -> 0.1 * cos(2π * x / 16000), b = (x, y, z) -> 0.005 * (1 + sin(2π * x / 16000)))
    return model
end

state(model) = Dict(string(k) => copy(Array(parent(f))) for (k, f) in pairs(prognostic_fields(model)))

function probe(name, velocity_scheme, tracer_scheme; open = true)
    A = build_model(velocity_scheme, tracer_scheme; open)
    run!(Simulation(A, Δt=10.0, stop_iteration=10, verbose=false))

    prefix = replace(name, r"[^A-Za-z]" => "_")
    B = build_model(velocity_scheme, tracer_scheme; open)
    simB = Simulation(B, Δt=10.0, stop_iteration=5, verbose=false)
    simB.output_writers[:checkpointer] = Checkpointer(B, schedule=IterationInterval(5), dir=DIR, prefix=prefix)
    run!(simB)

    C = build_model(velocity_scheme, tracer_scheme; open)
    simC = Simulation(C, Δt=10.0, stop_iteration=10, verbose=false)
    simC.output_writers[:checkpointer] = Checkpointer(C, schedule=IterationInterval(5), dir=DIR, prefix=prefix)
    set!(simC; checkpoint=:latest)
    run!(simC)

    a, c = state(A), state(C)
    diffs = [k => maximum(abs, a[k] .- c[k]) for k in sort(collect(keys(a)))]
    println(rpad("$TS  $name", 55), join(["$k=$(round(d, sigdigits=3))" for (k, d) in diffs], "  "))
end

probe("closed (no open boundaries)", nothing, nothing; open = false)
probe("NormalRadiation/NormalRadiation", NormalRadiation(inflow_timescale=100.0, outflow_timescale=1000.0),
                                         NormalRadiation(inflow_timescale=100.0, outflow_timescale=1000.0))
probe("ObliqueRadiation/TracerReservoir", ObliqueRadiation(inflow_timescale=100.0, outflow_timescale=1000.0),
                                          TracerReservoir(inflow_length_scale=5000.0, outflow_length_scale=2000.0))

# Right after restore + initialize!, before any step: do the boundary values match the run that wrote the checkpoint?
function probe_restore(name, velocity_scheme, tracer_scheme)
    prefix = replace(name, r"[^A-Za-z]" => "_") * "_restore"
    B = build_model(velocity_scheme, tracer_scheme)
    simB = Simulation(B, Δt=10.0, stop_iteration=5, verbose=false)
    simB.output_writers[:checkpointer] = Checkpointer(B, schedule=IterationInterval(5), dir=DIR, prefix=prefix)
    run!(simB)
    C = build_model(velocity_scheme, tracer_scheme)
    simC = Simulation(C, Δt=10.0, stop_iteration=10, verbose=false)
    simC.output_writers[:checkpointer] = Checkpointer(C, schedule=IterationInterval(5), dir=DIR, prefix=prefix)
    set!(simC; checkpoint=:latest)
    after_set = state(C)
    simC.initialized = false
    Oceananigans.initialize!(simC)
    after_init = state(C)
    b = state(B)
    d(s) = join(["$k=$(round(maximum(abs, b[k] .- s[k]), sigdigits=3))" for k in ("u", "b")], "  ")
    println(rpad("$TS  $name", 55), "after set!: ", d(after_set), "   after initialize!: ", d(after_init))
end
probe_restore("NormalRadiation/NormalRadiation", NormalRadiation(inflow_timescale=100.0, outflow_timescale=1000.0),
                                                 NormalRadiation(inflow_timescale=100.0, outflow_timescale=1000.0))
