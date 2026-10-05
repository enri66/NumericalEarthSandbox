# Same as mininumearth.jl, but the atmosphere temperature is a DISK-BACKED FieldTimeSeries opened with a
# windowed InMemory(N) backend (N < total slices) — the actual GLORYS/ERA5 mechanism (NumericalEarth's own
# DatasetBackend is the same idea, an AbstractInMemoryBackend with a `start` index that advances and reloads
# slices from a data source). This is the one thing not yet tried: the window has to advance mid-run, and its
# "which slices are currently loaded" bookkeeping is exactly the kind of state that might not survive a
# checkpoint/restart. Times are chosen so a window advance happens near the restart point (iteration 4/5).
using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: ObliqueRadiation, TracerReservoir, GravityWaveRadiationBoundaryCondition
using Oceananigans.OutputReaders: OnDisk, InMemory
using NumericalEarth
using NumericalEarth.EarthSystemModels: OceanOnlyModel
using NumericalEarth.Atmospheres: PrescribedAtmosphere
using NumericalEarth.Radiations: PrescribedRadiation
using Random

const DIR = mktempdir()
println("scratch dir: ", DIR)

const UNDERLYING_GRID = LatitudeLongitudeGrid(size = (12, 10, 8), halo = (7, 7, 7),
                                               longitude = (-71, -70), latitude = (37, 38), z = (-100, 0))
const BUMP(λ, φ) = -100 + 60 * exp(-((λ - (-70.5))^2 / 0.15^2 + (φ - 37.5)^2 / 0.15^2))
const GRID = ImmersedBoundaryGrid(UNDERLYING_GRID, PartialCellBottom(BUMP))

# Write a small multi-slice temperature series to disk once (shared by both A and B, like a shared GLORYS
# cache file), 11 slices at 60 s spacing spanning the whole 0..600 s test window, each with a distinct value.
const T_TIMES = 0:60:600.0
const T_PATH = joinpath(DIR, "T_series.jld2")
const U_PATH = joinpath(DIR, "u_series.jld2")
const V_PATH = joinpath(DIR, "v_series.jld2")
const Q_PATH = joinpath(DIR, "q_series.jld2")
const P_PATH = joinpath(DIR, "p_series.jld2")
function write_series(path, name, valuefn)
    f_tmp = Field{Center, Center, Nothing}(GRID)
    f = FieldTimeSeries{Center, Center, Nothing}(GRID, T_TIMES; backend = OnDisk(), path, name)
    for (n, t) in enumerate(T_TIMES)
        set!(f_tmp, valuefn(t))
        set!(f, f_tmp, n)
    end
end
write_series(T_PATH, "T", t -> 15 + t / 60)
write_series(U_PATH, "u", t -> -3.0 + t / 600)
write_series(V_PATH, "v", t -> 2.0 - t / 600)
write_series(Q_PATH, "q", t -> 0.008 + t / 6.0e5)
write_series(P_PATH, "p", t -> 101325.0 + t)

function build_simulation(tag)
    grid = GRID
    velocity_scheme = ObliqueRadiation(inflow_timescale = 600.0, outflow_timescale = 6000.0)
    tracer_scheme   = TracerReservoir(inflow_length_scale = 5000.0, outflow_length_scale = 2000.0)
    u_bcs = FieldBoundaryConditions(east = NormalFlowBoundaryCondition(0; scheme = velocity_scheme),
                                     west = NormalFlowBoundaryCondition(0; scheme = velocity_scheme))
    v_bcs = FieldBoundaryConditions(east = GradientBoundaryCondition(0), west = GradientBoundaryCondition(0))
    T_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(18.0; scheme = tracer_scheme),
                                     west = ValueBoundaryCondition(18.0; scheme = tracer_scheme))
    S_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(35.0; scheme = tracer_scheme),
                                     west = ValueBoundaryCondition(35.0; scheme = tracer_scheme))
    U_bcs = FieldBoundaryConditions(grid, (Face(), Center(), nothing),
                                     east = GravityWaveRadiationBoundaryCondition((0.0, 0.0)),
                                     west = GravityWaveRadiationBoundaryCondition((0.0, 0.0)))

    ocean = ocean_simulation(grid; Δt = 60,
                              boundary_conditions = (u = u_bcs, v = v_bcs, T = T_bcs, S = S_bcs, U = U_bcs))
    Random.seed!(1234)
    set!(ocean.model, T = (x, y, z) -> 20 + 0.01z + 1e-3 * randn(),
                       S = 35, u = (x, y, z) -> 1e-3 * randn(), v = (x, y, z) -> 1e-3 * randn())

    # Reopen the shared file with a windowed backend: only 3 of the 11 slices in memory at once, so the
    # window has to advance (reload from T_PATH) as the clock passes each slice's time. ALL atmosphere
    # fields (u, v, T, q, p) get the SAME multi-slice, same-window treatment — matching how a real ERA5
    # atmosphere is built (every field windowed in lockstep), to rule out a mismatched-backend artifact.
    T_fts = FieldTimeSeries(T_PATH, "T"; backend = InMemory(3), architecture = CPU())
    u_fts = FieldTimeSeries(U_PATH, "u"; backend = InMemory(3), architecture = CPU())
    v_fts = FieldTimeSeries(V_PATH, "v"; backend = InMemory(3), architecture = CPU())
    # Matches the REAL ERA5PrescribedAtmosphere pattern exactly: u/v/T/p get a windowed backend
    # (`era5_fts`, `time_indices_in_memory=24`), but specific humidity is built via the plain
    # `FieldTimeSeries{Center,Center,Nothing}(grid, times)` constructor with NO backend kwarg — which
    # resolves to `InMemory()` = `TotallyInMemory` (the whole record loaded, unwindowed).
    q_fts = FieldTimeSeries(Q_PATH, "q"; backend = InMemory(), architecture = CPU())
    p_fts = FieldTimeSeries(P_PATH, "p"; backend = InMemory(3), architecture = CPU())
    atmosphere = PrescribedAtmosphere(grid; velocities = (u = u_fts, v = v_fts),
                                       temperature = T_fts, specific_humidity = q_fts, pressure = p_fts)
    radiation = PrescribedRadiation(grid)
    set!(radiation; downwelling_shortwave = 100.0, downwelling_longwave = 300.0)

    model = OceanOnlyModel(ocean; atmosphere, radiation)
    simulation = Simulation(model; Δt = 60, stop_iteration = 4)
    simulation.output_writers[:checkpointer] = Checkpointer(model;
        schedule = IterationInterval(1), dir = DIR, prefix = tag, overwrite_files = false, cleanup = false)
    return simulation
end

# ---- state_snapshot / compare: same as mininumearth.jl ----
function state_snapshot(model)
    out = Dict{String, Any}()
    seen = IdDict{Any, Bool}()
    skipped = Union{Function, Symbol, AbstractString, Type, Nothing, Missing, Oceananigans.Grids.AbstractGrid}
    function walk(x, path, depth)
        depth > 14 && return
        try
            if x isa Number
                out[path] = x
            elseif x isa skipped
                return
            elseif x isa Oceananigans.Fields.Field
                walk(parent(x), path * ".data", depth + 1)
            elseif x isa AbstractArray{<:Number}
                length(x) < 10^7 && (out[path] = copy(Array(x)))
            elseif x isa Union{Tuple, NamedTuple}
                for (k, v) in pairs(x); walk(v, path * "." * string(k), depth + 1); end
            elseif isstructtype(typeof(x))
                if ismutable(x)
                    haskey(seen, x) && return
                    seen[x] = true
                end
                for name in fieldnames(typeof(x))
                    name in (:grid, :architecture) && continue
                    isdefined(x, name) || continue
                    walk(getfield(x, name), path * "." * string(name), depth + 1)
                end
            end
        catch e
            out[path * ".ERROR"] = first(sprint(showerror, e), 100)
        end
    end
    om = model.ocean.model
    walk((; om.velocities, om.tracers, om.free_surface, om.closure_fields,
          om.timestepper, om.pressure, om.clock,
          velocity_bcs = map(f -> f.boundary_conditions, om.velocities),
          tracer_bcs = map(f -> f.boundary_conditions, om.tracers),
          interfaces = model.interfaces, esm_clock = model.clock,
          atmosphere = model.atmosphere, radiation = model.radiation,
          transport_velocities = try om.transport_velocities catch; nothing end), "model", 0)
    return out
end

function compare(a, b)
    rows = Any[]; same = Ref(0)
    for k in sort(collect(union(keys(a), keys(b))))
        (haskey(a, k) && haskey(b, k)) || (push!(rows, (Inf, k, "only in one")); continue)
        x, y = a[k], b[k]
        if x isa Number && y isa Number
            isequal(x, y) ? (same[] += 1) : push!(rows, (abs(x - y), k, "$x vs $y"))
        elseif !(x isa AbstractArray{<:Number}) || !(y isa AbstractArray{<:Number})
            isequal(x, y) ? (same[] += 1) : push!(rows, (Inf, k, "$(first(string(x),60)) vs $(first(string(y),60))"))
        elseif size(x) == size(y)
            d = maximum(t -> isnan(t) ? 0.0 : abs(t), Float64.(x) .- Float64.(y); init = 0.0)
            n = count(.!isequal.(x, y))
            n == 0 ? (same[] += 1) : push!(rows, (d, k, "$n/$(length(x)) differ, max|Δ|=$(round(d, sigdigits=3))"))
        else
            push!(rows, (Inf, k, "size $(size(x)) vs $(size(y))"))
        end
    end
    println("$(same[]) identical, $(length(rows)) differ")
    for (d, k, msg) in sort(rows; by = r -> -r[1])
        println(rpad(k, 55), "  ", msg)
    end
end

using Oceananigans.OutputReaders: cpu_interpolating_time_indices
function probe_atmos(sim, label)
    m = sim.model
    a = m.atmosphere
    T_fts = a.temperature
    wanted = cpu_interpolating_time_indices(CPU(), T_fts.times, T_fts.time_indexing, m.clock.time)
    println(label, ": clock.time=", m.clock.time, "  atmos.clock.time=", a.clock.time,
            "  WANTED n1,n2=", (wanted.first_index, wanted.second_index),
            "  LOADED backend=", a.temperature.backend, " (start, length)",
            "  window[1]=", interior(a.temperature)[1,1,1,1], "  window[2]=", interior(a.temperature)[1,1,1,2],
            "  exchanger.T[1,1]=", interior(m.interfaces.exchanger.atmosphere.state.T)[1,1],
            "  q.backend=", a.specific_humidity.backend, "  exchanger.q[1,1]=", interior(m.interfaces.exchanger.atmosphere.state.q)[1,1])
end

println("\n=== A: step by step from iteration 0 ==="); simA = build_simulation("A_stepwise")
for it in 1:6
    simA.stop_iteration = it
    run!(simA)
    probe_atmos(simA, "A@iter$it")
end

println("\n=== A: uninterrupted, iterations 4/5/6"); simA = build_simulation("A")
run!(simA)                       # to iteration 4
snapA4 = state_snapshot(simA.model)
simA.stop_iteration = 5; run!(simA)   # one more step -> iteration 5 (crosses the 240s/300s slice boundary)
probe_atmos(simA, "A@iter5")
snapA5 = state_snapshot(simA.model)
simA.stop_iteration = 6; run!(simA)   # one more step -> iteration 6
snapA6 = state_snapshot(simA.model)

println("checkpoint files: ", readdir(DIR))
println("\n=== B: restart from A's iteration-4 checkpoint")
simB = build_simulation("B")
probe_atmos(simB, "B@iter4-pre-pickup(build)")
run!(simB; pickup = joinpath(DIR, "A_iteration4.jld2"))  # restore iteration 4, takes 1 step -> 5
probe_atmos(simB, "B@iter5-post-restart-step")
snapB5 = state_snapshot(simB.model)
simB.stop_iteration = 6; run!(simB)   # one more step -> iteration 6
snapB6 = state_snapshot(simB.model)

println("\n################ iteration 5: A (in-process) vs B (restarted from A's iteration-4 checkpoint, 1 step)")
compare(snapA5, snapB5)

println("\n################ iteration 6: A (in-process) vs B (restarted, 2 steps)")
compare(snapA6, snapB6)
