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
const RB = get(ENV, "RB_BACKEND", "3")
rb_backend() = RB == "all" ? InMemory() : InMemory(parse(Int, RB))

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
    T_fts = FieldTimeSeries(T_PATH, "T"; backend = rb_backend(), architecture = CPU())
    u_fts = FieldTimeSeries(U_PATH, "u"; backend = rb_backend(), architecture = CPU())
    v_fts = FieldTimeSeries(V_PATH, "v"; backend = rb_backend(), architecture = CPU())
    # Matches the REAL ERA5PrescribedAtmosphere pattern exactly: u/v/T/p get a windowed backend
    # (`era5_fts`, `time_indices_in_memory=24`), but specific humidity is built via the plain
    # `FieldTimeSeries{Center,Center,Nothing}(grid, times)` constructor with NO backend kwarg — which
    # resolves to `InMemory()` = `TotallyInMemory` (the whole record loaded, unwindowed).
    q_fts = FieldTimeSeries(Q_PATH, "q"; backend = rb_backend(), architecture = CPU())
    p_fts = FieldTimeSeries(P_PATH, "p"; backend = rb_backend(), architecture = CPU())
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

function compare_quiet(a, b)
    n = 0; worst = ("", 0.0)
    for k in intersect(keys(a), keys(b))
        x, y = a[k], b[k]
        (x isa AbstractArray{<:Number} && y isa AbstractArray{<:Number} && size(x) == size(y)) || continue
        occursin("tracers.T", k) || occursin("velocities.u", k) || continue
        d = maximum(abs, Float64.(x) .- Float64.(y); init = 0.0)
        d > 0 && (n += 1); d > worst[2] && (worst = (k, d))
    end
    println("differing T/u arrays: ", n, "   worst: ", worst)
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


# Copy one group of non-checkpointed arrays from the uninterrupted run at iteration 4 into the
# restored model, take one step, and compare with the uninterrupted run at iteration 5.
om(sim) = sim.model.ocean.model
groups = Dict(
    "none"         => m -> (),
    "timestepper"  => m -> (m.timestepper,),
    "filtered"     => m -> (m.free_surface.filtered_state,),
    "closure"      => m -> (m.closure_fields,),
    "transport"    => m -> (m.transport_velocities,),
    "pressure"     => m -> (m.pressure,),
    "all"          => m -> (m.timestepper, m.free_surface.filtered_state, m.closure_fields,
                            m.transport_velocities, m.pressure))

function refs(x, seen = IdDict())
    x isa Base.RefValue{<:Number} && return Any[x]
    (x isa Oceananigans.Grids.AbstractGrid || x isa AbstractArray || x isa Function) && return Any[]
    (x isa Union{Tuple, NamedTuple} || isstructtype(typeof(x))) || return Any[]
    if ismutable(x)
        haskey(seen, x) && return Any[]
        seen[x] = true
    end
    return reduce(vcat, [refs(getfield(x, n), seen) for n in fieldnames(typeof(x)) if isdefined(x, n)]; init = Any[])
end

function arrays(x)
    x isa Oceananigans.Fields.Field && return Any[parent(x)]
    x isa Oceananigans.Fields.AbstractField && return Any[]
    x isa AbstractArray{<:Number} && return Any[x]
    x isa Oceananigans.Grids.AbstractGrid && return Any[]
    (x isa Union{Tuple, NamedTuple} || isstructtype(typeof(x))) || return Any[]
    return reduce(vcat, [arrays(getfield(x, n)) for n in fieldnames(typeof(x)) if isdefined(x, n)]; init = Any[])
end

simA = build_simulation("A")
run!(simA)
saved = Dict(g => [copy(a) for a in arrays(f(om(simA)))] for (g, f) in groups)
saved_refs = [r[] for r in refs(simA.model)]
println("Refs in the model: ", length(saved_refs), "  values: ", saved_refs)
simA.stop_iteration = 5; run!(simA)
snapA5 = state_snapshot(simA.model)
simA.stop_iteration = 6; run!(simA)
snapA6_ref = state_snapshot(simA.model)



function full_snapshot(model)
    out = Dict{String, Any}()
    seen = IdDict{Any, Bool}()
    function walk(x, path, depth)
        depth > 16 && return
        try
            if x isa Number
                out[path] = x
            elseif x isa Base.RefValue{<:Number}
                out[path] = x[]
            elseif x isa Union{Function, Symbol, AbstractString, Type, Nothing, Missing, Oceananigans.Grids.AbstractGrid, Base.RefValue}
                return
            elseif x isa Oceananigans.Fields.Field
                out[path * ".data"] = copy(Array(parent(x)))
            elseif x isa Oceananigans.Fields.AbstractField
                return
            elseif x isa AbstractArray{<:Number}
                length(x) < 10^7 && (out[path] = copy(Array(x)))
            elseif x isa Union{Tuple, NamedTuple}
                for (k, v) in pairs(x); walk(v, path * "." * string(k), depth + 1); end
            elseif isstructtype(typeof(x))
                if ismutable(x)
                    haskey(seen, x) && return
                    seen[x] = true
                end
                for n in fieldnames(typeof(x))
                    n in (:grid, :architecture) && continue
                    isdefined(x, n) || continue
                    walk(getfield(x, n), path * "." * string(n), depth + 1)
                end
            end
        catch e
            out[path * ".ERROR"] = first(sprint(showerror, e), 80)
        end
    end
    walk(model, "esm", 0)
    return out
end
clocks(sim) = (esm = sim.model.clock, ocean = om(sim).clock, ocean_sim = sim.model.ocean.Δt,
               atmosphere = sim.model.atmosphere.clock, radiation = sim.model.radiation.clock)
showclock(c) = c isa Number ? string(c) : join(["$n=$(getfield(c, n))" for n in fieldnames(typeof(c))], " ")

simC = build_simulation("C")
simC.stop_iteration = 4; run!(simC)
println("A at iteration 4:")
foreach(((k, c),) -> println("   ", rpad(k, 11), showclock(c)), pairs(clocks(simC)))
simR = build_simulation("R")
set!(simR; checkpoint = joinpath(DIR, "A_iteration4.jld2"))
simR.initialized = false; Oceananigans.initialize!(simR)
fA4 = full_snapshot(simC.model)
fR4 = full_snapshot(simR.model)
println("\n######## full-model diff at iteration 4 (A after step 4 vs restored + initialized)")
compare(fA4, fR4)
println("restored from iteration 4 and initialized:")
foreach(((k, c),) -> println("   ", rpad(k, 11), showclock(c)), pairs(clocks(simR)))

groups["refs"] = m -> ()
for g in ("none", "all")
    for with_refs in (false, true)
        simB = build_simulation("B2_" * g * string(with_refs))
        set!(simB; checkpoint = joinpath(DIR, "A_iteration4.jld2"))
        for (dst, src) in zip(arrays(groups[g](om(simB))), get(saved, g, Any[]))
            size(dst) == size(src) && copyto!(dst, src)
        end
        with_refs && for (r, v) in zip(refs(simB.model), saved_refs); r[] = v; end
        simB.stop_iteration = 5; run!(simB)
        print(rpad("copied: $g, refs=$with_refs", 30)); compare_quiet(snapA5, state_snapshot(simB.model))
    end
end

for mark in (false, true)
    simB = build_simulation("B3_" * string(mark))
    set!(simB; checkpoint = joinpath(DIR, "A_iteration4.jld2"))
    mark && (simB.model.ocean.initialized = true)
    simB.stop_iteration = 5; run!(simB)
    print(rpad("ocean.initialized set=$mark", 30)); compare_quiet(snapA5, state_snapshot(simB.model))
    simB.stop_iteration = 6; run!(simB)
    print(rpad("   ... and one more step", 30)); compare_quiet(snapA6_ref, state_snapshot(simB.model))
end

simB = build_simulation("B4")
set!(simB; checkpoint = joinpath(DIR, "A_iteration4.jld2"))
simB.initialized = false; Oceananigans.initialize!(simB)
Oceananigans.initialize!(simB.model.ocean)
fB4 = full_snapshot(simB.model)
println("\n######## iteration 4: A after step 4 vs restored + ESM initialize! + ocean initialize!")
keep(d) = Dict(k => v for (k, v) in d if !occursin("timestepper", k) && !occursin("filtered_state", k) &&
                                        !occursin("callbacks", k) && !occursin("wall_time", k) && !occursin("initialized", k))
compare(keep(fA4), keep(fB4))

for g in ("all", "transport", "closure")
    simB = build_simulation("B5_" * g)
    set!(simB; checkpoint = joinpath(DIR, "A_iteration4.jld2"))
    for (dst, src) in zip(arrays(groups[g](om(simB))), saved[g])
        size(dst) == size(src) && copyto!(dst, src)
    end
    for (r, v) in zip(refs(simB.model), saved_refs); r[] = v; end
    simB.model.ocean.initialized = true
    simB.stop_iteration = 5; run!(simB)
    print(rpad("copied $g + refs, no ocean re-init", 40)); compare_quiet(snapA5, state_snapshot(simB.model))
end

Le_groups = Dict("Le"        => m -> (m.closure_fields.Le, m.closure_fields._tupled_implicit_linear_coefficients),
                 "transport" => m -> (m.transport_velocities,),
                 "both"      => m -> (m.closure_fields.Le, m.closure_fields._tupled_implicit_linear_coefficients,
                                      m.transport_velocities))
saved_A4 = Dict(g => [copy(a) for a in arrays(f(om(simC)))] for (g, f) in Le_groups)
for g in ("Le", "transport", "both")
    simB = build_simulation("B6_" * g)
    set!(simB; checkpoint = joinpath(DIR, "A_iteration4.jld2"))
    Oceananigans.initialize!(simB.model.ocean)       # normal restart path, done up front
    for (dst, src) in zip(arrays(Le_groups[g](om(simB))), saved_A4[g])
        size(dst) == size(src) && copyto!(dst, src)
    end
    simB.stop_iteration = 5; run!(simB)
    print(rpad("after ocean init, copied A's $g", 40)); compare_quiet(snapA5, state_snapshot(simB.model))
end

# Uninterrupted run whose ocean sub-simulation is re-initialized at the start of step 5, as a restart does
simD = build_simulation("D")
simD.stop_iteration = 4; run!(simD)
simD.model.ocean.initialized = false
simD.stop_iteration = 5; run!(simD)
snapD5 = state_snapshot(simD.model)
simD.stop_iteration = 6; run!(simD)
snapD6 = state_snapshot(simD.model)

simB = build_simulation("B7")
set!(simB; checkpoint = joinpath(DIR, "A_iteration4.jld2"))
simB.stop_iteration = 5; run!(simB)
snapB5n = state_snapshot(simB.model)
simB.stop_iteration = 6; run!(simB)
snapB6n = state_snapshot(simB.model)

println("\n######## iteration 5: uninterrupted run with ocean re-initialized at step 5 vs normal restart")
compare(snapD5, snapB5n)
println("\n######## iteration 6")
compare(snapD6, snapB6n)

using Oceananigans.TimeSteppers: time_step!
simBa = build_simulation("Ba"); set!(simBa; checkpoint = joinpath(DIR, "A_iteration4.jld2"))
simBa.initialized = false; Oceananigans.initialize!(simBa)
Oceananigans.initialize!(simBa.model.ocean)
simBb = build_simulation("Bb"); set!(simBb; checkpoint = joinpath(DIR, "A_iteration4.jld2"))
simBb.initialized = false; Oceananigans.initialize!(simBb)
time_step!(simBb.model.radiation, 60.0)
time_step!(simBb.model.atmosphere, 60.0)
Oceananigans.initialize!(simBb.model.ocean)
println("\n######## ocean initialized before vs after advancing radiation + atmosphere")
oceanonly(d) = Dict(k => v for (k, v) in d if startswith(k, "esm.ocean.model") && !occursin("callbacks", k))
compare(oceanonly(full_snapshot(simBa.model)), oceanonly(full_snapshot(simBb.model)))
