# Same exact-restart repro as minirestart.jl, but built through NumericalEarth's ocean_simulation and wrapped
# in OceanOnlyModel (EarthSystemModel with atmosphere=radiation=land=nothing) instead of a bare
# HydrostaticFreeSurfaceModel. Tests whether the EarthSystemModel/interfaces wrapping layer itself breaks
# exact restart, with no real atmosphere/GLORYS/ERA5 involved yet.
using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: ObliqueRadiation, TracerReservoir, GravityWaveRadiationBoundaryCondition
using NumericalEarth
using NumericalEarth.EarthSystemModels: OceanOnlyModel
using NumericalEarth.Atmospheres: PrescribedAtmosphere
using NumericalEarth.Radiations: PrescribedRadiation
using Random

const DIR = mktempdir()
println("scratch dir: ", DIR)

function build_simulation(tag)
    underlying_grid = LatitudeLongitudeGrid(size = (12, 10, 8), halo = (7, 7, 7),
                                             longitude = (-71, -70), latitude = (37, 38), z = (-100, 0))
    bump(λ, φ) = -100 + 60 * exp(-((λ - (-70.5))^2 / 0.15^2 + (φ - 37.5)^2 / 0.15^2))
    grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bump))

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

    # a trivial (non-ERA5, hand-built) prescribed atmosphere and radiation: no real data files, but the
    # temperature is a genuinely time-varying, multi-slice FieldTimeSeries (like GLORYS/ERA5's own
    # multi-time series, unlike the single-timestamp fields just tried), so time interpolation between
    # slices actually happens across the checkpoint/restart boundary.
    times = 0:60:600.0   # 11 slices spanning the whole 0..360 s test window
    T_fts = FieldTimeSeries{Center, Center, Nothing}(grid, times)
    for (n, t) in enumerate(times)
        interior(T_fts, :, :, :, n) .= 15 + t / 60   # 15, 16, 17, ... deg C — distinct per slice
    end
    atmosphere = PrescribedAtmosphere(grid; temperature = T_fts)
    set!(atmosphere; u = -3.0, v = 2.0, q = 0.008, p = 101325.0)
    radiation = PrescribedRadiation(grid)
    set!(radiation; downwelling_shortwave = 100.0, downwelling_longwave = 300.0)

    model = OceanOnlyModel(ocean; atmosphere, radiation)
    simulation = Simulation(model; Δt = 60, stop_iteration = 4)
    simulation.output_writers[:checkpointer] = Checkpointer(model;
        schedule = IterationInterval(1), dir = DIR, prefix = tag, overwrite_files = false, cleanup = false)
    return simulation
end

# ---- state_snapshot: same walker as before, rooted at the ocean sub-model plus the ESM's own top-level state ----
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

println("\n=== A: uninterrupted, iterations 4/5/6"); simA = build_simulation("A")
run!(simA)                       # to iteration 4
snapA4 = state_snapshot(simA.model)
simA.stop_iteration = 5; run!(simA)   # one more step -> iteration 5
snapA5 = state_snapshot(simA.model)
simA.stop_iteration = 6; run!(simA)   # one more step -> iteration 6
snapA6 = state_snapshot(simA.model)

println("checkpoint files: ", readdir(DIR))
println("\n=== B: restart from A's iteration-4 checkpoint")
simB = build_simulation("B")
run!(simB; pickup = joinpath(DIR, "A_iteration4.jld2"))  # restore iteration 4, unconditionally takes 1 step -> 5
snapB5 = state_snapshot(simB.model)
simB.stop_iteration = 6; run!(simB)   # one more step -> iteration 6
snapB6 = state_snapshot(simB.model)

println("\n################ iteration 5: A (in-process) vs B (restarted from A's iteration-4 checkpoint, 1 step)")
compare(snapA5, snapB5)

println("\n################ iteration 6: A (in-process) vs B (restarted, 2 steps)")
compare(snapA6, snapB6)
