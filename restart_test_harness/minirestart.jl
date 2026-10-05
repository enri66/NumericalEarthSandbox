# Minimal exact-restart repro: small HydrostaticFreeSurfaceModel + CATKE + SplitExplicitFreeSurface + WENO,
# no NumericalEarth/GLORYS/ERA5/tides/open-boundaries. Run to iteration 4, checkpoint every iteration to 10.
# Then restart from iteration 4 and compare iteration 5 against the uninterrupted run's iteration 5, using
# a whole-model snapshot (not just what the checkpointer saves) to see which arrays actually differ.
using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.BoundaryConditions: ObliqueRadiation, TracerReservoir, GravityWaveRadiationBoundaryCondition
using Oceananigans.OutputWriters: checkpoint
using Serialization, Random

const DIR = mktempdir()
println("scratch dir: ", DIR)
global_dir_check() = (println("checkpoint files: ", readdir(DIR)))

function build_simulation(tag)
    # MAB-like spherical grid instead of RectilinearGrid: a small lon/lat/depth box around 38N, matching
    # script 08's LatitudeLongitudeGrid + HydrostaticSphericalCoriolis (real f(φ), metric terms), still with
    # the PartialCellBottom bump from the previous step.
    underlying_grid = LatitudeLongitudeGrid(size = (12, 10, 8), halo = (4, 4, 4),
                                             longitude = (-71, -70), latitude = (37, 38), z = (-100, 0))
    bump(λ, φ) = -100 + 60 * exp(-((λ - (-70.5))^2 / 0.15^2 + (φ - 37.5)^2 / 0.15^2))
    grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bump))
    free_surface = SplitExplicitFreeSurface(grid; substeps = 30)
    closure = CATKEVerticalDiffusivity()
    coriolis = HydrostaticSphericalCoriolis()
    τx = FluxBoundaryCondition(-1e-4)

    # East/west open boundaries, the same ingredients MAB script 09a uses: ObliqueRadiation on the normal
    # velocity, GradientBoundaryCondition(0) on the tangential velocity, TracerReservoir on T/S, and
    # GravityWaveRadiation on the barotropic transport U (Flather). This is the one MAB ingredient not yet
    # tried in this cheap harness — the schemes actually radiating, not just their saved arrays.
    velocity_scheme = ObliqueRadiation(inflow_timescale = 600.0, outflow_timescale = 6000.0)
    tracer_scheme   = TracerReservoir(inflow_length_scale = 5000.0, outflow_length_scale = 2000.0)
    u_bcs = FieldBoundaryConditions(top = τx,
                                     east = NormalFlowBoundaryCondition(0; scheme = velocity_scheme),
                                     west = NormalFlowBoundaryCondition(0; scheme = velocity_scheme))
    v_bcs = FieldBoundaryConditions(east = GradientBoundaryCondition(0), west = GradientBoundaryCondition(0))
    T_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(18.0; scheme = tracer_scheme),
                                     west = ValueBoundaryCondition(18.0; scheme = tracer_scheme))
    S_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(35.0; scheme = tracer_scheme),
                                     west = ValueBoundaryCondition(35.0; scheme = tracer_scheme))
    U_bcs = FieldBoundaryConditions(grid, (Face(), Center(), nothing),
                                     east = GravityWaveRadiationBoundaryCondition((0.0, 0.0)),
                                     west = GravityWaveRadiationBoundaryCondition((0.0, 0.0)))

    model = HydrostaticFreeSurfaceModel(grid; free_surface, closure, coriolis,
                                          timestepper = :SplitRungeKutta3,  # matches NumericalEarth's ocean_simulation default
                                          tracers = (:T, :S),
                                          buoyancy = SeawaterBuoyancy(),
                                          momentum_advection = WENO(),
                                          tracer_advection = WENO(),
                                          boundary_conditions = (u = u_bcs, v = v_bcs, T = T_bcs, S = S_bcs, U = U_bcs))
    Random.seed!(1234)
    set!(model, T = (x, y, z) -> 20 + 0.01z + 1e-3 * randn(),
                S = 35, u = (x, y, z) -> 1e-3 * randn(), v = (x, y, z) -> 1e-3 * randn())
    simulation = Simulation(model; Δt = 60, stop_iteration = 4)
    simulation.output_writers[:checkpointer] = Checkpointer(model;
        schedule = IterationInterval(1), dir = DIR, prefix = tag, overwrite_files = false, cleanup = false)
    return simulation
end

# ---- state_snapshot: same walker as the MAB debug hook, trimmed ----
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
    walk((; model.velocities, model.tracers, model.free_surface, model.closure_fields,
          model.timestepper, model.pressure, model.clock,
          velocity_bcs = map(f -> f.boundary_conditions, model.velocities),
          tracer_bcs = map(f -> f.boundary_conditions, model.tracers),
          transport_velocities = try model.transport_velocities catch; nothing end), "model", 0)
    return out
end

function compare(a, b)
    rows = Any[]; same = 0
    for k in sort(collect(union(keys(a), keys(b))))
        (haskey(a, k) && haskey(b, k)) || (push!(rows, (Inf, k, "only in one")); continue)
        x, y = a[k], b[k]
        if x isa Number && y isa Number
            isequal(x, y) ? (same += 1) : push!(rows, (abs(x - y), k, "$x vs $y"))
        elseif x isa AbstractArray && y isa AbstractArray && size(x) == size(y)
            d = maximum(t -> isnan(t) ? 0.0 : abs(t), Float64.(x) .- Float64.(y); init = 0.0)
            n = count(.!isequal.(x, y))
            n == 0 ? (same += 1) : push!(rows, (d, k, "$n/$(length(x)) differ, max|Δ|=$(round(d, sigdigits=3))"))
        else
            isequal(x, y) ? (same += 1) : push!(rows, (Inf, k, "$x vs $y"))
        end
    end
    println("$same identical, $(length(rows)) differ")
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
# run!(...; pickup=...) unconditionally executes at least one step (same quirk seen in the MAB debugging),
# so this pickup call restores iteration 4 AND takes the first post-restart step, landing at iteration 5.
run!(simB; pickup = joinpath(DIR, "A_iteration4.jld2"))
snapB5 = state_snapshot(simB.model)   # first step after restart
simB.stop_iteration = 6; run!(simB)   # one more step -> iteration 6 (ordinary run!, respects stop_iteration)
snapB6 = state_snapshot(simB.model)

println("\n################ iteration 5: A (in-process) vs B (restarted from A's iteration-4 checkpoint, 1 step)")
compare(snapA5, snapB5)

println("\n################ iteration 6: A (in-process) vs B (restarted, 2 steps)")
compare(snapA6, snapB6)
