# Standalone, data-free reproducer for the T,S GLORYS-restoring-sponge crash found 2026-09-27 in
# scripts/04_mab_glorys_tides_reservoirs.jl (see triton-cluster memory for the full triage): with the
# T,S sponge on, the run blows up (DomainError, sqrt of a negative real, inside TEOS10 via CATKE) after
# ~0.75 sim-days on today's rebuilt Oceananigans `everything` branch, but the SAME sponge code runs fine
# on the pre-rebuild branch. This script copies the sponge machinery (sponge_masks, sponge_parameters,
# sponge_restoring, T_sponge/S_sponge) verbatim from script 04, but replaces GLORYS with a synthetic,
# analytic "target" profile — no network, no Copernicus/CDS accounts, no cached data — so it can run
# anywhere Oceananigans + NumericalEarth's fork branch are installed, including a cloud sandbox.
#
# Run:  julia --project=<env with the fork's Oceananigans> sponge_crash_repro.jl
# Exit code 1 + "CRASHED" if it reproduces the DomainError; 0 + "survived" otherwise.

using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: NormalRadiation
using Oceananigans.Grids: inactive_cell
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
const SWP = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "SeawaterPolynomials")]  # loaded by Oceananigans
const TEOS10EquationOfState = SWP.TEOS10.TEOS10EquationOfState
using Printf

const Nx, Ny, Nz = 32, 16, 20
const Lx, Ly, Lz = 200_000.0, 100_000.0, 500.0
const Δt = 5minutes
const stop_days = parse(Float64, get(ENV, "REPRO_DAYS", "2"))

grid = RectilinearGrid(CPU(), size=(Nx, Ny, Nz), x=(0, Lx), y=(0, Ly), z=(-Lz, 0),
                       topology=(Bounded, Periodic, Bounded))

# ---------------- a smooth, GLORYS-like analytic T,S,u,v target ----------------
# Roughly seasonal-shelf-like: warmer/fresher near the surface, cooling and salting with depth,
# with a synoptic-scale (~50 km) horizontal wiggle so the sponge sees a non-trivial baroclinic signal.
T_target(x, y, z) = 18 + 6 * (z / Lz) + 0.5 * sin(2π * x / 60_000) * cos(2π * y / 40_000)
S_target(x, y, z) = 34 + 1.5 * (z / Lz) + 0.2 * sin(2π * x / 60_000 + 1) * cos(2π * y / 40_000)
u_target(x, y, z) = 0.1 * cos(2π * y / Ly) * (1 + 0.5 * z / Lz)
v_target(x, y, z) = 0.05 * sin(2π * x / Lx) * (1 + 0.5 * z / Lz)

# Two identical frames a day apart so `frame()`'s time interpolation (unmodified from script 04) has
# a well-defined window; the target is constant in time, only varying in space.
times = [0.0, 1days]
function target_fts(loc, f)
    fts = FieldTimeSeries{loc...}(grid, times)
    field = Field{loc...}(grid)
    set!(field, f)
    for n in eachindex(times)
        set!(fts[n], field)
    end
    return fts
end
fts_T = target_fts((Center, Center, Center), T_target)
fts_S = target_fts((Center, Center, Center), S_target)
fts_u = target_fts((Face, Center, Center), u_target)
fts_v = target_fts((Center, Face, Center), v_target)

# ---------------- sponge (copied verbatim from scripts/04_mab_glorys_tides_reservoirs.jl) ----------------
const SPONGE_WIDTH = 8
const SPONGE_TAU = 0.25days

@inline function frame(times, t)
    n = searchsortedlast(times, t)
    n = clamp(n, 1, length(times) - 1)
    n1, n2 = n, n + 1
    w = (t - times[n1]) / (times[n2] - times[n1])
    return n1, n2, clamp(w, 0.0, 1.0)
end

function sponge_masks(wet, W)
    Nx, Ny = size(wet)
    μx = zeros(Nx, Ny)
    for j in 1:Ny, i in 1:Nx
        dwest, deast = i - 1, Nx - i
        d = min(dwest, deast)
        μx[i, j] = d < W ? cos(π/2 * d/W)^2 : 0.0
    end
    return μx
end

@inline wet_node(i, j, k, grid, ::Center, ::Center) = !inactive_cell(i, j, k, grid)
@inline wet_node(i, j, k, grid, ::Face, ::Center) = !inactive_cell(i, j, k, grid) & !inactive_cell(i - 1, j, k, grid)
@inline wet_node(i, j, k, grid, ::Center, ::Face) = !inactive_cell(i, j, k, grid) & !inactive_cell(i, j - 1, k, grid)

const wet_all = trues(Nx, Ny)
const μ = sponge_masks(wet_all, SPONGE_WIDTH)
const Δz_col = [Lz / Nz for _ in 1:Nz]

function sponge_parameters(fts, μ, LX, LY; baroclinic)
    columns = [(i, j) for j in axes(μ, 2), i in axes(μ, 1) if μ[i, j] > 0]
    column_index = zeros(Int, size(μ))
    weights = zeros(length(columns), Nz)
    for (c, (i, j)) in enumerate(columns)
        column_index[i, j] = c
        wet = [wet_node(i, j, k, grid, LX, LY) for k in 1:Nz]
        H = sum(Δz_col[wet]; init = 0.0)
        H > 0 && (weights[c, :] .= wet .* Δz_col ./ H)
    end

    ts = collect(times)
    target = zeros(length(columns), Nz, length(ts))
    loc = Oceananigans.instantiated_location(fts)
    profile = zeros(Nz)
    for n in eachindex(ts)
        frame_field = fts[n]
        for (c, (i, j)) in enumerate(columns)
            for k in 1:Nz
                profile[k] = weights[c, k] > 0 ?
                    Oceananigans.Fields.interpolate(Oceananigans.Grids.node(i, j, k, grid, LX, LY, Center()), frame_field, loc, fts.grid) : 0.0
            end
            mean_profile = baroclinic ? sum(weights[c, k] * profile[k] for k in 1:Nz) : 0.0
            for k in 1:Nz
                target[c, k, n] = weights[c, k] > 0 ? profile[k] - mean_profile : 0.0
            end
        end
    end

    return (; column_index, μ = [μ[i, j] for (i, j) in columns], weights, target, times = ts, rate = 1 / SPONGE_TAU, baroclinic)
end

@inline function sponge_restoring(i, j, k, grid, clock, ψ, p)
    ci = clamp(i, 1, size(p.column_index, 1))
    cj = clamp(j, 1, size(p.column_index, 2))
    c = @inbounds p.column_index[ci, cj]
    (c == 0 || @inbounds(p.weights[c, k]) == 0) && return zero(eltype(grid))
    ψ̄ = zero(eltype(grid))
    if p.baroclinic
        @inbounds for k′ in axes(p.weights, 2)
            ψ̄ += p.weights[c, k′] * ψ[i, j, k′]
        end
    end
    n₁, n₂, w = frame(p.times, clock.time)
    ψᴳ = @inbounds (1 - w) * p.target[c, k, n₁] + w * p.target[c, k, n₂]
    return @inbounds p.rate * p.μ[c] * (ψᴳ - (ψ[i, j, k] - ψ̄))
end

@inline T_sponge(i, j, k, grid, clock, fields, p) = sponge_restoring(i, j, k, grid, clock, fields.T, p)
@inline S_sponge(i, j, k, grid, clock, fields, p) = sponge_restoring(i, j, k, grid, clock, fields.S, p)
@inline u_sponge(i, j, k, grid, clock, fields, p) = sponge_restoring(i, j, k, grid, clock, fields.u, p)
@inline v_sponge(i, j, k, grid, clock, fields, p) = sponge_restoring(i, j, k, grid, clock, fields.v, p)

T_forcing = Forcing(T_sponge; discrete_form = true, parameters = sponge_parameters(fts_T, μ, Center(), Center(); baroclinic = false))
S_forcing = Forcing(S_sponge; discrete_form = true, parameters = sponge_parameters(fts_S, μ, Center(), Center(); baroclinic = false))

# ---------------- open boundary (east/west), matching script 04's velocity/tracer scheme choice ----------------
scheme = NormalRadiation(inflow_timescale = 3days, outflow_timescale = 360days)

u_bcs = FieldBoundaryConditions(east = NormalFlowBoundaryCondition(0; scheme), west = NormalFlowBoundaryCondition(0; scheme))
T_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(20.0; scheme), west = ValueBoundaryCondition(20.0; scheme))
S_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(34.5; scheme), west = ValueBoundaryCondition(34.5; scheme))

model = HydrostaticFreeSurfaceModel(grid;
    buoyancy = SeawaterBuoyancy(equation_of_state = TEOS10EquationOfState()),
    tracers = (:T, :S),
    closure = CATKEVerticalDiffusivity(),
    free_surface = SplitExplicitFreeSurface(grid; substeps = 10),
    forcing = (T = T_forcing, S = S_forcing),
    boundary_conditions = (u = u_bcs, T = T_bcs, S = S_bcs))

set!(model, T = T_target, S = S_target, u = u_target, v = v_target, e = 1e-6)

simulation = Simulation(model, Δt = Δt, stop_time = stop_days * days)

function progress(sim)
    T, S = sim.model.tracers.T, sim.model.tracers.S
    @printf("%s  T∈[%.2f,%.2f]  S∈[%.2f,%.2f]  |u|=%.3f\n",
            prettytime(sim.model.clock.time), minimum(T), maximum(T), minimum(S), maximum(S),
            maximum(abs, sim.model.velocities.u))
end
add_callback!(simulation, progress, TimeInterval(3hours))

crashed = false
try
    run!(simulation)
catch err
    crashed = err isa Union{DomainError} || (err isa TaskFailedException)
    println("CAUGHT: ", sprint(showerror, err) |> x -> first(x, 300))
end

if crashed
    println("\n🔴 CRASHED at t = $(prettytime(model.clock.time)) — reproduced")
    exit(1)
else
    println("\n🟢 survived to t = $(prettytime(model.clock.time))")
    exit(0)
end
