# Which model feature makes a distributed HydrostaticFreeSurfaceModel differ from the serial one? A MAB-like box built
# serially and on an x-partition of all ranks, stepped NSTEPS times, compared cell by cell. Features are switched on
# with FEATURES (comma-separated): latlon, immersed, teos10, catke, weno, coriolis, rk3, open, avid (WENO advection with
# AdaptiveVerticallyImplicitDiscretization, as in NumericalEarth's ocean_simulation), drag (semi-implicit quadratic bottom
# drag on the bottom and immersed boundaries, whose implicit coefficient reads the other velocity component), mabgrid
# (latitude-longitude bounds on multiples of 1/12 degree, which are not exact binary fractions, as in the MAB runs),
# cooling (a surface heat loss varying in x and y, so that CATKE's averaged surface buoyancy flux is not zero), oblique
# (with open: ObliqueRadiation on the normal v and the tangential u at the south and north boundaries). Without features the
# model is the plain closed box of split_explicit_open_y_distributed.jl.
#
#   FEATURES=latlon,immersed,teos10 julia --project=<env> launch_mpi.jl 2 distributed_features.jl
using MPI
MPI.Init()

using Oceananigans
using Oceananigans.Units
using Oceananigans.DistributedComputations: partition, cpu_architecture
using Oceananigans.BoundaryConditions: NormalRadiation, ObliqueRadiation, GravityWaveRadiationBoundaryCondition
using Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using Printf
using Profile    # lets SIGINFO / SIGUSR1 print the stacks of a hung rank

const FEATURES = Set(Symbol.(filter(!isempty, split(get(ENV, "FEATURES", ""), ","))))
has(f) = f in FEATURES
const NSTEPS = parse(Int, get(ENV, "NSTEPS", "2"))

function make_grid(arch)
    sz = (parse(Int, get(ENV, "NX", "40")), parse(Int, get(ENV, "NY", "20")), 8)
    z = (-1000, 0)
    underlying = has(:latlon) ?
        LatitudeLongitudeGrid(arch; size = sz, z, halo = (7, 7, 7),
                              longitude = has(:mabgrid) ? (-76, -76 + 40 / 12) : (-74, -70),
                              latitude  = has(:mabgrid) ? (34, 34 + 20 / 12)   : (36, 38)) :
        RectilinearGrid(arch; size = sz, x = (0, 400kilometers), y = (0, 200kilometers), z, halo = (7, 7, 7),
                        topology = (Bounded, Bounded, Bounded))
    has(:immersed) || return underlying
    x₀, y₀ = has(:mabgrid) ? (-74.3, 34.8) : has(:latlon) ? (-72, 37) : (200e3, 100e3)
    bump(x, y) = -1000 + 600 * exp(-((x - x₀)^2 + (y - y₀)^2) /
                                     (has(:latlon) ? 0.5 : 5e9))
    return ImmersedBoundaryGrid(underlying, GridFittedBottom(bump))
end

# with cooling, a surface mixed layer 300 m deep, which the cooling makes convective
Tᵢ(x, y, z) = 15 + 10 * exp((has(:cooling) ? min(z + 300, 0) : z) / 300) +
              0.5 * sin(2π * x / (has(:latlon) ? 4 : 400e3)) * cos(2π * y / (has(:latlon) ? 2 : 200e3))
Sᵢ(x, y, z) = 35 - 0.5 * exp(z / 200)
uᵢ(x, y, z) = 0.1 * cos(2π * y / (has(:latlon) ? 2 : 200e3)) * exp(z / 500)

@inline uspeed(i, j, k, grid, Φ) = @inbounds sqrt(Φ.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, Φ.v)^2 + 1e-4)
@inline vspeed(i, j, k, grid, Φ) = @inbounds sqrt(Φ.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, Φ.u)^2 + 1e-4)
@inline u_bottom_λ(i, j, grid, clock, Φ, μ) = - μ * uspeed(i, j, 1, grid, Φ)
@inline v_bottom_λ(i, j, grid, clock, Φ, μ) = - μ * vspeed(i, j, 1, grid, Φ)
@inline u_immersed_λ(i, j, k, grid, clock, Φ, μ) = - μ * uspeed(i, j, k, grid, Φ)
@inline v_immersed_λ(i, j, k, grid, clock, Φ, μ) = - μ * vspeed(i, j, k, grid, Φ)

function drag_bcs(λ_bottom, λ_immersed)
    bottom = IMEXFluxBoundaryCondition(0, λ_bottom; discrete_form = true, parameters = 3e-3)
    immersed = ImmersedBoundaryCondition(bottom = IMEXFluxBoundaryCondition(0, λ_immersed; discrete_form = true, parameters = 3e-3))
    return has(:immersed) ? (; bottom, immersed) : (; bottom)
end

dbg(msg) = haskey(ENV, "DEBUG_TAGS") && (println("rank ", MPI.Comm_rank(MPI.COMM_WORLD), ": ", msg); flush(stdout))

function build(grid)
    tracers = has(:teos10) ? (:T, :S) : ()
    buoyancy = has(:teos10) ? SeawaterBuoyancy(equation_of_state = TEOS10EquationOfState()) : nothing
    closure = has(:catke) ? CATKEVerticalDiffusivity() : nothing
    avid = AdaptiveVerticallyImplicitDiscretization(cfl = 0.5)
    momentum_advection = has(:avid) ? WENOVectorInvariant(time_discretization = avid) : has(:weno) ? WENOVectorInvariant() : nothing
    tracer_advection = has(:avid) ? WENO(order = 7, time_discretization = avid) : has(:weno) ? WENO(order = 7) : Centered()
    coriolis = has(:coriolis) ? HydrostaticSphericalCoriolis() : nothing
    timestepper = has(:rk3) ? :SplitRungeKutta3 : :QuasiAdamsBashforth2
    # Open south and north boundaries: the normal velocity v radiated, or with oblique radiation also the tangential u
    oblique = ObliqueRadiation(inflow_timescale = 1days, outflow_timescale = 360days)
    v_scheme = has(:oblique) ? oblique : NormalRadiation(outflow_timescale = 1days)
    v_sides = has(:open) ? (south = NormalFlowBoundaryCondition(0; scheme = v_scheme),
                            north = NormalFlowBoundaryCondition(0; scheme = v_scheme)) : (;)
    u_sides = has(:open) && has(:oblique) ? (south = ValueBoundaryCondition(0; scheme = oblique),
                                             north = ValueBoundaryCondition(0; scheme = oblique)) : (;)
    # with east: the east boundary open too, u normal and (with oblique) v tangential, as in the MAB
    if has(:open) && has(:east)
        u_sides = merge(u_sides, (east = NormalFlowBoundaryCondition(0; scheme = v_scheme),))
        has(:oblique) && (v_sides = merge(v_sides, (east = ValueBoundaryCondition(0; scheme = oblique),)))
    end
    u_drag = has(:drag) ? drag_bcs(u_bottom_λ, u_immersed_λ) : (;)
    v_drag = has(:drag) ? drag_bcs(v_bottom_λ, v_immersed_λ) : (;)
    boundary_conditions = (u = FieldBoundaryConditions(; u_sides..., u_drag...),
                           v = FieldBoundaryConditions(; v_sides..., v_drag...))
    if has(:open)
        V = FieldBoundaryConditions(grid, (Center(), Face(), nothing);
                                    south = GravityWaveRadiationBoundaryCondition((0.0, 0.0)),
                                    north = GravityWaveRadiationBoundaryCondition((0.0, 0.0)))
        boundary_conditions = merge(boundary_conditions, (; V))
        if has(:east)
            U = FieldBoundaryConditions(grid, (Face(), Center(), nothing); east = GravityWaveRadiationBoundaryCondition((0.0, 0.0)))
            boundary_conditions = merge(boundary_conditions, (; U))
        end
    end
    if has(:cooling)
        Jᵀ(x, y, t) = 1e-3 * (1 + 0.5 * sin(2π * x / (has(:latlon) ? 1 : 100e3)) * cos(2π * y / (has(:latlon) ? 1 : 100e3)))
        boundary_conditions = merge(boundary_conditions, (T = FieldBoundaryConditions(top = FluxBoundaryCondition(Jᵀ)),))
    end
    dbg("boundary conditions built")
    free_surface = SplitExplicitFreeSurface(grid; substeps = 10)
    dbg("free surface built, pid $(getpid())")
    model = HydrostaticFreeSurfaceModel(grid; free_surface, tracers, buoyancy, closure, momentum_advection, tracer_advection,
                                        coriolis, timestepper, boundary_conditions)
    dbg("model built")
    has(:teos10) ? set!(model, T = Tᵢ, S = Sᵢ, u = uᵢ) : set!(model, u = uᵢ)
    dbg("model set")
    a = Oceananigans.architecture(grid)
    if a isa Distributed && haskey(ENV, "DEBUG_TAGS")
        println("rank ", MPI.Comm_rank(MPI.COMM_WORLD), " index ", a.local_index, " distributed fields ", a.field_count[])
        flush(stdout)
    end
    Oceananigans.initialize!(model)    # as a Simulation does before its first step
    for _ in 1:NSTEPS
        time_step!(model, 60)
    end
    return model
end

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
const PX = parse(Int, get(ENV, "PX", string(MPI.Comm_size(comm))))
const PY = parse(Int, get(ENV, "PY", "1"))
arch = Distributed(CPU(); partition = Partition(PX, PY))
grid = make_grid(arch)
mp = build(grid)
ms = build(make_grid(CPU()))
cpu_arch = cpu_architecture(arch)
names = has(:teos10) ? (:u, :v, :η, :T, :S) : (:u, :v, :η)
rank == 0 && @printf("FEATURES = %-45s NSTEPS = %d  ranks = %d\n", join(sort(collect(string.(FEATURES))), ","), NSTEPS, MPI.Comm_size(comm))
for name in names
    pick(m) = name == :η ? m.free_surface.displacement : name in (:T, :S) ? getproperty(m.tracers, name) : getproperty(m.velocities, name)
    p = interior(pick(mp)); s = partition(interior(pick(ms)), cpu_arch, size(p))
    d, idx = findmax(abs.(p .- s))
    Δ = MPI.Allreduce(d, max, comm); M = MPI.Allreduce(maximum(abs, s), max, comm)
    ri, rj, _ = arch.local_index
    at = d == Δ && Δ > 0 ? (idx[1] + (ri - 1) * size(grid, 1), idx[2] + (rj - 1) * size(grid, 2), idx[3]) : (0, 0, 0)
    at = (MPI.Allreduce(at[1], max, comm), MPI.Allreduce(at[2], max, comm), MPI.Allreduce(at[3], max, comm))
    rank == 0 && @printf("   %s: max|Δ| = %.2e of max %.2e%s\n", name, Δ, M, Δ > 0 ? " at (i, j, k) = $(at)" : "")
    if haskey(ENV, "SHOWDIFF")
        bad = findall(abs.(p .- s) .> 1e-12)
        if !isempty(bad)
            is = [b[1] for b in bad] .+ (ri - 1) * size(grid, 1); js = [b[2] for b in bad] .+ (rj - 1) * size(grid, 2)
            ks = [b[3] for b in bad]
            println("     rank $rank $name: $(length(bad)) cells, i ∈ $(extrema(is)), j ∈ $(extrema(js)), k ∈ $(extrema(ks))")
        end
        MPI.Barrier(comm)
    end
end
