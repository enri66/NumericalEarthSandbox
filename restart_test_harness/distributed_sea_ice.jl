# Does the coupled ocean–sea-ice model (NumericalEarth's ocean_simulation + sea_ice_simulation with its default EVP
# dynamics, coupled to a prescribed atmosphere) give the same answer on a PX × PY partition as serially? A cold box near
# 50°N, an ice patch that the wind and the ocean push across the rank edges, built serially and distributed, stepped
# NSTEPS coupled steps of DT seconds and compared cell by cell (ocean u, v, η, T, S; ice h, ℵ, u, v).
# EVP_SUBSTEPS sets the EVP solver's substeps (default 8). FEATURES (comma-separated): immersed (a land strip along the
# west wall and a seamount), noice (the same ocean and atmosphere without sea ice), noadvection (no ice advection,
# as in sea_ice_simulation's default), nodynamics (thermodynamics only), supercooled (a near-surface patch below the
# freezing point across the rank edges: clipped when the coupled model is built, then frazil forms there).
# FILLICE=true fills the ice thickness and concentration halos after set! (a control for the unfilled halos set!
# leaves). SHOWDIFF=true lists the differing cells. The net ocean surface fluxes (τˣ, τʸ, Jᵀ, Jˢ, Jʷ) are compared
# too; NSTEPS=0 compares the state the coupled model starts from.
#
#   PX=2 PY=1 FEATURES=immersed mpiexec -n 2 julia --project=<env> distributed_sea_ice.jl
using MPI
MPI.Init()

using Oceananigans
using Oceananigans.Units
using Oceananigans.DistributedComputations: partition, cpu_architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using NumericalEarth
using NumericalEarth.SeaIces: sea_ice_dynamics
using ClimaSeaIce.SeaIceDynamics: SplitExplicitSolver
using Printf

const FEATURES = Set(Symbol.(filter(!isempty, split(get(ENV, "FEATURES", ""), ","))))
has(f) = f in FEATURES
const NSTEPS = parse(Int, get(ENV, "NSTEPS", "4"))
const DT = parse(Float64, get(ENV, "DT", "300"))
const NX = parse(Int, get(ENV, "NX", "80"))
const NY = parse(Int, get(ENV, "NY", "40"))
const FILLICE = get(ENV, "FILLICE", "false") == "true"
const SHOWDIFF = get(ENV, "SHOWDIFF", "false") == "true"
# On a distributed grid the split-explicit EVP solver widens the halos to 2 substeps + 3 cells, which must fit in each
# rank's slab; NumericalEarth's default with an RK3 ocean is 240 substeps (483-cell halos)
const EVP_SUBSTEPS = parse(Int, get(ENV, "EVP_SUBSTEPS", "8"))
# Barotropic substeps, the same for both builds: NumericalEarth's default free surface fixes them from an estimated
# maximum Δt on a distributed grid but adapts them to Δt serially (9 are enough here at DT = 300 s)
const BAROTROPIC_SUBSTEPS = parse(Int, get(ENV, "BAROTROPIC_SUBSTEPS", "20"))
const λ₀, φ₀ = -60.0, 48.0          # the south-west corner; 1/12° cells as in the regional runs

function make_grid(arch)
    underlying = LatitudeLongitudeGrid(arch; size = (NX, NY, 10), halo = (7, 7, 7), z = (-500, 0),
                                       longitude = (λ₀, λ₀ + NX / 12), latitude = (φ₀, φ₀ + NY / 12))
    has(:immersed) || return underlying
    λc, φc = λ₀ + NX / 24, φ₀ + NY / 24
    bottom(λ, φ) = λ < λ₀ + 3 / 12 ? 10.0 : -500 + 300 * exp(-((λ - λc)^2 + (φ - φc)^2) / 0.05)
    return ImmersedBoundaryGrid(underlying, GridFittedBottom(bottom))
end

# A cold, windy atmosphere on its own serial grid, wider than the ocean box, constant in time
function make_atmosphere()
    Lλ, Lφ = NX / 12 + 2, NY / 12 + 2          # the ocean box with a 1° margin, at 1/4°
    grid = LatitudeLongitudeGrid(CPU(); size = (ceil(Int, 4Lλ), ceil(Int, 4Lφ)), longitude = (λ₀ - 1, λ₀ - 1 + ceil(Int, 4Lλ) / 4),
                                 latitude = (φ₀ - 1, φ₀ - 1 + ceil(Int, 4Lφ) / 4), topology = (Bounded, Bounded, Flat))
    times = [0.0, 10days]
    atmosphere = PrescribedAtmosphere(grid, times)
    radiation = PrescribedRadiation(grid, times)
    for n in 1:2
        set!(atmosphere.velocities.u[n], (λ, φ) -> 8 + 4 * sin(2π * (λ - λ₀)) * cos(2π * (φ - φ₀) / 1.5))
        set!(atmosphere.velocities.v[n], (λ, φ) -> -6 + 3 * cos(2π * (λ - λ₀) / 1.5))
        set!(atmosphere.temperature[n], (λ, φ) -> 253 + 3 * sin(2π * (φ - φ₀)))
        set!(atmosphere.specific_humidity[n], 5e-4)
        set!(radiation.downwelling_shortwave[n], 50)
        set!(radiation.downwelling_longwave[n], 230)
    end
    return atmosphere, radiation
end

Tᵢ(λ, φ, z) = -1.2 + 0.2 * sin(2π * (λ - λ₀) * 2) * cos(2π * (φ - φ₀) * 2) - 2z / 500
Sᵢ(λ, φ, z) = 32 - 0.5 * exp(z / 100)
cold(λ, φ) = exp(-((λ - λ₀ - NX / 24)^2 / 1.0^2 + (φ - φ₀ - NY / 24)^2 / 0.6^2))
Tˢᶜ(λ, φ, z) = Tᵢ(λ, φ, z) - 0.8 * exp(z / 60) * cold(λ, φ)
uᵢ(λ, φ, z) = 0.1 * cos(2π * (φ - φ₀) * 2) * exp(z / 200)
# a patch of ice centred in the box, so that it crosses the rank edges of a 2 × 1 and a 1 × 2 partition
patch(λ, φ) = exp(-((λ - λ₀ - NX / 24)^2 / 0.6^2 + (φ - φ₀ - NY / 24)^2 / 0.4^2))
hᵢ(λ, φ) = 0.8 * patch(λ, φ)
ℵᵢ(λ, φ) = min(1, 1.2 * patch(λ, φ))

function build(grid)
    free_surface = SplitExplicitFreeSurface(grid; substeps = BAROTROPIC_SUBSTEPS)
    ocean = ocean_simulation(grid; Δt = DT, free_surface)
    set!(ocean.model, T = has(:supercooled) ? Tˢᶜ : Tᵢ, S = Sᵢ, u = uᵢ)
    atmosphere, radiation = make_atmosphere()
    if has(:noice)
        coupled = OceanOnlyModel(ocean; atmosphere, radiation)
        sea_ice = nothing
    else
        kw = has(:noadvection) ? (;) : (; advection = WENO(order = 7))
        dynamics = has(:nodynamics) ? nothing :
                   sea_ice_dynamics(grid, ocean; solver = SplitExplicitSolver(grid; substeps = EVP_SUBSTEPS))
        sea_ice = sea_ice_simulation(grid, ocean; Δt = DT, dynamics, kw...)
        set!(sea_ice.model, h = hᵢ, ℵ = ℵᵢ)
        FILLICE && foreach(fill_halo_regions!, (sea_ice.model.ice_thickness, sea_ice.model.ice_concentration))
        coupled = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
    end
    simulation = Simulation(coupled; Δt = DT, stop_iteration = NSTEPS, verbose = false)
    run!(simulation)
    return ocean.model, isnothing(sea_ice) ? nothing : sea_ice.model, coupled.interfaces.net_fluxes.ocean
end

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
const PX = parse(Int, get(ENV, "PX", string(MPI.Comm_size(comm))))
const PY = parse(Int, get(ENV, "PY", "1"))
arch = Distributed(CPU(); partition = Partition(PX, PY))
grid = make_grid(arch)
op, ip, Jp = build(grid)
os, is, Js = build(make_grid(CPU()))
cpu_arch = cpu_architecture(arch)

ocean_fields(o) = (u = o.velocities.u, v = o.velocities.v, η = o.free_surface.displacement, T = o.tracers.T, S = o.tracers.S)
fields(o, ::Nothing) = ocean_fields(o)
fields(o, i) = merge(ocean_fields(o), (h = i.ice_thickness, ℵ = i.ice_concentration, uᵢ = i.velocities.u, vᵢ = i.velocities.v))
rank == 0 && @printf("FEATURES = %-30s NSTEPS = %d  DT = %g  EVP substeps %d  partition %d × %d\n",
                     join(sort(collect(string.(FEATURES))), ","), NSTEPS, DT, EVP_SUBSTEPS, PX, PY)
flux_fields(J) = (τˣ = J.u, τʸ = J.v, Jᵀ = J.T, Jˢ = J.S, Jʷ = J.η)
fp, fs = merge(fields(op, ip), flux_fields(Jp)), merge(fields(os, is), flux_fields(Js))
for name in keys(fp)
    p = interior(fp[name]); s = partition(interior(fs[name]), cpu_arch, size(p))
    d, idx = findmax(abs.(p .- s))
    Δ = MPI.Allreduce(d, max, comm); M = MPI.Allreduce(maximum(abs, s), max, comm)
    ri, rj, _ = arch.local_index
    at = d == Δ && Δ > 0 ? (idx[1] + (ri - 1) * size(p, 1), idx[2] + (rj - 1) * size(p, 2), idx[3]) : (0, 0, 0)
    at = (MPI.Allreduce(at[1], max, comm), MPI.Allreduce(at[2], max, comm), MPI.Allreduce(at[3], max, comm))
    rank == 0 && @printf("   %-3s max|Δ| = %.2e of max %.2e%s\n", name, Δ, M, Δ > 0 ? " at (i, j, k) = $(at)" : "")
    if SHOWDIFF
        bad = findall(abs.(p .- s) .> 0)
        if !isempty(bad)
            bi = [b[1] for b in bad] .+ (ri - 1) * size(p, 1); js = [b[2] for b in bad] .+ (rj - 1) * size(p, 2)
            ks = [b[3] for b in bad]
            println("     rank $rank $name: $(length(bad)) cells, i ∈ $(extrema(bi)), j ∈ $(extrema(js)), k ∈ $(extrema(ks)); i = $(sort(unique(bi))[1:min(end, 8)])")
        end
        MPI.Barrier(comm)
    end
end
