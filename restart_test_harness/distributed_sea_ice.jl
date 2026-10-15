# Does the coupled ocean–sea-ice model (NumericalEarth's ocean_simulation + sea_ice_simulation with its default EVP
# dynamics, coupled to a prescribed atmosphere) give the same answer on a PX × PY partition as serially? A cold box near
# 50°N, an ice patch that the wind and the ocean push across the rank edges, built serially and distributed, stepped
# NSTEPS coupled steps of DT seconds and compared cell by cell (ocean u, v, η, T, S; ice h, ℵ, u, v).
# FEATURES (comma-separated): immersed (a land strip along the west wall and a seamount), noadvection (no ice advection,
# as in sea_ice_simulation's default), nodynamics (thermodynamics only).
#
#   PX=2 PY=1 FEATURES=immersed mpiexec -n 2 julia --project=<env> distributed_sea_ice.jl
using MPI
MPI.Init()

using Oceananigans
using Oceananigans.Units
using Oceananigans.DistributedComputations: partition, cpu_architecture
using NumericalEarth
using Printf

const FEATURES = Set(Symbol.(filter(!isempty, split(get(ENV, "FEATURES", ""), ","))))
has(f) = f in FEATURES
const NSTEPS = parse(Int, get(ENV, "NSTEPS", "4"))
const DT = parse(Float64, get(ENV, "DT", "300"))
const NX = parse(Int, get(ENV, "NX", "40"))
const NY = parse(Int, get(ENV, "NY", "20"))
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
    grid = LatitudeLongitudeGrid(CPU(); size = (32, 20), longitude = (λ₀ - 1, λ₀ + 3), latitude = (φ₀ - 1, φ₀ + 1.5),
                                 topology = (Bounded, Bounded, Flat))
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
uᵢ(λ, φ, z) = 0.1 * cos(2π * (φ - φ₀) * 2) * exp(z / 200)
# a patch of ice centred in the box, so that it crosses the rank edges of a 2 × 1 and a 1 × 2 partition
patch(λ, φ) = exp(-((λ - λ₀ - NX / 24)^2 / 0.6^2 + (φ - φ₀ - NY / 24)^2 / 0.4^2))
hᵢ(λ, φ) = 0.8 * patch(λ, φ)
ℵᵢ(λ, φ) = min(1, 1.2 * patch(λ, φ))

function build(grid)
    ocean = ocean_simulation(grid; Δt = DT)
    set!(ocean.model, T = Tᵢ, S = Sᵢ, u = uᵢ)
    kw = has(:noadvection) ? (;) : (; advection = WENO(order = 7))
    kw = has(:nodynamics) ? merge(kw, (; dynamics = nothing)) : kw
    sea_ice = sea_ice_simulation(grid, ocean; Δt = DT, kw...)
    set!(sea_ice.model, h = hᵢ, ℵ = ℵᵢ)
    atmosphere, radiation = make_atmosphere()
    coupled = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
    simulation = Simulation(coupled; Δt = DT, stop_iteration = NSTEPS, verbose = false)
    run!(simulation)
    return ocean.model, sea_ice.model
end

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
const PX = parse(Int, get(ENV, "PX", string(MPI.Comm_size(comm))))
const PY = parse(Int, get(ENV, "PY", "1"))
arch = Distributed(CPU(); partition = Partition(PX, PY))
grid = make_grid(arch)
op, ip = build(grid)
os, is = build(make_grid(CPU()))
cpu_arch = cpu_architecture(arch)

fields(o, i) = (u = o.velocities.u, v = o.velocities.v, η = o.free_surface.displacement, T = o.tracers.T, S = o.tracers.S,
                h = i.ice_thickness, ℵ = i.ice_concentration, uᵢ = i.velocities.u, vᵢ = i.velocities.v)
rank == 0 && @printf("FEATURES = %-30s NSTEPS = %d  DT = %g  partition %d × %d\n",
                     join(sort(collect(string.(FEATURES))), ","), NSTEPS, DT, PX, PY)
fp, fs = fields(op, ip), fields(os, is)
for name in keys(fp)
    p = interior(fp[name]); s = partition(interior(fs[name]), cpu_arch, size(p))
    d, idx = findmax(abs.(p .- s))
    Δ = MPI.Allreduce(d, max, comm); M = MPI.Allreduce(maximum(abs, s), max, comm)
    ri, rj, _ = arch.local_index
    at = d == Δ && Δ > 0 ? (idx[1] + (ri - 1) * size(p, 1), idx[2] + (rj - 1) * size(p, 2), idx[3]) : (0, 0, 0)
    at = (MPI.Allreduce(at[1], max, comm), MPI.Allreduce(at[2], max, comm), MPI.Allreduce(at[3], max, comm))
    rank == 0 && @printf("   %-3s max|Δ| = %.2e of max %.2e%s\n", name, Δ, M, Δ > 0 ? " at (i, j, k) = $(at)" : "")
end
