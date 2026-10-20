using MPI
MPI.Init()

using Oceananigans
using Oceananigans.Units
using Oceananigans.DistributedComputations: partition, cpu_architecture
using NumericalEarth
using Printf

λ₀, φ₀, NX, NY, Δt = -60.0, 48.0, 80, 40, 300.0

# An ocean at 1/12° just above the freezing point, with a near-surface patch below it across the rank edge
Tᵢ(λ, φ, z) = -1.2 + 0.2 * sin(2π * (λ - λ₀) * 2) * cos(2π * (φ - φ₀) * 2) - 2z / 500
cold(λ, φ) = exp(-((λ - λ₀ - NX / 24)^2 + ((φ - φ₀ - NY / 24) / 0.6)^2))
Tˢᶜ(λ, φ, z) = Tᵢ(λ, φ, z) - 1.5 * exp(z / 60) * cold(λ, φ)
Sᵢ(λ, φ, z) = 32 - 0.5 * exp(z / 100)
uᵢ(λ, φ, z) = 0.1 * cos(2π * (φ - φ₀) * 2) * exp(z / 200)

# A cold, windy atmosphere on its own serial grid
function atmosphere_and_radiation()
    grid = LatitudeLongitudeGrid(CPU(); size = (4 * (NX ÷ 12 + 3), 4 * (NY ÷ 12 + 3)), longitude = (λ₀ - 1, λ₀ - 1 + NX ÷ 12 + 3),
                                 latitude = (φ₀ - 1, φ₀ - 1 + NY ÷ 12 + 3), topology = (Bounded, Bounded, Flat))
    atmosphere = PrescribedAtmosphere(grid, [0.0, 10days])
    radiation = PrescribedRadiation(grid, [0.0, 10days])
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

function build(arch)
    grid = LatitudeLongitudeGrid(arch; size = (NX, NY, 10), halo = (7, 7, 7), z = (-500, 0),
                                 longitude = (λ₀, λ₀ + NX / 12), latitude = (φ₀, φ₀ + NY / 12))
    ocean = ocean_simulation(grid; Δt, free_surface = SplitExplicitFreeSurface(grid; substeps = 20))
    set!(ocean.model, T = Tˢᶜ, S = Sᵢ, u = uᵢ)
    atmosphere, radiation = atmosphere_and_radiation()
    coupled = OceanOnlyModel(ocean; atmosphere, radiation)
    run!(Simulation(coupled; Δt, stop_iteration = 4, verbose = false))
    return ocean.model
end

arch = Distributed(CPU(); partition = Partition(2, 1))
distributed, serial = build(arch), build(CPU())

for name in (:u, :v, :T, :S)
    field(m) = name in (:T, :S) ? m.tracers[name] : m.velocities[name]
    p = interior(field(distributed))
    s = partition(interior(field(serial)), cpu_architecture(arch), size(p))
    Δ = MPI.Allreduce(maximum(abs, p .- s), max, MPI.COMM_WORLD)
    MPI.Comm_rank(MPI.COMM_WORLD) == 0 && @printf("%s  max |distributed - serial| = %.2e\n", name, Δ)
end
