# Distributed SplitExplicitFreeSurface with open boundaries in y and ranks split in x:
# the distributed model departs from the serial one at the open boundary next to each rank interface
# when `extend_halos = true`, and matches it when `extend_halos = false`.
#
#   mpiexec -n 2 julia --project repro.jl     (also -n 4)
#
# Same construction as the "Distributed open boundary conditions" testset of test/mpi/split_explicit_boundaries.jl,
# with the open boundaries moved from x to y: there they are on the west/east faces and the x-partition never
# extends halos across them.

using MPI
MPI.Init()

using Oceananigans
using Oceananigans.DistributedComputations: reconstruct_global_grid, partition, cpu_architecture
using Oceananigans.BoundaryConditions: NormalRadiation, GravityWaveRadiationBoundaryCondition
using Printf

ηᵢ(x, y, z) = 0.01 * exp(-(y - 0.5)^2 / 0.08) * (1 + 0.2 * cos(2π * x))

function build_open(grid; extend_halos)
    v_bcs = FieldBoundaryConditions(south = NormalFlowBoundaryCondition(0; scheme = NormalRadiation(outflow_timescale = 100.0)),
                                    north = NormalFlowBoundaryCondition(0; scheme = NormalRadiation(outflow_timescale = 100.0)))

    V_bcs = FieldBoundaryConditions(grid, (Center(), Face(), nothing);
                                    south = GravityWaveRadiationBoundaryCondition((0.0, 0.0)),
                                    north = GravityWaveRadiationBoundaryCondition((0.0, 0.0)))

    free_surface = SplitExplicitFreeSurface(grid; substeps=8, extend_halos)

    model = HydrostaticFreeSurfaceModel(grid; free_surface,
                                        boundary_conditions = (v = v_bcs, V = V_bcs),
                                        momentum_advection = nothing,
                                        buoyancy = nothing,
                                        tracers = ())
    set!(model, η = ηᵢ)

    for _ in 1:10
        time_step!(model, 5e-3)
    end

    return model
end

comm  = MPI.COMM_WORLD
rank  = MPI.Comm_rank(comm)
arch  = Distributed(CPU(); partition = Partition(x = MPI.Comm_size(comm)))

for topology in ((Bounded, Bounded, Bounded), (Periodic, Bounded, Bounded))
    grid = RectilinearGrid(arch, size=(40, 20, 2), x=(0, 1), y=(0, 1), z=(-1, 0), halo=(4, 4, 2); topology)
    global_grid = reconstruct_global_grid(grid)

    for extend_halos in (true, false)
        mp = build_open(grid; extend_halos)         # partitioned
        ms = build_open(global_grid; extend_halos)  # serial reference, on every rank

        cpu_arch = cpu_architecture(arch)
        result = map((:u, :v, :η)) do name
            fp = name == :η ? mp.free_surface.displacement : getproperty(mp.velocities, name)
            fs = name == :η ? ms.free_surface.displacement : getproperty(ms.velocities, name)
            p  = interior(fp)
            s  = partition(interior(fs), cpu_arch, size(p))
            d, idx = findmax(abs.(p .- s))
            Δ  = MPI.Allreduce(d, max, comm)
            M  = MPI.Allreduce(maximum(abs, s), max, comm)
            # whole-domain (i, j) of the largest difference: this rank's offset in x is rank × its interior size
            at = d == Δ && Δ > 0 ? (idx[1] + rank * size(grid, 1), idx[2]) : (0, 0)
            at = (MPI.Allreduce(at[1], max, comm), MPI.Allreduce(at[2], max, comm))
            (name, Δ, M, at)
        end

        if rank == 0
            @printf("%-9s x, %d ranks, extend_halos = %-5s", string(topology[1]), MPI.Comm_size(comm), extend_halos)
            for (name, Δ, M, at) in result
                @printf("   %s: max|Δ| = %.2e of max %.2e", name, Δ, M)
                Δ > 0 && @printf(" at (i, j) = %s", string(at))
            end
            println()
        end
    end
end
