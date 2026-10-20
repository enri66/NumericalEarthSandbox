# Runs the "Distributed bounds-preserving tracer advection" test of Oceananigans' test/mpi/hydrostatic_model.jl (the block
# below is copied from it) on 4 ranks, outside Oceananigans' test environment:
#   mpirun -np 4 julia --project=<env with Oceananigans and MPI> bounds_limiter_mpi_test.jl
using MPI
MPI.Init()
using Oceananigans, Test
using Oceananigans.Grids: MutableVerticalDiscretization
using Oceananigans.Architectures: on_architecture
using Oceananigans.DistributedComputations: partition, cpu_architecture, child_architecture, Fractional

archs = (Distributed(CPU(); synchronized_communication = false, partition = Partition(4)),
         Distributed(CPU(); synchronized_communication = false, partition = Partition(1, 4)),
         Distributed(CPU(); synchronized_communication = false, partition = Partition(2, 2)))

# A bounded tracer with a sharp front, advected across the rank edges by a uniform flow
function bounded_tracer_advection_test(grid)
    model = HydrostaticFreeSurfaceModel(grid; tracers = :c, buoyancy = nothing,
                                        free_surface = SplitExplicitFreeSurface(grid; substeps = 4),
                                        tracer_advection = WENO(order = 5, bounds = (0, 1)))

    set!(model, u = 0.1, v = 0.1, c = (x, y, z) -> ifelse(abs(x - 5e4) < 2e4 && abs(y - 5e4) < 2e4, 1, 0))

    for _ in 1:10
        time_step!(model, 600)
    end

    return model
end

@testset "Distributed bounds-preserving tracer advection" begin
    for arch in archs
        (arch.partition.x isa Fractional || arch.partition.y isa Fractional) && continue

        bounded_grid(arch) = RectilinearGrid(arch, size = (32, 32, 4), x = (0, 1e5), y = (0, 1e5), z = (-100, 0),
                                             halo = (4, 4, 4), topology = (Periodic, Periodic, Bounded))

        ms = bounded_tracer_advection_test(bounded_grid(child_architecture(arch)))
        mp = bounded_tracer_advection_test(bounded_grid(arch))

        cpu_arch = cpu_architecture(arch)
        p = interior(on_architecture(cpu_arch, mp.tracers.c))
        s = partition(interior(on_architecture(CPU(), ms.tracers.c)), cpu_arch, size(p))
        @test all(isapprox(p, s; atol = eps(eltype(p)), rtol = sqrt(eps(eltype(p)))))
    end
end
