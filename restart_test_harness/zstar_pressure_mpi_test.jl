# Runs the "Distributed z-star hydrostatic pressure" test of Oceananigans' test/mpi/hydrostatic_model.jl (the block
# below is copied from it) on 4 ranks, outside Oceananigans' test environment:
#   mpirun -np 4 julia --project=<env with Oceananigans and MPI> zstar_pressure_mpi_test.jl
using MPI
MPI.Init()
using Oceananigans, Test
using Oceananigans.Grids: MutableVerticalDiscretization
using Oceananigans.Architectures: on_architecture
using Oceananigans.DistributedComputations: partition, cpu_architecture, child_architecture, Fractional

archs = (Distributed(CPU(); synchronized_communication = false, partition = Partition(4)),
         Distributed(CPU(); synchronized_communication = false, partition = Partition(1, 4)),
         Distributed(CPU(); synchronized_communication = false, partition = Partition(2, 2)))

# A z-star grid with buoyancy varying in x and y, so that the hydrostatic pressure gradient includes the grid slope
function zstar_hydrostatic_pressure_test(grid)
    free_surface = SplitExplicitFreeSurface(grid; substeps = 4)
    model = HydrostaticFreeSurfaceModel(grid; free_surface, buoyancy = BuoyancyTracer(), tracers = :b)

    set!(model, b = (x, y, z) -> 1e-5 * z + 1e-3 * sin(2π * x / 1e5) * sin(2π * y / 1e5))

    for _ in 1:10
        time_step!(model, 60)
    end

    return model
end

@testset "Distributed z-star hydrostatic pressure" begin
    for arch in archs
        (arch.partition.x isa Fractional || arch.partition.y isa Fractional) && continue

        zstar_grid(arch) = RectilinearGrid(arch, size = (32, 32, 8), x = (0, 1e5), y = (0, 1e5),
                                           z = MutableVerticalDiscretization((-200, 0)),
                                           halo = (3, 3, 3), topology = (Bounded, Bounded, Bounded))

        ms = zstar_hydrostatic_pressure_test(zstar_grid(child_architecture(arch)))
        mp = zstar_hydrostatic_pressure_test(zstar_grid(arch))

        cpu_arch = cpu_architecture(arch)

        for (fp, fs) in ((mp.velocities.u, ms.velocities.u), (mp.velocities.v, ms.velocities.v),
                         (mp.tracers.b, ms.tracers.b), (mp.free_surface.displacement, ms.free_surface.displacement))
            p = interior(on_architecture(cpu_arch, fp))
            s = partition(interior(on_architecture(CPU(), fs)), cpu_arch, size(p))
            @test all(isapprox(p, s; atol = eps(eltype(fp)), rtol = sqrt(eps(eltype(fp)))))
        end
    end
end
