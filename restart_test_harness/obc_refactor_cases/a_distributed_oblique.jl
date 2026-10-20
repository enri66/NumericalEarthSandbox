using MPI; MPI.Init()
using Test, Oceananigans
using Oceananigans.Architectures: on_architecture
using Oceananigans.DistributedComputations: child_architecture, cpu_architecture, partition, reconstruct_global_grid
XT = get(ENV, "OPEN", "") == "true" ? Open : Bounded
YT = XT
using Oceananigans.BoundaryConditions: NormalRadiation, ObliqueRadiation, GravityWaveRadiationBoundaryCondition
child_arch = CPU()
# # Distributed oblique radiation
#
# A tracer leaving through south and north boundaries with `ObliqueRadiation`, varying along them; `x` is
# periodic, so the tangential differences cross the rank edges and the periodic ends.

function build_oblique(grid)
    scheme = ObliqueRadiation(inflow_timescale = 100.0, outflow_timescale = 100.0)

    v_bcs = FieldBoundaryConditions(south = NormalFlowBoundaryCondition(0; scheme),
                                    north = NormalFlowBoundaryCondition(0; scheme))

    V_bcs = FieldBoundaryConditions(grid, (Center(), Face(), nothing);
                                    south = GravityWaveRadiationBoundaryCondition((0.0, 0.0)),
                                    north = GravityWaveRadiationBoundaryCondition((0.0, 0.0)))

    c_bcs = FieldBoundaryConditions(south = ValueBoundaryCondition(0; scheme),
                                    north = ValueBoundaryCondition(0; scheme))

    free_surface = SplitExplicitFreeSurface(grid; substeps=8, extend_halos=false)

    model = HydrostaticFreeSurfaceModel(grid; free_surface,
                                        boundary_conditions = (v = v_bcs, V = V_bcs, c = c_bcs),
                                        momentum_advection = nothing,
                                        buoyancy = nothing,
                                        tracers = :c)

    set!(model, η = (x, y, z) -> 0.01 * exp(-(y - 0.5)^2 / 0.08) * (1 + 0.2 * cos(2π * x)),
                c = (x, y, z) -> (1 + 0.5 * sin(2π * x)) * y)

    for _ in 1:60
        time_step!(model, 5e-3)
    end

    return model
end

@testset "Distributed oblique radiation" begin
    oblique_archs = (Distributed(child_arch; synchronized_communication=false, partition=Partition(4)),
                     Distributed(child_arch; synchronized_communication=false, partition=Partition(2, 2)))

    for arch in (get(ENV, "ONLY", "") == "22" ? oblique_archs[2:2] : oblique_archs)
        cpu_arch = cpu_architecture(arch)

        grid = RectilinearGrid(arch, size=(40, 20, 2), x=(0, 1), y=(0, 1), z=(-1, 0),
                               halo=(4, 4, 2), topology=(Periodic, YT, Bounded))

        mp = build_oblique(grid)                           # partitioned
        ms = build_oblique(reconstruct_global_grid(grid))  # serial reference

        cp = interior(on_architecture(cpu_arch, mp.tracers.c))
        cs = partition(interior(on_architecture(CPU(), ms.tracers.c)), cpu_arch, size(cp))

        d = MPI.Allreduce(maximum(abs, Array(cp) .- Array(cs)), max, MPI.COMM_WORLD)
        MPI.Comm_rank(MPI.COMM_WORLD) == 0 && println("PARTITION ", arch.partition.x, "x", arch.partition.y, "  max|distributed - serial| = ", d)
        @testset "oblique radiation [$(typeof(arch.partition))]" begin
            @test all(isapprox.(cp, cs))
        end
    end
end
