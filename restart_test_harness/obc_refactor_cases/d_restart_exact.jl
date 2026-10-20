using Oceananigans
using Oceananigans.BoundaryConditions: NormalRadiation, ObliqueRadiation, TracerReservoir, GravityWaveRadiationBoundaryCondition
XT = get(ENV, "OPEN", "") == "true" ? Open : Bounded
dir = mktempdir()

function make(timestepper, scheme, stop)
    grid = RectilinearGrid(size=(16, 8, 4), x=(0, 1e5), y=(0, 5e4), z=(-100, 0), topology=(XT, Bounded, Bounded))
    obc = NormalFlowBoundaryCondition(0; scheme)
    u_bcs = FieldBoundaryConditions(west=obc, east=obc)
    U_bcs = FieldBoundaryConditions(grid, (Face(), Center(), nothing);
                                    east = GravityWaveRadiationBoundaryCondition((0.0, 0.0)),
                                    west = GravityWaveRadiationBoundaryCondition((0.0, 0.0)))
    rsv = ValueBoundaryCondition(0; scheme=TracerReservoir(inflow_length_scale=1e4))
    b_bcs = FieldBoundaryConditions(west=rsv, east=rsv)
    model = HydrostaticFreeSurfaceModel(grid; timestepper, buoyancy=BuoyancyTracer(), tracers=:b,
                                        free_surface=SplitExplicitFreeSurface(grid; substeps=10),
                                        coriolis=FPlane(f=1e-4), boundary_conditions=(u=u_bcs, U=U_bcs, b=b_bcs))
    set!(model, b=(x, y, z) -> 1e-5 * z + 1e-4 * exp(-((x - 5e4)^2 + (y - 2.5e4)^2) / 1e8),
                u=(x, y, z) -> 0.05 * sin(2π * x / 1e5))
    return Simulation(model; Δt=60, stop_iteration=stop, verbose=false)
end

for timestepper in (:QuasiAdamsBashforth2, :SplitRungeKutta3), scheme in (NormalRadiation(), ObliqueRadiation())
    prefix = "r_$(timestepper)_$(nameof(typeof(scheme)))"
    a = make(timestepper, scheme, 20); run!(a)
    b = make(timestepper, scheme, 10)
    b.output_writers[:c] = Checkpointer(b.model; schedule=IterationInterval(10), dir, prefix)
    run!(b)
    c = make(timestepper, scheme, 20)
    c.output_writers[:c] = Checkpointer(c.model; schedule=IterationInterval(1000), dir, prefix)
    set!(c; checkpoint=:latest)
    run!(c)
    d = maximum(maximum(abs, parent(getproperty(a.model, s)[n]) .- parent(getproperty(c.model, s)[n]))
                for (s, n) in ((:velocities, :u), (:velocities, :v), (:tracers, :b)))
    println(rpad("$timestepper $(nameof(typeof(scheme)))", 45), "max |continuous - restarted| = ", d)
end
