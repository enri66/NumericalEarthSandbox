using Oceananigans
XT = get(ENV, "OPEN", "") == "true" ? Open : Bounded
YT = XT
using Oceananigans.BoundaryConditions: TracerReservoir
# Uniform flow u = 1 through a channel; the inflow value at the west boundary oscillates in time, cᵂ(t) = sin(2π t / 0.5),
# so the interior near the boundary must follow the boundary value, not an extrapolation of the interior.
grid = RectilinearGrid(size = (64, 4, 4), x = (0, 1), y = (0, 1), z = (0, 1), halo = (4, 4, 4), topology = (XT, Periodic, Bounded))
cᵂ(y, z, t) = sin(2π * t / 0.5)
for (name, west) in (("Value", ValueBoundaryCondition(cᵂ)),
                     ("Value+TracerReservoir", ValueBoundaryCondition(cᵂ; scheme = TracerReservoir(inflow_length_scale = 0))))
    for order in (5, 7)
        c_bcs = FieldBoundaryConditions(west = west, east = ValueBoundaryCondition(0; scheme = TracerReservoir()))
        model = HydrostaticFreeSurfaceModel(grid; velocities = PrescribedVelocityFields(u = 1),
                                            tracer_advection = WENO(; order), tracers = :c, buoyancy = nothing,
                                            boundary_conditions = (; c = c_bcs))
        for _ in 1:730
            time_step!(model, 0.001)
        end
        c = Array(interior(model.tracers.c, :, 2, 2))
        exact = [sin(2π * (model.clock.time - (i - 0.5) / 64) / 0.5) for i in 1:16]
        println(rpad("$name WENO($order)", 30), "c[1:6] = ", round.(c[1:6]; digits = 3), "  exact ", round.(exact[1:6]; digits = 3),
                "  max|c - exact| over 16 cells ", round(maximum(abs, c[1:16] .- exact); digits = 3))
    end
end
