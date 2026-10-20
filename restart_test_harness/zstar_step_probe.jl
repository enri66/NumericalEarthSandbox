# Step a z-star case NSTEPS times on the serial and the distributed model (distributed_features.jl), then take the
# next Adams-Bashforth step one operation at a time on both, comparing every internal field after each (with NSTEPS=0
# the first step: the model's first update and an Euler step).
#   PROBE=1 NSTEPS=1 SHOWDIFF=1 FEATURES=latlon,teos10,zstar mpirun -np 2 julia --project=<env> zstar_step_probe.jl
include("distributed_features.jl")

using Oceananigans.Models.HydrostaticFreeSurfaceModels: update_transport_velocities!, compute_momentum_flux_bcs!,
    ab2_step_velocities!, mask_immersed_horizontal_velocities!, compute_free_surface_tendency!, step_free_surface!,
    compute_transport_velocities!, fill_velocity_halos_read_by_boundary_conditions!, compute_tracer_tendencies!,
    ab2_step_grid!, correct_barotropic_mode!, ab2_step_tracers!, cache_previous_tendencies!
using Oceananigans.Models.HydrostaticFreeSurfaceModels: compute_hydrostatic_momentum_tendencies!,
    complete_communication_and_compute_momentum_buffer!, compute_coriolis_prognostic_tendencies!
using Oceananigans.Models: update_model_field_time_series!, surface_kernel_parameters, volume_kernel_parameters,
    interior_tendency_kernel_parameters
using Oceananigans.Models.NonhydrostaticModels: update_hydrostatic_pressure!
using Oceananigans.BuoyancyFormulations: compute_buoyancy_gradients!
using Oceananigans.TurbulenceClosures: compute_closure_fields!
using Oceananigans.Forcings: compute_forcing!
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using Oceananigans.Utils: get_active_cells_map
using Oceananigans.Operators: Δx⁻¹ᶠᶜᶜ, ∂xᶠᶜᶜ
using Oceananigans.TimeSteppers: tick!, update_state!, maybe_prepare_first_time_step!
using Oceananigans.BoundaryConditions: fill_halo_regions!, update_boundary_conditions!
using Oceananigans.Models.HydrostaticFreeSurfaceModels: update_vertical_velocities!

const Δt = 60.0
const first_step = NSTEPS == 0
χ(m) = convert(eltype(m.grid), first_step ? -0.5 : m.timestepper.χ)

function diagnostics!(m)
    g = m.grid
    compute_buoyancy_gradients!(m.buoyancy, g, m.tracers, parameters = volume_kernel_parameters(g))
    update_vertical_velocities!(m.velocities, g, m, parameters = surface_kernel_parameters(g))
    update_hydrostatic_pressure!(m.pressure.pHY′, Oceananigans.architecture(g), g, m.buoyancy, m.tracers, parameters = surface_kernel_parameters(g))
    compute_closure_fields!(m.closure_fields, m.closure, m, parameters = surface_kernel_parameters(g))
    fill_halo_regions!(m.closure_fields; only_local_halos = true)
    fill_halo_regions!(m.pressure.pHY′; only_local_halos = true)
end

interior_momentum_tendencies!(m) =
    compute_hydrostatic_momentum_tendencies!(m, m.velocities, interior_tendency_kernel_parameters(Oceananigans.architecture(m.grid), m.grid);
                                             active_cells_map = get_active_cells_map(m.grid, Val(:core)))

operations = [
    ("first update", m -> first_step && maybe_prepare_first_time_step!(m, Δt, [])),
    ("velocities", m -> (update_transport_velocities!(m.transport_velocities, m.velocities, m.free_surface);
                         compute_momentum_flux_bcs!(m); ab2_step_velocities!(m.velocities, m, Δt, χ(m));
                         mask_immersed_horizontal_velocities!(m.velocities))),
    ("free surface tendency", m -> compute_free_surface_tendency!(m.grid, m, m.free_surface, Δt)),
    ("free surface step", m -> step_free_surface!(m.free_surface, m, m.timestepper, Δt)),
    ("transport velocities", m -> (compute_transport_velocities!(m, m.free_surface);
                                   fill_velocity_halos_read_by_boundary_conditions!(m.velocities);
                                   fill_halo_regions!((m.velocities.u, m.velocities.v), m.clock, Oceananigans.fields(m); async = true))),
    ("tracer tendencies", m -> compute_tracer_tendencies!(m)),
    ("grid", m -> ab2_step_grid!(m.grid, m, m.vertical_coordinate, Δt, χ(m))),
    ("barotropic correction, tracers", m -> (correct_barotropic_mode!(m, Δt); ab2_step_tracers!(m.tracers, m, Δt, χ(m)))),
    ("cache, tick", m -> (cache_previous_tendencies!(m); tick!(m.clock, Δt); m.clock.last_Δt = Δt)),
    ("update state: boundary conditions", m -> (mask_immersed_field!(m.tracers); update_model_field_time_series!(m, m.clock);
                                                 compute_forcing!(m.forcing); update_boundary_conditions!(Oceananigans.fields(m), m))),
    ("update state: tracer halos (async)", m -> fill_halo_regions!(m.tracers, m.clock, Oceananigans.fields(m); async = true)),
    ("update state: diagnostics", diagnostics!),
    ("update state: interior momentum tendencies", interior_momentum_tendencies!),
    ("update state: buffer momentum tendencies", m -> complete_communication_and_compute_momentum_buffer!(m, m.grid, Oceananigans.architecture(m.grid)))]

for (label, op!) in operations
    op!(mp)
    op!(ms)
    rank == 0 && println("== after ", label)
    MPI.Barrier(comm)
    deep_compare(mp, ms, grid)
    MPI.Barrier(comm)
end

# PROBE_AT="i,j,k" (local index on the last rank): the values the u tendency there is computed from, on both models
if haskey(ENV, "PROBE_AT") && rank == MPI.Comm_size(comm) - 1
    i, j, k = parse.(Int, split(ENV["PROBE_AT"], r"[,:]"))
    gi = i + (PX - 1) * size(grid, 1)

    for (label, m, ii) in (("distributed", mp, i), ("serial", ms, gi))
        p = m.pressure.pHY′
        @printf("  %-11s pHY′[i-1] %.17g  pHY′[i] %.17g  Δx⁻¹ %.17g  ∂x p %.17g  Gu %.17g\n", label,
                p[ii-1, j, k], p[ii, j, k], Δx⁻¹ᶠᶜᶜ(ii, j, k, m.grid), ∂xᶠᶜᶜ(ii, j, k, m.grid, p), m.timestepper.Gⁿ.u[ii, j, k])
    end
end
