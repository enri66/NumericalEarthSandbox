# Included by distributed_sea_ice.jl when EVPPROBE=true: replaces ClimaSeaIce's split-explicit `time_step_momentum!`
# with a copy that snapshots the EVP fields (u, v, h, ℵ, the rheology auxiliaries, the external stresses and η, on the
# halo-extended velocity grid) at the start of each momentum step, after each stress computation and after each
# velocity update. `compare_evp_snapshots` then compares the distributed snapshots with the serial ones at the same
# global indices, interior and halos (the extended halos included), and prints where each field first departs.
#
# EVPFIX (fixes joined with "+") = stress refreshes the halo-extended copies of the external stresses (the ocean surface
# velocities in the ice-ocean drag) and of η once the sea ice is built, before the coupled model computes its initial
# fluxes: on a distributed grid NumericalEarth's initial ice-ocean stress otherwise reads their unfilled halos, while
# serially they are the live ocean fields. corners also fills, on the distributed grid, the wall (Bounded-side) halos of these fields across the
# connected direction's halos as well, which `fill_halo_regions!` leaves out: its south/north fill spans i ∈ 1:Nx only,
# so on a 2 × 1 partition the corner blocks (rank-edge halo columns × wall halo rows) are never refreshed while the
# substeps update the halo columns.
using Oceananigans.Fields: Field, instantiated_location
using Oceananigans.Grids: topology, RightConnected, LeftConnected, FullyConnected
using Oceananigans.BoundaryConditions: construct_boundary_conditions_kernels
using Oceananigans.DistributedComputations: Distributed
using ClimaSeaIce

const EVP_SNAPSHOTS = Dict{Symbol, Vector{Any}}(:distributed => Any[], :serial => Any[])
const EVPFIX = Set(split(get(ENV, "EVPFIX", "none"), "+"))

is_field(f) = f isa Field
connected(T) = T <: Union{RightConnected, LeftConnected, FullyConnected}

# Indices spanning the whole parent array along the connected directions, so that the fill kernels of the other
# (Bounded) sides cover the corners of the extended halo
function extended_indices(f)
    TX, TY, _ = topology(f.grid)
    ix = connected(TX) ? (first(axes(f.data, 1)):last(axes(f.data, 1))) : Colon()
    iy = connected(TY) ? (first(axes(f.data, 2)):last(axes(f.data, 2))) : Colon()
    return (ix, iy, Colon())
end

extended_bcs(f) = construct_boundary_conditions_kernels(f.boundary_conditions, f.data, f.grid, instantiated_location(f), extended_indices(f))

fill_corners_locally!(f) = fill_halo_regions!(f.data, extended_bcs(f), f.indices, instantiated_location(f), f.grid,
                                              f.communication_buffers; only_local_halos = true)

function refresh_external_stresses!(model)
    dynamics = model.dynamics
    grid = model.velocities.u.grid
    ClimaSeaIce.SeaIceDynamics.update_external_stress!(dynamics.external_momentum_stresses.top, grid)
    ClimaSeaIce.SeaIceDynamics.update_external_stress!(dynamics.external_momentum_stresses.bottom, grid)
    ClimaSeaIce.SeaIceDynamics.update_free_surface!(dynamics.free_surface)
    return nothing
end

function evp_fields(model, dynamics)
    pairs_ = Pair{Symbol, Any}[]
    for (k, f) in pairs(dynamics.auxiliaries.fields)
        push!(pairs_, k => f)
    end
    push!(pairs_, :u => model.velocities.u, :v => model.velocities.v, :h => model.ice_thickness, :ℵ => model.ice_concentration)
    for (side, s) in pairs(dynamics.external_momentum_stresses)
        for name in (:u, :v, :uₑ, :vₑ)
            hasproperty(s, name) && is_field(getproperty(s, name)) && push!(pairs_, Symbol(side, "_", name) => getproperty(s, name))
        end
    end
    fs = dynamics.free_surface
    hasproperty(fs, :η) && is_field(fs.η) && push!(pairs_, :η => fs.η)
    return [p for p in pairs_ if is_field(p.second)]
end

function evp_snapshot!(label, model, dynamics)
    key = model.grid.architecture isa Distributed ? :distributed : :serial
    push!(EVP_SNAPSHOTS[key], (label, Dict(k => copy(f.data) for (k, f) in evp_fields(model, dynamics))))
    return nothing
end

function fix_corners_before_substeps!(model, dynamics)
    model.grid.architecture isa Distributed || return nothing
    for (k, f) in evp_fields(model, dynamics)
        k in (:u, :v, :h, :ℵ, :σ₁₁, :σ₂₂, :σ₁₂, :top_u, :top_v, :bottom_uₑ, :bottom_vₑ, :η) && fill_corners_locally!(f)
    end
    return nothing
end

@eval ClimaSeaIce.SeaIceDynamics function time_step_momentum!(model, dynamics::SplitExplicitMomentumEquation, Δt)
    grid = model.velocities.u.grid
    arch = architecture(grid)
    rheology      = dynamics.rheology
    u, v          = model.velocities
    free_drift    = dynamics.free_drift
    clock         = model.clock
    coriolis      = dynamics.coriolis
    massmin       = dynamics.minimum_mass
    ℵmin          = dynamics.minimum_concentration
    u_forcing     = model.forcing.u
    v_forcing     = model.forcing.v
    u_immersed_bc = u.boundary_conditions.immersed
    v_immersed_bc = v.boundary_conditions.immersed
    top_stress    = dynamics.external_momentum_stresses.top
    basal_stress  = dynamics.basal_stress
    bottom_stress = dynamics.external_momentum_stresses.bottom
    model_fields  = merge(dynamics.auxiliaries.fields, model.velocities,
                          (; h = model.ice_thickness, ℵ = model.ice_concentration, ρ = model.sea_ice_density,
                             free_surface = dynamics.free_surface))

    fix = "corners" in Main.EVPFIX && arch isa Distributed

    reset_velocities!(u, v, model.timestepper)
    update_external_stress!(top_stress, grid)
    update_external_stress!(bottom_stress, grid)
    update_free_surface!(dynamics.free_surface)
    fix && Main.fix_corners_before_substeps!(model, dynamics)
    initialize_rheology!(model, dynamics.rheology)
    Main.evp_snapshot!("it$(clock.iteration) start", model, dynamics)

    params = dynamics.solver.kernel_parameters
    u_velocity_kernel!, _ = configure_kernel(arch, grid, params, _u_velocity_step!)
    v_velocity_kernel!, _ = configure_kernel(arch, grid, params, _v_velocity_step!)
    substeps = dynamics.solver.substeps

    u_args = (u, grid, Δt, substeps, rheology, model_fields, free_drift, clock, coriolis, massmin, ℵmin, u_immersed_bc, top_stress, bottom_stress, basal_stress, u_forcing)
    v_args = (v, grid, Δt, substeps, rheology, model_fields, free_drift, clock, coriolis, massmin, ℵmin, v_immersed_bc, top_stress, bottom_stress, basal_stress, v_forcing)
    u_bcs = fix ? Main.extended_bcs(u) : u.boundary_conditions
    v_bcs = fix ? Main.extended_bcs(v) : v.boundary_conditions
    u_halo = (u.data, u_bcs, u.indices, instantiated_location(u), grid, u.communication_buffers)
    v_halo = (v.data, v_bcs, v.indices, instantiated_location(v), grid, v.communication_buffers)
    stresses_args = (model_fields, grid, rheology, Δt, u_immersed_bc, v_immersed_bc)

    fill_halo_regions!(u_halo...; only_local_halos = true)
    fill_halo_regions!(v_halo...; only_local_halos = true)

    for substep in 1:substeps
        compute_stresses!(dynamics, stresses_args...)
        Main.evp_snapshot!("it$(clock.iteration) s$substep σ", model, dynamics)
        if iseven(substep)
            u_velocity_kernel!(u_args...); fill_halo_regions!(u_halo...; only_local_halos = true)
            Main.evp_snapshot!("it$(clock.iteration) s$substep u", model, dynamics)
            v_velocity_kernel!(v_args...); fill_halo_regions!(v_halo...; only_local_halos = true)
            Main.evp_snapshot!("it$(clock.iteration) s$substep v", model, dynamics)
        else
            v_velocity_kernel!(v_args...); fill_halo_regions!(v_halo...; only_local_halos = true)
            Main.evp_snapshot!("it$(clock.iteration) s$substep v", model, dynamics)
            u_velocity_kernel!(u_args...); fill_halo_regions!(u_halo...; only_local_halos = true)
            Main.evp_snapshot!("it$(clock.iteration) s$substep u", model, dynamics)
        end
    end

    finalize_rheology!(model_fields, rheology)
    return nothing
end

# Distributed snapshots against the serial ones at the same global indices. For each snapshot and field: the largest
# interior difference, and the halo difference nearest the rank's interior (distance in cells outside it).
function compare_evp_snapshots(grid, arch, rank, comm)
    dist, ser = EVP_SNAPSHOTS[:distributed], EVP_SNAPSHOTS[:serial]
    if length(dist) != length(ser)
        rank == 0 && println("EVP probe: $(length(dist)) distributed vs $(length(ser)) serial snapshots")
        return nothing
    end
    nx, ny = size(grid, 1), size(grid, 2)
    ri, rj, _ = arch.local_index
    ox, oy = (ri - 1) * nx, (rj - 1) * ny
    first_bad = Dict{Symbol, String}()
    lines = String[]
    for ((label, dp), (label_s, sp)) in zip(dist, ser)
        for name in sort(collect(keys(dp)))
            haskey(sp, name) || continue
            p, s = dp[name], sp[name]
            imax, iat, hnear, hat, hval = 0.0, (0, 0), typemax(Int), (0, 0), 0.0
            for j in axes(p, 2), i in axes(p, 1)
                gi, gj = i + ox, j + oy
                checkbounds(Bool, s, gi, gj, 1) || continue
                d = abs(p[i, j, 1] - s[gi, gj, 1])
                (isnan(d) || d > 0) || continue
                d = isnan(d) ? Inf : d
                inside = 1 <= i <= nx && 1 <= j <= ny
                if inside
                    d > imax && ((imax, iat) = (d, (gi, gj)))
                else
                    dist_out = max(1 - i, i - nx, 1 - j, j - ny, 0)
                    dist_out < hnear && ((hnear, hat, hval) = (dist_out, (gi, gj), d))
                end
            end
            (imax > 0 || hnear < typemax(Int)) || continue
            msg = @sprintf("rank %d %-16s %-6s interior %.2e at %-10s nearest halo diff %s", rank, label, name, imax,
                           imax > 0 ? string(iat) : "-", hnear < typemax(Int) ? @sprintf("%d out at %s (%.2e)", hnear, hat, hval) : "-")
            push!(lines, msg)
            imax > 0 && !haskey(first_bad, name) && (first_bad[name] = msg)
        end
    end
    for r in 0:MPI.Comm_size(comm)-1
        if r == rank
            println("EVP probe rank $rank: first interior departure per field")
            for msg in sort(collect(values(first_bad)), by = m -> findfirst(==(m), lines))
                println("  ", msg)
            end
            println("EVP probe rank $rank: first 60 lines with any difference")
            foreach(m -> println("  ", m), lines[1:min(end, 60)])
        end
        MPI.Barrier(comm)
    end
    return nothing
end
