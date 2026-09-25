# ==================================================================
# 04 — MAB: GLORYS open boundaries + tides + tracer reservoirs + de-tided output
#
# Combines, for the first time, everything this sandbox carries separately:
#   - 02_mab_glorys_obc.jl's GLORYS-driven open boundaries (real subtidal Uᵉˣᵗ/ηᵉˣᵗ,
#     including today's fix — see that script and NumericalEarth/CLAUDE.md for the
#     Flather Uᵉˣᵗ bug history).
#   - 03_mab_m2_tide_smoke.jl's tidal machinery (`earth_tidal_harmonics`, TPXO10Atlas,
#     `tidal_forcing` body force) — but ADDED to the GLORYS subtidal signal rather than
#     replacing it: `tidal_boundary_conditions` builds standalone (U,V,η) BCs on its own,
#     which would just overwrite the GLORYS ones, so this script instead computes the
#     same per-side Flather tidal (U,η) contribution `tidal_boundary_conditions` does
#     internally and sums it into the GLORYS `U_west`/`U_east`/`V_south`/`V_north`
#     discrete functions. Physically: total sea level = subtidal (GLORYS) + tide; same
#     for the barotropic transport. Tangential/tracer BCs are untouched — Flather/Chapman
#     is still the ONLY entry point for the barotropic mode (see 02's header note), so
#     adding the tide there is enough to carry it into `u`/`v` too via the split-explicit
#     corrector.
#   - `TracerReservoir` (Oceananigans PR #5964) on the T/S boundaries, replacing 02's
#     `NormalRadiation` tracer scheme.
#   - `FilteredTimeInterval` (Oceananigans PR #5971) for de-tided daily (and 5-day) SSH and
#     surface-field output, alongside raw hourly output for animation.
#
# Run:  julia -t 8 --project=. scripts/04_mab_glorys_tides_reservoirs.jl
# ==================================================================

using NumericalEarth
using NumericalEarth.DataWrangling: Metadata, MetadataSet, BoundingBox
using NumericalEarth.NestedModels: Interpolated
using CopernicusMarine              # activates the GLORYS download backend
using CopernicusClimateDataStore    # activates the ERA5 download backend
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: ExponentialDiscretization, znodes, λnodes, φnodes
using Oceananigans.BoundaryConditions: PerturbationAdvection, NormalRadiation, ObliqueRadiation, TracerReservoir,
                                       GravityWaveRadiationBoundaryCondition,
                                       SurfaceWaveRadiationBoundaryCondition
using Dates, Printf, Statistics

include(joinpath(@__DIR__, "glorys_bathymetry.jl"))

# See 02_mab_glorys_obc.jl for why this defaults to the shared out-of-Dropbox cache.
const DATA_DIR   = get(ENV, "MAB_DATA_DIR", joinpath(homedir(), "Data", "NumericalEarth"))
const resolution = 1 / 12
const Nz = 40

const n_pad      = 2
const data_λ     = (-76.0, -64.0)          # GLORYS box (what is on disk)
const data_φ     = ( 34.0,  42.0)
const λ_bounds   = (data_λ[1] + n_pad*resolution, data_λ[2] - n_pad*resolution)   # model
const φ_bounds   = (data_φ[1] + n_pad*resolution, data_φ[2] - n_pad*resolution)
const Nλ = round(Int, (λ_bounds[2] - λ_bounds[1]) / resolution)
const Nφ = round(Int, (φ_bounds[2] - φ_bounds[1]) / resolution)
const Nλ_src = round(Int, (data_λ[2] - data_λ[1]) / resolution)
const Nφ_src = round(Int, (data_φ[2] - data_φ[1]) / resolution)

const start_date = DateTime(2019, 4, 1)
const sim_days   = parse(Int, get(ENV, "MAB_DAYS", "14"))
const stop_date  = start_date + Day(sim_days)

# variant switches (see 02_mab_glorys_obc.jl for UEXT/TANGENTIAL — both default to the
# values that script validated: native Uᵉˣᵗ, radiated tangential)
const TANGENTIAL  = get(ENV, "MAB_TANGENTIAL", "radiated")
# 3D u/v open-boundary scheme: "oblique" (default, ObliqueRadiation on normal and tangential components)
# or "legacy" (PerturbationAdvection normal + NormalRadiation tangential, as before 2026-09-24)
const VELOCITY_SCHEME = get(ENV, "MAB_VELOCITY_SCHEME", "oblique")
# ObliqueRadiation nudging toward the exterior (GLORYS + tide) and phase-speed averaging weight
const OBLIQUE_TAU_IN  = parse(Float64, get(ENV, "MAB_OBLIQUE_TAU_IN", "3")) * days
const OBLIQUE_TAU_OUT = parse(Float64, get(ENV, "MAB_OBLIQUE_TAU_OUT", "360")) * days
const PHASE_SPEED_WEIGHT = parse(Float64, get(ENV, "MAB_PHASE_SPEED_WEIGHT", "0.3"))
const TAU_IN      = parse(Float64, get(ENV, "MAB_TAU_IN", "1")) * days
const TAG         = get(ENV, "MAB_TAG", "mab_glorys_tides")
const MATCH_BATHY = get(ENV, "MAB_MATCH_BATHY", "true") == "true"
const N_MATCH     = parse(Int, get(ENV, "MAB_N_MATCH", "4"))
const UEXT_MODE   = get(ENV, "MAB_UEXT", "native")
const SPIKE_THRESHOLD = parse(Float64, get(ENV, "MAB_SPIKE_THRESHOLD", "2.0"))

# tides
const TIDE_CONSTITUENTS = Symbol.(split(get(ENV, "MAB_TIDE_CONSTITUENTS", "M2,S2,N2,K2,K1,O1,P1,Q1,Mf,Mm"), ","))
const TIDE_RAMP  = parse(Float64, get(ENV, "MAB_TIDE_RAMP", "1")) * days
const TPXO_DIR   = get(ENV, "TPXO_DIR", joinpath(homedir(), "Data", "TPXO10_atlas_v2_nc"))

# tracer reservoirs (MOM6-scale defaults: relax over 20 km on inflow, memoryless on outflow)
const RESERVOIR_L_IN  = parse(Float64, get(ENV, "MAB_RESERVOIR_L_IN", "20000"))
const RESERVOIR_L_OUT = parse(Float64, get(ENV, "MAB_RESERVOIR_L_OUT", "0"))
# T/S open-boundary scheme: "reservoir" (default), "radiation" (NormalRadiation, as script 02) or "oblique"
const TRACER_SCHEME = get(ENV, "MAB_TRACER_SCHEME", "reservoir")
# "true" (default): make the 3D normal velocity at each open face integrate to the barotropic exterior transport
# used by the Flather condition (see ConsistentNormalFlow below)
const CONSISTENT_UBC = get(ENV, "MAB_CONSISTENT_UBC", "true") == "true"
# GLORYS restoring sponge along the open boundaries: the variables to restore ("" = no sponge, "T,S" or
# "T,S,u,v", where u and v restore only their baroclinic part), the band width in cells and the restoring
# timescale at the boundary (see sponge_masks below)
const SPONGE_VARS  = Symbol.(filter(!isempty, split(get(ENV, "MAB_SPONGE_VARS", ""), ",")))
const SPONGE_WIDTH = parse(Int, get(ENV, "MAB_SPONGE_WIDTH", "8"))
const SPONGE_TAU   = parse(Float64, get(ENV, "MAB_SPONGE_TAU", "0.25")) * days
# Surface salinity restoring toward GLORYS as a salt flux with this piston velocity (m/day; 0 turns it off). Over a
# mixed layer of depth h it damps salinity differences in h / piston: ~30 days for a 15 m summer mixed layer,
# ~120 days for 60 m, slow enough to keep the model's eddies and fast enough to hold the seasonal cycle.
const SSS_PISTON   = parse(Float64, get(ENV, "MAB_SSS_PISTON", "0.5")) / days

# checkpoint/restart — PICKUP itself is parsed further down, right before it's used, since it
# can be a Bool, an iteration number, or a filepath (see the comment there)
const CHECKPOINT_EVERY = parse(Float64, get(ENV, "MAB_CHECKPOINT_EVERY", "5")) * days

mkpath(DATA_DIR)

# ---------------- grid + bathymetry (identical to script 02) ----------------
z = ExponentialDiscretization(Nz, -4000, 0; scale = 1400)
grid = LatitudeLongitudeGrid(CPU(); size = (Nλ, Nφ, Nz),
                             longitude = λ_bounds, latitude = φ_bounds, z, halo = (7, 7, 7))
bottom_height = regrid_bathymetry(grid; dataset = ETOPO2022(), height_above_water = 1,
                                  minimum_depth = 10, major_basins = 1, interpolation_passes = 10)

if MATCH_BATHY
    match_boundary_bathymetry!(bottom_height, grid, DATA_DIR;
                               region = BoundingBox(longitude = data_λ, latitude = data_φ),
                               n_match = N_MATCH, minimum_depth = 10)
end

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))
@show grid

# ---------------- GLORYS as FieldTimeSeries on the model grid ----------------
region = BoundingBox(longitude = data_λ, latitude = data_φ)

src_grid = LatitudeLongitudeGrid(CPU(); size = (Nλ_src, Nφ_src, Nz),
                                 longitude = data_λ, latitude = data_φ, z, halo = (7, 7, 7))
dates  = start_date : Day(1) : stop_date
glorys = GLORYSDaily()
meta(name) = Metadata(name; dataset = glorys, dates, dir = DATA_DIR, region)

@info "building GLORYS FieldTimeSeries on the model grid (downloads/inpaints as needed)…"
fts_u = FieldTimeSeries(meta(:u_velocity),  src_grid)
fts_v = FieldTimeSeries(meta(:v_velocity),  src_grid)
fts_T = FieldTimeSeries(meta(:temperature), src_grid)
fts_S = FieldTimeSeries(meta(:salinity),    src_grid)
fts_η = FieldTimeSeries(meta(:free_surface), src_grid)
@info "  done: $(length(fts_u.times)) times, $(Dates.format(start_date,"yyyy-mm-dd")) → $(Dates.format(stop_date,"yyyy-mm-dd"))"

# ---------------- the SUBTIDAL barotropic exterior, from GLORYS (see 02_mab_glorys_obc.jl) ----------------
const zc_src = collect(znodes(src_grid, Center()))
const Δz_src = diff(collect(znodes(src_grid, Face())))
const bh     = Array(interior(bottom_height))[:, :, 1]
const floors = ((bh[1, :], bh[end, :]), (bh[:, 1], bh[:, end]))

function wet_transport(col, floor_line)
    Ns = size(col, 1)
    T, Hwet = zeros(Ns), zeros(Ns)
    for s in 1:Ns
        b = floor_line[clamp(s - n_pad, 1, length(floor_line))]
        for k in 1:size(col, 2)
            if zc_src[k] > b
                T[s]    += Δz_src[k] * col[s, k]
                Hwet[s] += Δz_src[k]
            end
        end
    end
    return T, Hwet
end

function boundary_columns(A, dim; Nz_expected = Nz)
    size(A) == (Nλ_src, Nφ_src, Nz_expected) ||
        error("unexpected GLORYS FieldTimeSeries frame size $(size(A)); expected $((Nλ_src, Nφ_src, Nz_expected))")
    p = n_pad + 1
    dim == 1 ? (A[p, :, :], A[end-n_pad, :, :]) : (A[:, p, :], A[:, end-n_pad, :])
end

function boundary_transport_slabs(fts, dim)
    times = collect(fts.times)
    slabs = map(eachindex(times)) do n
        lo, hi = boundary_columns(Array(interior(fts[n])), dim)
        (wet_transport(lo, floors[dim][1])[1], wet_transport(hi, floors[dim][2])[1])
    end
    return times, first.(slabs), last.(slabs)
end

z_native = NumericalEarth.DataWrangling.z_interfaces(meta(:u_velocity))
Nz_native = length(z_native) - 1
src_grid_native_z = LatitudeLongitudeGrid(CPU(); size = (Nλ_src, Nφ_src, Nz_native),
                                          longitude = data_λ, latitude = data_φ,
                                          z = z_native, halo = (7, 7, 7))
@info "building the native-resolution GLORYS u/v for the true depth mean… ($Nz_native levels vs the model's $Nz)"
fts_u_native = FieldTimeSeries(meta(:u_velocity), src_grid_native_z)
fts_v_native = FieldTimeSeries(meta(:v_velocity), src_grid_native_z)
Hg_src = glorys_deptho_on_grid(src_grid, DATA_DIR; region)
true_floors = ((-Hg_src[1, :], -Hg_src[end, :]), (-Hg_src[:, 1], -Hg_src[:, end]))

function true_depth_mean(col, zc, Δz, floor_line)
    Ns = size(col, 1)
    ubar = zeros(Ns)
    for s in 1:Ns
        b = floor_line[clamp(s - n_pad, 1, length(floor_line))]
        T = H = 0.0
        for k in 1:size(col, 2)
            v = col[s, k]
            if zc[k] > b && isfinite(v)
                T += Δz[k] * v
                H += Δz[k]
            end
        end
        ubar[s] = H > 0 ? T / H : 0.0
    end
    return ubar
end

function boundary_transport_slabs_native(fts, floor_pair, dim)
    zc = collect(znodes(fts.grid, Center()))
    Δz = diff(collect(znodes(fts.grid, Face())))
    times = collect(fts.times)
    slabs = map(eachindex(times)) do n
        lo, hi = boundary_columns(Array(interior(fts[n])), dim; Nz_expected = length(zc))
        ū_lo = true_depth_mean(lo, zc, Δz, floor_pair[1])
        ū_hi = true_depth_mean(hi, zc, Δz, floor_pair[2])
        Hwet_lo = last(wet_transport(zeros(size(lo, 1), Nz), floors[dim][1]))
        Hwet_hi = last(wet_transport(zeros(size(hi, 1), Nz), floors[dim][2]))
        (ū_lo .* Hwet_lo, ū_hi .* Hwet_hi)
    end
    return times, first.(slabs), last.(slabs)
end

if UEXT_MODE == "native"
    tU, U_lo, U_hi = boundary_transport_slabs_native(fts_u_native, true_floors[1], 1)
    tV, V_lo, V_hi = boundary_transport_slabs_native(fts_v_native, true_floors[2], 2)
else
    tU, U_lo, U_hi = boundary_transport_slabs(fts_u, 1)
    tV, V_lo, V_hi = boundary_transport_slabs(fts_v, 2)
end

η_slabs = map(eachindex(fts_η.times)) do n
    A = Array(interior(fts_η[n]))[:, :, 1]
    p = n_pad + 1
    (A[p, :], A[end-n_pad, :], A[:, p], A[:, end-n_pad])
end

@inline function frame(times, t)
    t <= times[1]   && return (1, 1, 0.0)
    t >= times[end] && return (length(times), length(times), 0.0)
    n = searchsortedlast(times, t)
    return (n, n + 1, (t - times[n]) / (times[n+1] - times[n]))
end
@inline function lerp_slab(times, slabs, i, t)
    n1, n2, w = frame(times, t)
    a = slabs[n1]; b = slabs[n2]
    ii = clamp(i + n_pad, 1, length(a))
    return (1 - w) * a[ii] + w * b[ii]
end
ηs(k) = [s[k] for s in η_slabs]
η_w, η_e, η_s, η_n = ηs(1), ηs(2), ηs(3), ηs(4)

# ---------------- the TIDAL barotropic exterior, from TPXO10 ----------------
# Mirrors `tidal_boundary_conditions`'s own node choice and `tidal_transport_and_elevation`'s
# math exactly (see Oceananigans.BoundaryConditions), so that ADDING this to the GLORYS
# subtidal (U,η) above is equivalent to what `tidal_boundary_conditions` would compute if it
# were driving the barotropic mode alone — just summed with a second (subtidal) contribution
# instead of being the only one. `tidal_boundary_conditions` itself is not called here for
# exactly that reason: it returns a complete, standalone (U,V,η), which would overwrite the
# GLORYS one rather than add to it.
harmonics = earth_tidal_harmonics(start_date; constituents = TIDE_CONSTITUENTS, ramp_time = TIDE_RAMP)
@info "tides: $(harmonics)"

λᶜ, λᶠ = λnodes(grid, Center()), λnodes(grid, Face())
φᶜ, φᶠ = φnodes(grid, Center()), φnodes(grid, Face())

tidal_constants(field, nodes) =
    reduce(hcat, getproperty.(map(name -> tidal_atlas_constants(TPXO10Atlas(), nodes, name; dir = TPXO_DIR),
                                  harmonics.constituents), field))

west_U_const  = tidal_constants(:eastward_transport,  [(first(λᶠ), φ) for φ in φᶜ])
west_η_const  = tidal_constants(:sea_surface_height,  [(first(λᶜ), φ) for φ in φᶜ])
east_U_const  = tidal_constants(:eastward_transport,  [(last(λᶠ),  φ) for φ in φᶜ])
east_η_const  = tidal_constants(:sea_surface_height,  [(last(λᶜ),  φ) for φ in φᶜ])
south_V_const = tidal_constants(:northward_transport, [(λ, first(φᶠ)) for λ in λᶜ])
south_η_const = tidal_constants(:sea_surface_height,  [(λ, first(φᶜ)) for λ in λᶜ])
north_V_const = tidal_constants(:northward_transport, [(λ, last(φᶠ))  for λ in λᶜ])
north_η_const = tidal_constants(:sea_surface_height,  [(λ, last(φᶜ))  for λ in λᶜ])

@inline function tidal_UV_eta(Uconst, ηconst, harmonics, t, i)
    U = η = 0.0
    @inbounds for n in eachindex(harmonics.frequencies)
        rot = harmonics.nodal_factors[n] * cis(harmonics.frequencies[n] * t + harmonics.phases[n])
        U += real(Uconst[i, n] * rot)
        η += real(ηconst[i, n] * rot)
    end
    ramp = harmonics.ramp_time > 0 ? tanh(t / harmonics.ramp_time) : 1.0
    return ramp * U, ramp * η
end

function U_west(j, k, grid, clock, f)
    Us, ηsub = lerp_slab(tU, U_lo, j, clock.time), lerp_slab(tU, η_w, j, clock.time)
    Ut, ηt = tidal_UV_eta(west_U_const, west_η_const, harmonics, clock.time, j)
    return (Us + Ut, ηsub + ηt)
end
function U_east(j, k, grid, clock, f)
    Us, ηsub = lerp_slab(tU, U_hi, j, clock.time), lerp_slab(tU, η_e, j, clock.time)
    Ut, ηt = tidal_UV_eta(east_U_const, east_η_const, harmonics, clock.time, j)
    return (Us + Ut, ηsub + ηt)
end
function V_south(i, k, grid, clock, f)
    Vs, ηsub = lerp_slab(tV, V_lo, i, clock.time), lerp_slab(tV, η_s, i, clock.time)
    Vt, ηt = tidal_UV_eta(south_V_const, south_η_const, harmonics, clock.time, i)
    return (Vs + Vt, ηsub + ηt)
end
function V_north(i, k, grid, clock, f)
    Vs, ηsub = lerp_slab(tV, V_hi, i, clock.time), lerp_slab(tV, η_n, i, clock.time)
    Vt, ηt = tidal_UV_eta(north_V_const, north_η_const, harmonics, clock.time, i)
    return (Vs + Vt, ηsub + ηt)
end

# ---------------- 3D normal velocity consistent with the barotropic exterior (MAB_CONSISTENT_UBC) ----------------
# The Flather condition drives the barotropic transport toward `U_west/U_east/V_south/V_north` (GLORYS
# depth mean times the model's wet depth, plus the tide), while the 3D normal velocity is the GLORYS
# profile interpolated onto the model levels, without the tide. Their depth integrals differ, and in the
# last interior column that difference appears as a spurious surface vertical velocity. This wrapper adds a
# depth-uniform velocity to the interpolated profile so that its integral over the wet cells equals the
# barotropic target exactly.
import Oceananigans.BoundaryConditions: regularize_boundary_condition, getbc

struct ConsistentNormalFlow{B, U, T}
    base   :: B   # the Interpolated GLORYS velocity
    Utotal :: U   # (s, k, grid, clock, fields) -> (barotropic transport, η), as for the Flather functions
    table  :: T   # times, depth integral of the interpolated profile per time and boundary point, wet depth, floor
end

regularize_boundary_condition(c::ConsistentNormalFlow, grid, loc, dim, Side, args...) =
    ConsistentNormalFlow(regularize_boundary_condition(c.base, grid, loc, dim, Side, args...), c.Utotal, c.table)

@inline function getbc(c::ConsistentNormalFlow, s::Integer, k::Integer, grid::Oceananigans.Grids.AbstractGrid, clock = nothing, args...)
    u = getbc(c.base, s, k, grid, clock, args...)
    tb = c.table
    1 <= s <= length(tb.floor) || return u
    Oceananigans.Grids.znode(1, 1, k, grid, Center(), Center(), Center()) > tb.floor[s] || return u
    t = isnothing(clock) ? 0.0 : clock.time
    n1, n2, w = frame(tb.times, t)
    Uint = (1 - w) * tb.Uint[n1, s] + w * tb.Uint[n2, s]
    return u + (c.Utotal(s, k, grid, (; time = t), nothing)[1] - Uint) / tb.Hwet[s]
end

function make_consistency_table(floor_line)
    Hwet = [sum(Δz_src[k] for k in 1:Nz if zc_src[k] > b; init = 0.0) for b in floor_line]
    return (; times = collect(fts_u.times), Uint = zeros(length(fts_u.times), length(floor_line)),
            Hwet, floor = collect(floor_line))
end

consistency_tables = (u_west  = make_consistency_table(floors[1][1]), u_east  = make_consistency_table(floors[1][2]),
                      v_south = make_consistency_table(floors[2][1]), v_north = make_consistency_table(floors[2][2]))

normal_condition(fts, Ufun, table) = CONSISTENT_UBC ? ConsistentNormalFlow(Interpolated(fts), Ufun, table) : Interpolated(fts)

# ---------------- the boundary conditions ----------------
oblique_scheme = ObliqueRadiation(inflow_timescale = OBLIQUE_TAU_IN, outflow_timescale = OBLIQUE_TAU_OUT,
                                  phase_speed_weight = PHASE_SPEED_WEIGHT)
normal_scheme  = VELOCITY_SCHEME == "oblique" ? oblique_scheme :
                 VELOCITY_SCHEME == "legacy"  ? PerturbationAdvection(inflow_timescale = TAU_IN, outflow_timescale = Inf) :
                 error("MAB_VELOCITY_SCHEME must be oblique or legacy, got $VELOCITY_SCHEME")
# MAB_TANGENTIAL: "radiated" (the velocity scheme's own radiation), "normal" (NormalRadiation with the
# 1-day inflow nudging, whatever the velocity scheme), "oblique" (ObliqueRadiation, whatever the velocity
# scheme), or "prescribed" (the GLORYS value, no radiation)
tangential_sch = TANGENTIAL == "oblique" ? oblique_scheme :
                 TANGENTIAL == "normal" || VELOCITY_SCHEME != "oblique" ?
                 NormalRadiation(inflow_timescale = TAU_IN, outflow_timescale = Inf) : oblique_scheme
tracer_scheme  = TRACER_SCHEME == "reservoir" ? TracerReservoir(inflow_length_scale = RESERVOIR_L_IN, outflow_length_scale = RESERVOIR_L_OUT) :
                 TRACER_SCHEME == "radiation" ? NormalRadiation(inflow_timescale = TAU_IN, outflow_timescale = Inf) :
                 TRACER_SCHEME == "oblique"   ? ObliqueRadiation(inflow_timescale = TAU_IN, outflow_timescale = Inf) :
                 error("MAB_TRACER_SCHEME must be reservoir, radiation or oblique, got $TRACER_SCHEME")

tangential(fts) = TANGENTIAL == "prescribed" ?
    ValueBoundaryCondition(Interpolated(fts)) :
    ValueBoundaryCondition(Interpolated(fts); scheme = tangential_sch)

u_bcs = FieldBoundaryConditions(
    west  = NormalFlowBoundaryCondition(normal_condition(fts_u, U_west, consistency_tables.u_west); scheme = normal_scheme),
    east  = NormalFlowBoundaryCondition(normal_condition(fts_u, U_east, consistency_tables.u_east); scheme = normal_scheme),
    south = tangential(fts_u),
    north = tangential(fts_u))

v_bcs = FieldBoundaryConditions(
    south = NormalFlowBoundaryCondition(normal_condition(fts_v, V_south, consistency_tables.v_south); scheme = normal_scheme),
    north = NormalFlowBoundaryCondition(normal_condition(fts_v, V_north, consistency_tables.v_north); scheme = normal_scheme),
    west  = tangential(fts_v),
    east  = tangential(fts_v))

tracer_bcs(fts) = FieldBoundaryConditions(
    west  = ValueBoundaryCondition(Interpolated(fts); scheme = tracer_scheme),
    east  = ValueBoundaryCondition(Interpolated(fts); scheme = tracer_scheme),
    south = ValueBoundaryCondition(Interpolated(fts); scheme = tracer_scheme),
    north = ValueBoundaryCondition(Interpolated(fts); scheme = tracer_scheme))

U_bcs = FieldBoundaryConditions(grid, (Face(), Center(), nothing);
    west = GravityWaveRadiationBoundaryCondition(U_west; discrete_form = true),
    east = GravityWaveRadiationBoundaryCondition(U_east; discrete_form = true))
V_bcs = FieldBoundaryConditions(grid, (Center(), Face(), nothing);
    south = GravityWaveRadiationBoundaryCondition(V_south; discrete_form = true),
    north = GravityWaveRadiationBoundaryCondition(V_north; discrete_form = true))
η_bcs = FieldBoundaryConditions(grid, (Center(), Center(), Face());
    west = SurfaceWaveRadiationBoundaryCondition(), east = SurfaceWaveRadiationBoundaryCondition(),
    south = SurfaceWaveRadiationBoundaryCondition(), north = SurfaceWaveRadiationBoundaryCondition())

boundary_conditions = (u = u_bcs, v = v_bcs, T = tracer_bcs(fts_T), S = tracer_bcs(fts_S),
                       U = U_bcs, V = V_bcs, η = η_bcs)

# ---------------- GLORYS restoring sponge (MAB_SPONGE_VARS) ----------------
# The mask falls from ≈1 at an open boundary to 0 at SPONGE_WIDTH cells from it as cos²(π d / 2W), where d is
# the distance in cells from the boundary face. A cell belongs to a side's band only if it is connected to that
# side's wet boundary cell along its row (west/east) or column (south/north), so bays and sounds behind land
# are left alone. Near a corner the mask is the larger of the two sides' values, i.e. it follows the smaller
# distance, so the seam between the two bands is the diagonal into the corner.
function sponge_masks(wet, W)
    Nx, Ny = size(wet)
    d = fill(Inf, Nx, Ny)
    for j in 1:Ny
        for i in 1:min(W, Nx)
            wet[i, j] || break
            d[i, j] = min(d[i, j], i - 0.5)
        end
        for i in Nx:-1:max(Nx - W + 1, 1)
            wet[i, j] || break
            d[i, j] = min(d[i, j], Nx - i + 0.5)
        end
    end
    for i in 1:Nx
        for j in 1:min(W, Ny)
            wet[i, j] || break
            d[i, j] = min(d[i, j], j - 0.5)
        end
        for j in Ny:-1:max(Ny - W + 1, 1)
            wet[i, j] || break
            d[i, j] = min(d[i, j], Ny - j + 0.5)
        end
    end
    μᶜ = [d[i, j] < W ? cos(π * d[i, j] / 2W)^2 : 0.0 for i in 1:Nx, j in 1:Ny]
    # u and v points: the mean of the two neighbouring cell centres (the edge value at the boundary faces)
    μᵘ = [(μᶜ[clamp(i - 1, 1, Nx), j] + μᶜ[clamp(i, 1, Nx), j]) / 2 for i in 1:Nx+1, j in 1:Ny]
    μᵛ = [(μᶜ[i, clamp(j - 1, 1, Ny)] + μᶜ[i, clamp(j, 1, Ny)]) / 2 for i in 1:Nx, j in 1:Ny+1]
    return (; μᶜ, μᵘ, μᵛ)
end

const sponge = sponge_masks(bh .< 0, SPONGE_WIDTH)

# The sponge restores toward the same GLORYS series the boundary conditions use. GLORYS is interpolated onto the
# model's sponge cells once per GLORYS frame, before the run; during the run the target is only interpolated in time.
# T and S restore as they are: r μ (ψᴳ - ψ). u and v restore only their baroclinic part: r μ [(ψᴳ - ψ̄ᴳ) - (ψ - ψ̄)],
# where the overbar is the mean over the model's wet column. That integrates to zero over the column, so the
# barotropic flow (and the tide) is left to the Flather condition. The model's depth mean is computed in the
# forcing itself, so it is always the current one.
using Oceananigans.Grids: inactive_cell

@inline wet_node(i, j, k, grid, ::Center, ::Center) = !inactive_cell(i, j, k, grid)
@inline wet_node(i, j, k, grid, ::Face, ::Center) = !inactive_cell(i, j, k, grid) & !inactive_cell(i - 1, j, k, grid)
@inline wet_node(i, j, k, grid, ::Center, ::Face) = !inactive_cell(i, j, k, grid) & !inactive_cell(i, j - 1, k, grid)

function sponge_parameters(fts, μ, LX, LY; baroclinic)
    columns = [(i, j) for j in axes(μ, 2), i in axes(μ, 1) if μ[i, j] > 0]
    column_index = zeros(Int, size(μ))
    weights = zeros(length(columns), Nz)
    for (c, (i, j)) in enumerate(columns)
        column_index[i, j] = c
        wet = [wet_node(i, j, k, grid, LX, LY) for k in 1:Nz]
        H = sum(Δz_src[wet]; init = 0.0)
        H > 0 && (weights[c, :] .= wet .* Δz_src ./ H)
    end

    times = collect(fts.times)
    target = zeros(length(columns), Nz, length(times))
    loc = Oceananigans.instantiated_location(fts)
    profile = zeros(Nz)
    for n in eachindex(times)
        frame_field = fts[n]
        for (c, (i, j)) in enumerate(columns)
            for k in 1:Nz
                profile[k] = weights[c, k] > 0 ?
                    Oceananigans.Fields.interpolate(Oceananigans.Grids.node(i, j, k, grid, LX, LY, Center()), frame_field, loc, fts.grid) : 0.0
            end
            mean_profile = baroclinic ? sum(weights[c, k] * profile[k] for k in 1:Nz) : 0.0
            for k in 1:Nz
                target[c, k, n] = weights[c, k] > 0 ? profile[k] - mean_profile : 0.0
            end
        end
    end

    return (; column_index, μ = [μ[i, j] for (i, j) in columns], weights, target, times, rate = 1 / SPONGE_TAU, baroclinic)
end

@inline function sponge_restoring(i, j, k, grid, clock, ψ, p)
    ci = clamp(i, 1, size(p.column_index, 1))
    cj = clamp(j, 1, size(p.column_index, 2))
    c = @inbounds p.column_index[ci, cj]
    (c == 0 || @inbounds(p.weights[c, k]) == 0) && return zero(eltype(grid))
    ψ̄ = zero(eltype(grid))
    if p.baroclinic
        @inbounds for k′ in axes(p.weights, 2)
            ψ̄ += p.weights[c, k′] * ψ[i, j, k′]
        end
    end
    n₁, n₂, w = frame(p.times, clock.time)
    ψᴳ = @inbounds (1 - w) * p.target[c, k, n₁] + w * p.target[c, k, n₂]
    return @inbounds p.rate * p.μ[c] * (ψᴳ - (ψ[i, j, k] - ψ̄))
end

@inline T_sponge(i, j, k, grid, clock, fields, p) = sponge_restoring(i, j, k, grid, clock, fields.T, p)
@inline S_sponge(i, j, k, grid, clock, fields, p) = sponge_restoring(i, j, k, grid, clock, fields.S, p)
@inline u_sponge(i, j, k, grid, clock, fields, p) = sponge_restoring(i, j, k, grid, clock, fields.u, p)
@inline v_sponge(i, j, k, grid, clock, fields, p) = sponge_restoring(i, j, k, grid, clock, fields.v, p)

function sponge_forcing_for(name)
    name == :T && return Forcing(T_sponge; discrete_form = true,
                                 parameters = sponge_parameters(fts_T, sponge.μᶜ, Center(), Center(); baroclinic = false))
    name == :S && return Forcing(S_sponge; discrete_form = true,
                                 parameters = sponge_parameters(fts_S, sponge.μᶜ, Center(), Center(); baroclinic = false))
    name == :u && return Forcing(u_sponge; discrete_form = true,
                                 parameters = sponge_parameters(fts_u, sponge.μᵘ, Face(), Center(); baroclinic = true))
    name == :v && return Forcing(v_sponge; discrete_form = true,
                                 parameters = sponge_parameters(fts_v, sponge.μᵛ, Center(), Face(); baroclinic = true))
    error("MAB_SPONGE_VARS entries must be T, S, u or v, got $name")
end

sponge_forcing = NamedTuple(name => sponge_forcing_for(name) for name in SPONGE_VARS)

tides   = tidal_forcing(harmonics)
forcing = merge(sponge_forcing,
                (u = :u in SPONGE_VARS ? (tides.u, sponge_forcing.u) : tides.u,
                 v = :v in SPONGE_VARS ? (tides.v, sponge_forcing.v) : tides.v))
isempty(SPONGE_VARS) || @info @sprintf("GLORYS sponge on %s: %d cells, τ = %.2f days at the boundary, %d wet cells with μ > 0.01",
                                       join(SPONGE_VARS, ","), SPONGE_WIDTH, SPONGE_TAU / days, count(>(0.01), sponge.μᶜ))

# ---------------- surface salinity restoring toward GLORYS (MAB_SSS_PISTON) ----------------
# A surface salt flux J = -w (Sᴳ - S) on the top cell, with w the piston velocity, added to the bulk-formula fluxes.
# GLORYS salinity is interpolated onto the model's top cells once per GLORYS frame, before the run.
struct SurfaceSalinityRestoring{A, T}
    target :: A
    times  :: T
    piston :: Float64
end

function SurfaceSalinityRestoring(fts, piston)
    times = collect(fts.times)
    loc = Oceananigans.instantiated_location(fts)
    target = zeros(Nλ, Nφ, length(times))
    for n in eachindex(times)
        frame_field = fts[n]
        for j in 1:Nφ, i in 1:Nλ
            X = Oceananigans.Grids.node(i, j, Nz, grid, Center(), Center(), Center())
            target[i, j, n] = Oceananigans.Fields.interpolate(X, frame_field, loc, fts.grid)
        end
    end
    return SurfaceSalinityRestoring(target, times, piston)
end

@inline function getbc(r::SurfaceSalinityRestoring, i::Integer, j::Integer, grid::Oceananigans.Grids.AbstractGrid, clock, fields)
    ci = clamp(i, 1, size(r.target, 1))
    cj = clamp(j, 1, size(r.target, 2))
    n₁, n₂, w = frame(r.times, clock.time)
    Sᴳ = @inbounds (1 - w) * r.target[ci, cj, n₁] + w * r.target[ci, cj, n₂]
    return @inbounds - r.piston * (Sᴳ - fields.S[i, j, grid.Nz])
end

sss_restoring = SSS_PISTON > 0 ? SurfaceSalinityRestoring(fts_S, SSS_PISTON) : nothing
additional_surface_fluxes = isnothing(sss_restoring) ? NamedTuple() : (; S = sss_restoring)
isnothing(sss_restoring) || @info @sprintf("surface salinity restoring toward GLORYS: piston velocity %.2f m/day", SSS_PISTON * days)

# ---------------- ocean simulation: GLORYS boundaries + equilibrium tidal body force ----------------
ocean = ocean_simulation(grid; boundary_conditions, forcing, additional_surface_fluxes)

if CONSISTENT_UBC
    vel = ocean.model.velocities
    mgrid = ocean.model.grid
    # (boundary conditions, side, table, source series, normal dimension, boundary face index)
    faces = ((vel.u.boundary_conditions, :west,  :u_west,  fts_u, 1, 1),
             (vel.u.boundary_conditions, :east,  :u_east,  fts_u, 1, mgrid.Nx + 1),
             (vel.v.boundary_conditions, :south, :v_south, fts_v, 2, 1),
             (vel.v.boundary_conditions, :north, :v_north, fts_v, 2, mgrid.Ny + 1))
    for (bcs, side, name, fts, dim, iface) in faces
        c = getproperty(bcs, side).condition
        c isa ConsistentNormalFlow || error("expected a ConsistentNormalFlow on $name, found $(typeof(c))")
        tb = consistency_tables[name]
        # Integrate each source frame exactly as the boundary condition samples it (the same node and
        # spatial interpolation), from the frame fields themselves: querying the time series at
        # arbitrary times here would read slices that are not in memory yet.
        loc = Oceananigans.instantiated_location(fts)
        for n in eachindex(tb.times), s in eachindex(tb.Hwet)
            acc = 0.0
            for k in 1:Nz
                zc_src[k] > tb.floor[s] || continue
                X = dim == 1 ? Oceananigans.Grids.node(iface, s, k, mgrid, Face(), Center(), Center()) :
                               Oceananigans.Grids.node(s, iface, k, mgrid, Center(), Face(), Center())
                acc += Δz_src[k] * Oceananigans.Fields.interpolate(X, fts[n], loc, fts.grid)
            end
            tb.Uint[n, s] = acc
        end
        δ = [abs(c.Utotal(s, 1, mgrid, (; time = 0.0), nothing)[1] - tb.Uint[1, s]) / tb.Hwet[s]
             for s in eachindex(tb.Hwet) if tb.Hwet[s] > 0]
        @info @sprintf("consistent normal velocity, %s: depth-uniform correction at t=0 has max %.2e m/s, mean %.2e m/s",
                       name, maximum(δ), mean(δ))
    end
end

set!(ocean.model, MetadataSet(:temperature, :salinity, :u_velocity, :v_velocity;
                              dataset = glorys, date = start_date, dir = DATA_DIR, region))
# Start at GLORYS's own sea level: the Flather exterior values are GLORYS's zos, so removing the domain mean
# here would make the boundaries import or export that mean in the first few hours.
set!((; free_surface = ocean.model.free_surface.displacement),
     MetadataSet(:free_surface; dataset = glorys, date = start_date, dir = DATA_DIR, region))

atmosphere = ERA5PrescribedAtmosphere(; start_date, end_date = stop_date, region,
                                      dir = joinpath(DATA_DIR, "era5"))
radiation  = ERA5PrescribedRadiation(;  start_date, end_date = stop_date, region,
                                      dir = joinpath(DATA_DIR, "era5"))
model = OceanOnlyModel(ocean; atmosphere, radiation)
simulation = Simulation(model; Δt = 5minutes, stop_time = sim_days * days)

wall = Ref(time())
function progress(sim)
    u, v, w = sim.model.ocean.model.velocities
    η = sim.model.ocean.model.free_surface.displacement
    date = start_date + Second(round(Int, sim.model.clock.time))
    @printf("%s  |u|=%.2f |v|=%.2f  η∈[%+.2f,%+.2f]  (%.0f s)\n",
            Dates.format(date, "yyyy-mm-dd HH:MM"), maximum(abs, u), maximum(abs, v),
            minimum(filter(isfinite, interior(η))), maximum(filter(isfinite, interior(η))),
            time() - wall[]); wall[] = time()
    return nothing
end
add_callback!(simulation, progress, TimeInterval(6hours))

const λf = λnodes(grid, Face());   const λc = λnodes(grid, Center())
const φf = φnodes(grid, Face());   const φc = φnodes(grid, Center())
const zc = znodes(grid, Center())
function report_velocity_spike!(sim)
    u, v, w = sim.model.ocean.model.velocities
    ua, va = Array(interior(u)), Array(interior(v))

    iu = argmax(abs.(ua))
    if abs(ua[iu]) > SPIKE_THRESHOLD
        i, j, k = iu.I
        @printf("  SPIKE u=%+.2f m/s at (i=%d,j=%d,k=%d) λ=%.3f φ=%.3f z=%.1f m  t=%s\n",
                ua[iu], i, j, k, λf[i], φc[j], zc[k], prettytime(sim))
    end

    iv = argmax(abs.(va))
    if abs(va[iv]) > SPIKE_THRESHOLD
        i, j, k = iv.I
        @printf("  SPIKE v=%+.2f m/s at (i=%d,j=%d,k=%d) λ=%.3f φ=%.3f z=%.1f m  t=%s\n",
                va[iv], i, j, k, λc[i], φf[j], zc[k], prettytime(sim))
    end
    return nothing
end
add_callback!(simulation, report_velocity_spike!, TimeInterval(1hours))

oc = ocean.model
outfile = joinpath(@__DIR__, "..", TAG * ".jld2")
save_fields = (u = Field(oc.velocities.u; indices = (:, :, grid.Nz)),
               v = Field(oc.velocities.v; indices = (:, :, grid.Nz)),
               T = Field(oc.tracers.T;   indices = (:, :, grid.Nz)),
               S = Field(oc.tracers.S;   indices = (:, :, grid.Nz)))
# The surface salinity restoring flux (psu m/s, positive out of the ocean), saved with the surface fields
@inline sss_restoring_flux(i, j, k, grid, r, clock, S) = getbc(r, i, j, grid, clock, (; S))
if !isnothing(sss_restoring)
    sss_flux = Field(KernelFunctionOperation{Center, Center, Nothing}(sss_restoring_flux, grid, sss_restoring, oc.clock, oc.tracers.S))
    save_fields = merge(save_fields, (; sss_flux))
end
volume_fields = (u = oc.velocities.u, v = oc.velocities.v, w = oc.velocities.w,
                 T = oc.tracers.T, S = oc.tracers.S)
η_out = (; η = oc.free_surface.displacement)

# Checkpoint/restart. IMPORTANT for the FilteredTimeInterval writers below: nothing about the filter's
# internal running window is checkpointed (see NumericalEarth/CLAUDE.md's 2026-09-11 evening
# session), so a run picked up from a checkpoint starts those filters from scratch and needs
# `window/2` (2.5 days, for the 5-day daily window) of fresh integration before their output is
# valid again. To keep the daily/pentad output CONTINUOUS across a restart, don't pick up from a
# checkpoint taken at the exact split time — pick up from one taken `window/2` EARLIER than that,
# so the filter has rebuilt its window by the time simulated time reaches the split. Concretely:
# if segment 1 stops at t=S, checkpoint at t=S-2.5days is the one to resume segment 2 from (its
# valid daily output then starts at (S-2.5days)+2.5days = S, exactly where segment 1's own valid
# output ends). `MAB_PICKUP`: "false" (default, fresh start) | "true" (latest checkpoint) |
# an iteration number | a checkpoint filepath.
checkpoint_dir = isabspath(TAG) ? dirname(TAG) : joinpath(@__DIR__, "..")
pickup_raw = get(ENV, "MAB_PICKUP", "false")
PICKUP = pickup_raw == "false" ? false :
         pickup_raw == "true"  ? true  :
         occursin(r"^\d+$", pickup_raw) ? parse(Int, pickup_raw) :
         pickup_raw   # a checkpoint filepath

# A restart wipes and restarts every writer's file unless told not to — `!fresh_start` keeps
# a picked-up run's earlier segment instead of overwriting it.
fresh_start = PICKUP == false

# Raw (tidal) hourly output for the SSH animation, plus de-tided daily/pentad output
# (FilteredTimeInterval, then named LowPassFilter — see NumericalEarth/CLAUDE.md's 2026-09-11 evening session for the original
# 5-day-window/40h-cutoff daily + 10-day-window/10-day-cutoff pentad spec).
#
# Daily uses window=6days (not Oceananigans' own 5-day default) so its half-window (3 days) is
# an EXACT multiple of the 1-day interval: `FilteredTimeInterval`'s first-valid-frame time is
# `ceil(Int, (t + window/2) / interval) * interval` (see low_pass_filter.jl) — with the 5-day
# default that's `ceil(2.5) = 3`, a day later than the naive `window/2 = 2.5`; with 6 days it's
# `ceil(3) = 3` too (same first frame, same 2×3=6-day restart-continuity offset — see
# NumericalEarth/CLAUDE.md's mab_1month_continuous session), but now by exact arithmetic rather
# than a rounding coincidence. Cutoff (40h default) is unchanged, so what gets filtered out
# doesn't change, just a slightly wider Lanczos taper. Scoped to our own calls, not a change to
# Oceananigans' own default. Pentad's window=10days/interval=5days is already exact
# (half-window 5 = 1×interval) and doesn't need this.
simulation.output_writers[:surface] = JLD2Writer(oc, save_fields;
    filename = outfile, schedule = TimeInterval(3hours), overwrite_files = fresh_start)
simulation.output_writers[:surface_daily] = JLD2Writer(oc, save_fields;
    filename = joinpath(@__DIR__, "..", TAG * "_surface_daily.jld2"),
    schedule = FilteredTimeInterval(LanczosKernel(6days; cutoff = 40hours); interval = 1days), overwrite_files = fresh_start)
# Full-depth u/v/T/S, same de-tided daily cadence as surface_daily — for transects and other
# uses that need more than the top level (e.g. plot_temperature_transects.jl, which otherwise
# has to fall back to pulling a full 3D snapshot out of a checkpoint instead).
simulation.output_writers[:volume_daily] = JLD2Writer(oc, volume_fields;
    filename = joinpath(@__DIR__, "..", TAG * "_volume_daily.jld2"),
    schedule = FilteredTimeInterval(LanczosKernel(6days; cutoff = 40hours); interval = 1days), overwrite_files = fresh_start)
simulation.output_writers[:eta] = JLD2Writer(oc, η_out;
    filename = joinpath(@__DIR__, "..", TAG * "_eta.jld2"),
    schedule = TimeInterval(1hours), overwrite_files = fresh_start)
simulation.output_writers[:eta_daily] = JLD2Writer(oc, η_out;
    filename = joinpath(@__DIR__, "..", TAG * "_eta_daily.jld2"),
    schedule = FilteredTimeInterval(LanczosKernel(6days; cutoff = 40hours); interval = 1days), overwrite_files = fresh_start)
simulation.output_writers[:eta_pentad] = JLD2Writer(oc, η_out;
    filename = joinpath(@__DIR__, "..", TAG * "_eta_pentad.jld2"),
    schedule = FilteredTimeInterval(LanczosKernel(10days; cutoff = 10days); interval = 5days), overwrite_files = fresh_start)

simulation.output_writers[:checkpointer] = Checkpointer(model;
    schedule = TimeInterval(CHECKPOINT_EVERY), dir = checkpoint_dir,
    prefix = basename(TAG) * "_checkpoint", overwrite_files = false, cleanup = false)

@info "running: tangential=$TANGENTIAL  Uᵉˣᵗ=$UEXT_MODE  τ_in=$(TAU_IN/86400) day(s)  " *
      "tides=$(join(TIDE_CONSTITUENTS, ",")) velocity=$VELOCITY_SCHEME$(VELOCITY_SCHEME == "oblique" ? "(τ_in=$(OBLIQUE_TAU_IN/days)d, τ_out=$(OBLIQUE_TAU_OUT/days)d, w=$PHASE_SPEED_WEIGHT)" : "") consistent_ubc=$CONSISTENT_UBC tracers=$TRACER_SCHEME reservoir(L_in=$RESERVOIR_L_IN, L_out=$RESERVOIR_L_OUT)  " *
      "$(sim_days) days  pickup=$PICKUP"
run!(simulation; pickup = PICKUP, checkpoint_at_end = true)
println("\n✅ done — $(TAG)")
