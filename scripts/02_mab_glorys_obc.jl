# ==================================================================
# 02 — MAB with GLORYS-driven OPEN BOUNDARIES (replacing the edge restoring)
#
# Script 01 nudges the interior toward GLORYS in a 6-cell edge sponge. That
# works, but it is not an open boundary: the domain edge is still effectively
# closed and the reanalysis is pulled in by a body force. This script instead
# drives the boundaries themselves from GLORYS, and removes the restoring, so
# the boundary condition is the only thing communicating with the outside.
#
# WHY THIS IS THE RIGHT NEXT STEP. Re-running the bare-Oceananigans MAB
# (script 9a) with oblique radiation cut the OUTFLOW-face residual by a factor
# of 76 and largely removed the edge jet — but the inflow/outflow attribution
# then flipped from 0.18 to 4.97, i.e. essentially all the remaining error sat
# where water ENTERS. That configuration prescribes an exterior state of exactly
# zero, so the boundary's instruction on inflow is "arrive at rest" while the
# real MAB inflow is nothing of the kind. No radiation scheme can invent the
# Gulf Stream. The fix is data, and this is the data.
#
# THE FULL TRIAD, with real exterior values for the first time:
#
#   u, v NORMAL      NormalFlowBoundaryCondition(Interpolated(GLORYS); scheme)
#   u, v TANGENTIAL  ValueBoundaryCondition(Interpolated(GLORYS); scheme)   ← see note
#   T, S             ValueBoundaryCondition(Interpolated(GLORYS); scheme)
#   U, V barotropic  Flather, fed the DEPTH-INTEGRATED GLORYS transport and zos
#   η                Chapman (radiates only; Flather is the data entry point)
#
# The barotropic pair is not optional. The split-explicit corrector replaces the
# depth mean of the boundary velocity with whatever the barotropic solver says,
# so driving `u` at the boundary while leaving `U` unconditioned gets the inflow
# wiped — measured exactly, and it produced *precisely zero* response until `U`
# was given a condition too.
#
# NOTE ON THE TANGENTIAL COMPONENT. NumericalEarth's `parent_boundary_conditions`
# builds this whole structure from a parent FieldTimeSeries, but applies `schemes`
# only on `NormalFlow` sides — its docstring recommends a stiff Dirichlet for the
# tangential component. In nested tests that was the WORST of the three options
# (enstrophy 1.611 against a truth of 1.000, versus 0.330 for `Gradient(0)` and
# 0.923 radiated): imposing the tangential velocity puts a discontinuity straight
# into ∂v/∂x, and since ζ = ∂v/∂x − ∂u/∂y that appears as a vorticity band sitting
# on the boundary. So the boundary conditions here are hand-built from
# `Interpolated` in order to give the tangential component a scheme too.
#
# Run:  julia -t 8 --project=. scripts/02_mab_glorys_obc.jl
# ==================================================================

using NumericalEarth
using NumericalEarth.DataWrangling: Metadata, MetadataSet, BoundingBox
using NumericalEarth.NestedModels: Interpolated
using CopernicusMarine              # activates the GLORYS download backend
using CopernicusClimateDataStore    # activates the ERA5 download backend
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: ExponentialDiscretization, znodes, λnodes, φnodes
using Oceananigans.BoundaryConditions: PerturbationAdvection, NormalRadiation,
                                       GravityWaveRadiationBoundaryCondition,
                                       SurfaceWaveRadiationBoundaryCondition
using Dates, Printf, Statistics

include(joinpath(@__DIR__, "glorys_bathymetry.jl"))

# Copied from ../NumericalEarth/scripts/02_mab_glorys_obc.jl (keep in sync there — see that
# repo's CLAUDE.md for the session history). Only change from the original: DATA_DIR defaults
# to NumericalEarth's own cache (sibling repo) instead of a fresh, empty one here, so this
# doesn't re-download the ~380 MB of GLORYS/ETOPO/ERA5 data that's already sitting there.
# Override with MAB_DATA_DIR if that layout doesn't hold (e.g. NumericalEarth isn't a sibling).
const DATA_DIR   = get(ENV, "MAB_DATA_DIR", joinpath(@__DIR__, "..", "..", "NumericalEarth", "data"))
const resolution = 1 / 12
const Nz = 40

# `Interpolated` evaluates the parent at the child's boundary FACE, which lies half
# a cell outside the outermost parent cell CENTRE — so the parent source must
# BRACKET the child grid, not merely coincide with it. Building the GLORYS series
# on the model grid fails validation for exactly that reason.
#
# The cached GLORYS covers −76…−64, 34…42, so rather than re-download a larger box
# the MODEL domain is inset by `n_pad` cells and the GLORYS source keeps the full
# extent. Costs ~0.17° of domain on each side and nothing scientifically.
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
const sim_days   = parse(Int, get(ENV, "MAB_DAYS", "10"))
const stop_date  = start_date + Day(sim_days)

# variant switches
const TANGENTIAL = get(ENV, "MAB_TANGENTIAL", "radiated")   # "radiated" | "prescribed"
const TAU_IN     = parse(Float64, get(ENV, "MAB_TAU_IN", "1")) * days   # inflow nudging
const TAG        = get(ENV, "MAB_TAG", "mab_obc_glorys")
const MATCH_BATHY = get(ENV, "MAB_MATCH_BATHY", "true") == "true"   # blend model H → GLORYS deptho near the edges
const N_MATCH     = parse(Int, get(ENV, "MAB_N_MATCH", "4"))
const UEXT_MODE    = get(ENV, "MAB_UEXT", "native")   # Flather Uᵉˣᵗ: "native" (true depth mean from GLORYS's own
                                                   #             vertical grid × the model's wet depth — see below)
                                                   #             | "masked" (transport through the model's wet
                                                   #             column, but on GLORYS resampled to the MODEL's z)
                                                   #             | "legacy" (unmasked full-depth Integral, for A/B)
const SPIKE_THRESHOLD = parse(Float64, get(ENV, "MAB_SPIKE_THRESHOLD", "2.0"))   # m/s, see report_velocity_spike!

mkpath(DATA_DIR)

# ---------------- grid + bathymetry (identical to script 01) ----------------
z = ExponentialDiscretization(Nz, -4000, 0; scale = 1400)
grid = LatitudeLongitudeGrid(CPU(); size = (Nλ, Nφ, Nz),
                             longitude = λ_bounds, latitude = φ_bounds, z, halo = (7, 7, 7))
bottom_height = regrid_bathymetry(grid; dataset = ETOPO2022(), height_above_water = 1,
                                  minimum_depth = 10, major_basins = 1, interpolation_passes = 10)

# Match the model bathymetry to GLORYS `deptho` in a band near every open
# boundary. U = ∫u dz, so a depth mismatch there is a transport error: the
# unmatched western shelf differs from GLORYS by up to ~60%, which drives a
# spurious inflow jet (the SW-corner max|u| = 3.5 m/s artefact).
if MATCH_BATHY
    match_boundary_bathymetry!(bottom_height, grid, DATA_DIR;
                               region = BoundingBox(longitude = data_λ, latitude = data_φ),
                               n_match = N_MATCH, minimum_depth = 10)
end

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))
@show grid

# ---------------- GLORYS as FieldTimeSeries on the model grid ----------------
region = BoundingBox(longitude = data_λ, latitude = data_φ)

# The GLORYS source grid: full data extent, so it brackets the model on every side.
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

# ---------------- the barotropic exterior, from GLORYS ----------------
# Flather is Uᵇ = Uᵉˣᵗ ± √(gH)(ηᵇ − ηᵉˣᵗ) and is the ONLY data entry point for the
# barotropic mode (Chapman on η radiates but carries no data).
#
# Uᵉˣᵗ is the GLORYS transport through the MODEL's wet column at each boundary
# cell: Σₖ Δzₖ uₖ over the cells the GridFittedBottom leaves wet (cell centre
# above the floor), so that Uᵉˣᵗ/H_wet is the depth mean of the GLORYS profile
# over the column the model actually has.
#
# It is deliberately NOT `Integral(fts[n], dims = 3)`. `src_grid` is a plain
# LatitudeLongitudeGrid with no immersed bottom, so that integral runs over the
# whole ~4000 m column, including the inpainted velocity that fills the cells
# below GLORYS's real sea floor (drawn from offshore water). On the shelf that
# transport is then divided by a column tens of metres deep. Emulated from the
# day-1 GLORYS files, at the west face, model row 7 (the (1,7) spike): implied
# depth-mean Uᵉˣᵗ/H ≈ −3.4 m/s against a true depth mean of +0.2 m/s. Matching
# the boundary bathymetry cannot fix that — it only changes H, not the numerator.
# `report_uext` below prints the real numbers; `MAB_UEXT=legacy` restores the
# old behaviour for an A/B run.
@info "GLORYS transports for the Flather conditions (Uᵉˣᵗ mode = $UEXT_MODE)…"

const zc_src = collect(znodes(src_grid, Center()))     # source-grid layers = model layers
const Δz_src = diff(collect(znodes(src_grid, Face())))
const bh     = Array(interior(bottom_height))[:, :, 1]  # model floor height (<0 ocean), after boundary matching
const floors = ((bh[1, :], bh[end, :]),                 # (west, east) rows  — u boundaries
                (bh[:, 1], bh[:, end]))                 # (south, north) columns — v boundaries

# `col` is (source index along the boundary, Nz). Source index s ↔ model index
# s − n_pad (the model is inset by n_pad cells), clamped at the ends. Returns the
# wet-column transport and the wet column depth H_wet = Σₖ Δzₖ (0 where dry).
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

# boundary columns of a (Nλ_src, Nφ_src, Nz_expected) array: dim 1 → (west, east), 2 → (south, north).
# `Nz_expected` defaults to the model's `Nz`; the native-resolution path below passes GLORYS's
# own (larger) level count instead.
function boundary_columns(A, dim; Nz_expected = Nz)
    size(A) == (Nλ_src, Nφ_src, Nz_expected) ||
        error("unexpected GLORYS FieldTimeSeries frame size $(size(A)); expected $((Nλ_src, Nφ_src, Nz_expected))")
    p = n_pad + 1
    dim == 1 ? (A[p, :, :], A[end-n_pad, :, :]) : (A[:, p, :], A[:, end-n_pad, :])
end

function boundary_transport_slabs(fts, dim)
    times = collect(fts.times)
    slabs = map(eachindex(times)) do n
        if UEXT_MODE == "legacy"
            F = Field(Integral(fts[n], dims = 3))
            compute!(F)
            I = Array(interior(F))
            p = n_pad + 1
            dim == 1 ? (I[p, :, 1], I[end-n_pad, :, 1]) : (I[:, p, 1], I[:, end-n_pad, 1])
        else
            lo, hi = boundary_columns(Array(interior(fts[n])), dim)
            (wet_transport(lo, floors[dim][1])[1], wet_transport(hi, floors[dim][2])[1])
        end
    end
    return times, first.(slabs), last.(slabs)
end

# ---------------- "native" mode: true depth mean from GLORYS's OWN vertical grid ----------------
# `wet_transport` above already fixes the wrong DENOMINATOR (integrating past the model's
# floor). But its numerator still isn't quite raw GLORYS either: `fts_u`/`fts_v` live on
# `src_grid`, which shares the MODEL's stretched z, so by the time `wet_transport` sees the
# profile it has already been vertically resampled onto the model's handful of coarse shelf
# layers — a partial-cell error of up to half a model layer wherever the (boundary-matched)
# model floor doesn't land exactly on a model level.
#
# `native` decouples the two: ū is the true depth mean of `uo`/`vo` over GLORYS's OWN valid
# levels — its real bottom, from the static `deptho` field, not any model floor — computed on
# GLORYS's own (finer) vertical grid. Only H_wet — the model's own wet depth, already computed
# above by `wet_transport` — is model-derived, because that's the H the model's Flather/
# barotropic mode actually needs. `z_interfaces` (GLORYS.jl) approximates each native level's
# thickness from its midpoint spacing; we don't have GLORYS's exact partial-cell thickness
# field, so this, like the model-z version, is still an approximation — just a much finer one.
#
# `FieldTimeSeries(metadata, CPU())` (no target grid) would build GLORYS's own auto-snapped
# `native_grid`, but that pads/snaps its horizontal box independently of ours (measured: 146×98
# vs our 144×96 — NOT a clean n_pad offset, since the snap depends on where our box edges fall
# inside a native cell). Horizontally we want EXACTLY `src_grid`'s lattice, so the existing
# n_pad-offset indexing stays valid; only the vertical grid should be GLORYS's own. So: build a
# grid identical to `src_grid` except for `z`, using GLORYS's real depth interfaces there — this
# routes through the same horizontal-interpolation-onto-src_grid path already proven correct by
# `fts_u`/`fts_v`, just without a vertical resample.
z_native = NumericalEarth.DataWrangling.z_interfaces(meta(:u_velocity))
Nz_native = length(z_native) - 1
src_grid_native_z = LatitudeLongitudeGrid(CPU(); size = (Nλ_src, Nφ_src, Nz_native),
                                          longitude = data_λ, latitude = data_φ,
                                          z = z_native, halo = (7, 7, 7))
@info "building the native-resolution GLORYS u/v for the true depth mean… ($Nz_native levels vs the model's $Nz)"
fts_u_native = FieldTimeSeries(meta(:u_velocity), src_grid_native_z)
fts_v_native = FieldTimeSeries(meta(:v_velocity), src_grid_native_z)
Hg_src = glorys_deptho_on_grid(src_grid, DATA_DIR; region)   # GLORYS's own depth (m, positive down), NaN on land
true_floors = ((-Hg_src[1, :], -Hg_src[end, :]),             # (west, east) rows  — u boundaries
               (-Hg_src[:, 1], -Hg_src[:, end]))             # (south, north) columns — v boundaries

# Depth mean of a native-resolution boundary column over GLORYS's OWN valid levels. `b` is
# `NaN` where GLORYS itself has no ocean; `zc[k] > NaN` is `false` for every `k`, so those
# columns fall out naturally (H stays 0 ⇒ ubar 0) without a separate land check.
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
        # H_wet doesn't depend on GLORYS's values, only on the model's own floor — reuse
        # `wet_transport`'s H_wet return with a throwaway (model-shaped) column of zeros.
        Hwet_lo = last(wet_transport(zeros(size(lo, 1), Nz), floors[dim][1]))
        Hwet_hi = last(wet_transport(zeros(size(hi, 1), Nz), floors[dim][2]))
        (ū_lo .* Hwet_lo, ū_hi .* Hwet_hi)
    end
    return times, first.(slabs), last.(slabs)
end

# Diagnostic (frame 1): implied depth-mean boundary velocity Uᵉˣᵗ/H_wet from all three modes,
# per side. legacy tracks the spike; masked and native should both track GLORYS, with native
# the more precise of the two (see the comment block above).
function report_uext(fts, fts_native, dim, names)
    A = Array(interior(fts[1]))
    F = Field(Integral(fts[1], dims = 3)); compute!(F)
    I = Array(interior(F)); p = n_pad + 1
    cols = boundary_columns(A, dim)
    legs = dim == 1 ? (I[p, :, 1], I[end-n_pad, :, 1]) : (I[:, p, 1], I[:, end-n_pad, 1])

    zc_n = collect(znodes(fts_native.grid, Center()))
    Δz_n = diff(collect(znodes(fts_native.grid, Face())))
    cols_n = boundary_columns(Array(interior(fts_native[1])), dim; Nz_expected = length(zc_n))

    for side in 1:2
        T, Hwet = wet_transport(cols[side], floors[dim][side])
        wet = findall(Hwet .> 0)
        isempty(wet) && continue
        ū_leg = legs[side][wet] ./ Hwet[wet]
        ū_msk = T[wet]          ./ Hwet[wet]
        ū_nat = true_depth_mean(cols_n[side], zc_n, Δz_n, true_floors[dim][side])[wet]
        iL, iM, iN = argmax(abs.(ū_leg)), argmax(abs.(ū_msk)), argmax(abs.(ū_nat))
        @printf("  Uᵉˣᵗ/H_wet, %-5s frame 1: legacy max|ū| = %6.3f m/s at %3d   masked max|ū| = %6.3f m/s at %3d   native max|ū| = %6.3f m/s at %3d\n",
                names[side], abs(ū_leg[iL]), wet[iL] - n_pad, abs(ū_msk[iM]), wet[iM] - n_pad, abs(ū_nat[iN]), wet[iN] - n_pad)
    end
end
report_uext(fts_u, fts_u_native, 1, ("west", "east"))
report_uext(fts_v, fts_v_native, 2, ("south", "north"))

if UEXT_MODE == "native"
    tU, U_lo, U_hi = boundary_transport_slabs_native(fts_u_native, true_floors[1], 1)   # west, east
    tV, V_lo, V_hi = boundary_transport_slabs_native(fts_v_native, true_floors[2], 2)   # south, north
else
    tU, U_lo, U_hi = boundary_transport_slabs(fts_u, 1)   # west, east
    tV, V_lo, V_hi = boundary_transport_slabs(fts_v, 2)   # south, north
end

η_slabs = map(eachindex(fts_η.times)) do n
    A = Array(interior(fts_η[n]))[:, :, 1]
    p = n_pad + 1
    (A[p, :], A[end-n_pad, :], A[:, p], A[:, end-n_pad])   # west, east, south, north
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

U_west(j, k, grid, clock, f)  = (lerp_slab(tU, U_lo, j, clock.time), lerp_slab(tU, η_w, j, clock.time))
U_east(j, k, grid, clock, f)  = (lerp_slab(tU, U_hi, j, clock.time), lerp_slab(tU, η_e, j, clock.time))
V_south(i, k, grid, clock, f) = (lerp_slab(tV, V_lo, i, clock.time), lerp_slab(tV, η_s, i, clock.time))
V_north(i, k, grid, clock, f) = (lerp_slab(tV, V_hi, i, clock.time), lerp_slab(tV, η_n, i, clock.time))

# ---------------- the boundary conditions ----------------
normal_scheme    = PerturbationAdvection(inflow_timescale = TAU_IN, outflow_timescale = Inf)
tangential_sch   = NormalRadiation(inflow_timescale = TAU_IN, outflow_timescale = Inf)
tracer_scheme    = NormalRadiation(inflow_timescale = TAU_IN, outflow_timescale = Inf)

tangential(fts) = TANGENTIAL == "prescribed" ?
    ValueBoundaryCondition(Interpolated(fts)) :
    ValueBoundaryCondition(Interpolated(fts); scheme = tangential_sch)

u_bcs = FieldBoundaryConditions(
    west  = NormalFlowBoundaryCondition(Interpolated(fts_u); scheme = normal_scheme),
    east  = NormalFlowBoundaryCondition(Interpolated(fts_u); scheme = normal_scheme),
    south = tangential(fts_u),
    north = tangential(fts_u))

v_bcs = FieldBoundaryConditions(
    south = NormalFlowBoundaryCondition(Interpolated(fts_v); scheme = normal_scheme),
    north = NormalFlowBoundaryCondition(Interpolated(fts_v); scheme = normal_scheme),
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

# ---------------- ocean simulation: NO restoring forcing ----------------
ocean = ocean_simulation(grid; boundary_conditions)

set!(ocean.model, MetadataSet(:temperature, :salinity, :u_velocity, :v_velocity;
                              dataset = glorys, date = start_date, dir = DATA_DIR, region))
set!((; free_surface = ocean.model.free_surface.displacement),
     MetadataSet(:free_surface; dataset = glorys, date = start_date, dir = DATA_DIR, region))
let η = ocean.model.free_surface.displacement
    η .-= mean(filter(isfinite, interior(η)))
end

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

# Ported from 03_mab_m2_tide_smoke.jl's report_velocity_spike! (commit bc6eac4): the
# 6-hourly progress line only reports the domain max, not where it is, which isn't
# enough to tell a real boundary artefact from a healthy Gulf-Stream max. Logs the
# location of any hourly |u|/|v| above MAB_SPIKE_THRESHOLD (default 2.0 m/s).
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
# Output the fields only. The GLORYS-backed BCs carry `Metadata` objects whose
# types do not exist outside the run session, so `FieldTimeSeries` read-back of a
# file that serialized them emits a wall of reconstruction warnings. Passing the
# raw interior arrays via `FunctionField`-free `Field`s and letting JLD2Writer
# serialize just the data avoids it.
save_fields = (u = Field(oc.velocities.u; indices = (:, :, grid.Nz)),
               v = Field(oc.velocities.v; indices = (:, :, grid.Nz)),
               T = Field(oc.tracers.T;   indices = (:, :, grid.Nz)),
               S = Field(oc.tracers.S;   indices = (:, :, grid.Nz)))
simulation.output_writers[:surface] = JLD2Writer(oc, save_fields;
    filename = outfile, schedule = TimeInterval(3hours), overwrite_files = true)
simulation.output_writers[:eta] = JLD2Writer(oc, (; η = oc.free_surface.displacement);
    filename = joinpath(@__DIR__, "..", TAG * "_eta.jld2"),
    schedule = TimeInterval(3hours), overwrite_files = true)

@info "running: tangential=$TANGENTIAL  τ_in=$(TAU_IN/86400) day(s)  $(sim_days) days"
run!(simulation)
println("\n✅ done — $(TAG)")
