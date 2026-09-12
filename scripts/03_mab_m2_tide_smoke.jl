# ==================================================================
# 03 — M2 TIDE SMOKE TEST for the Mid-Atlantic Bight
#
# The first run that puts both halves of the tidal forcing together:
#
#   BODY FORCE      the M2 equilibrium tide, +g∇η_eq  (TidalBodyForce)
#   BOUNDARY FORCE  TPXO10 M2 transport and sea level, through Flather
#                   on U/V with Chapman on η  (tidal_boundary_conditions)
#
# both built from ONE `TidalHarmonics`, so they share (ω, f, Θ) by
# construction and cannot drift out of phase.
#
# WHY M2 ALONE, AND WHY THAT IS ENOUGH TO VERIFY. The Rayleigh criterion
# says two constituents are separable only over a record longer than
# 2π/|ωᵢ−ωⱼ|. For M2 alone a few days suffice — a 6-day record fits it at
# condition number 1.4 — whereas the full semidiurnal group needs 183
# days to separate S2 from K2 (measured in verify_tidal_harmonics.jl).
# So a short single-constituent run can be checked QUANTITATIVELY, which
# a short ten-constituent run cannot. Ten constituents come later, with a
# year-long run.
#
# WHAT IS DELIBERATELY ABSENT. No stratification, no tracers, no GLORYS,
# no atmosphere. The barotropic M2 tide needs none of them, and their
# absence is what makes this a test of the tidal forcing rather than of
# everything at once. In particular it makes the comparison against TPXO
# apples-to-apples: TPXO is itself a barotropic tidal solution.
#
# With no stratification the 3-D velocity is essentially the barotropic
# velocity, so `u`/`v` get radiation conditions toward a ZERO exterior and
# the tide enters entirely through the barotropic pair. That is correct
# here rather than lazy: the split-explicit corrector overwrites the depth
# mean of the boundary velocity with the barotropic solver's value, so the
# tide arrives via Flather, and the baroclinic deviation radiating to zero
# is right when there is no baroclinic tide to carry.
#
# THE TEST. Forced at its boundaries by TPXO, the model must REPRODUCE
# TPXO in the interior. Harmonically analyse the model's free surface at
# every wet cell, compare amplitude and phase against the atlas, and
# report the complex RMS difference — the standard tidal-model skill
# metric. Secondarily, compare against NOAA gauge constants where the
# model has a wet cell.
#
# Env knobs:
#   MAB_DAYS=10        run length
#   MAB_NZ=24          vertical levels (see the vertical-grid note below)
#   MAB_BODY_FORCE=1   set to 0 to isolate the boundary forcing
#   MAB_CLOSURE=catke  vertical closure: catke (default) | nuz (uniform ν_z probe) | kepsilon
#   MAB_STRAT=none     none | deep | shelf  (horizontally uniform GLORYS profile)
#   MAB_TAG=...        output name
#
# Run:  julia -t 8 --project=. scripts/03_mab_m2_tide_smoke.jl
# ==================================================================

using NumericalEarth
using NumericalEarth.DataWrangling: BoundingBox
using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization, CATKEVerticalDiffusivity,
                                       TKEDissipationVerticalDiffusivity
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEMixingLength
using Oceananigans.Grids: ExponentialDiscretization, static_column_depthᶜᶜᵃ
using Oceananigans.BoundaryConditions: NormalRadiation, ObliqueRadiation,
                                       GravityWaveRadiationBoundaryCondition,
                                       SurfaceWaveRadiationBoundaryCondition
using Dates, Printf, Statistics

# The tidal FORCING now comes from the NumericalEarth package (branch `tides`). These scripts stay
# for the SKILL ANALYSIS — `reconstruct`, `harmonic_analysis` and the atlas readers, which the package
# deliberately does not provide — inside a module so their `TidalHarmonics` does not collide with it.
module LocalTides
using Oceananigans, Oceananigans.Units, NumericalEarth
using Dates, Printf, Statistics, LinearAlgebra
include(joinpath(@__DIR__, "tidal_harmonics.jl"))
include(joinpath(@__DIR__, "tpxo.jl"))
include(joinpath(@__DIR__, "tpxo_boundaries.jl"))
end

using .LocalTides: ALL_CONSTITUENTS, TPXO_DIR, load_tpxo, tpxo_elevation, tpxo_depth,
                   amplitude_and_phase, reconstruct, reconstruction_parameters,
                   harmonic_analysis, constituent_period

include(joinpath(@__DIR__, "noaa_harcon_mab.jl"))
include(joinpath(@__DIR__, "variable_bottom_drag.jl"))
include(joinpath(@__DIR__, "kepsilon_tuple_closure_patch.jl"))
include(joinpath(@__DIR__, "glorys_profiles.jl"))

const OUT_DIR = get(ENV, "MAB_OUT", joinpath(homedir(), "Data", "mab_tides"))
mkpath(OUT_DIR)

# ---------------- configuration ----------------
# Grid identical to scripts 01/02 so that verified numbers transfer directly
# when the tide is folded into the realistic configuration.
const resolution = 1 / 12
const Nz         = parse(Int, get(ENV, "MAB_NZ", "24"))
const n_pad      = 2
const data_λ     = (-76.0, -64.0)
const data_φ     = ( 34.0,  42.0)
const λ_bounds   = (data_λ[1] + n_pad*resolution, data_λ[2] - n_pad*resolution)
const φ_bounds   = (data_φ[1] + n_pad*resolution, data_φ[2] - n_pad*resolution)
const Nλ = round(Int, (λ_bounds[2] - λ_bounds[1]) / resolution)
const Nφ = round(Int, (φ_bounds[2] - φ_bounds[1]) / resolution)

const start_date = DateTime(2019, 4, 1)
const sim_days   = parse(Int, get(ENV, "MAB_DAYS", "10"))
const BODY_FORCE = get(ENV, "MAB_BODY_FORCE", "1") == "1"
const TAG        = get(ENV, "MAB_TAG", "mab_m2_" * get(ENV, "MAB_DRAG", "const") * "_" * get(ENV, "MAB_CLOSURE", "catke") *
                                      (get(ENV, "MAB_STRAT", "none") == "none" ? "" : "_strat" * get(ENV, "MAB_STRAT", "none")) *
                                      (get(ENV, "MAB_CONSTITUENTS", "M2") == "all" ? "_all10" : ""))

# "M2" (default) — the single-constituent test, validated by harmonic constants.
# "all" — all ten constituents, validated in the TIME DOMAIN against the full TPXO
#         prediction: a short record cannot separate them (S2–K2 and K1–P1 need
#         183 days, M2–N2 28; see verify_tidal_harmonics.jl), but a direct η(t)
#         comparison needs no separation at all, because the forcing is known.
const MULTI = get(ENV, "MAB_CONSTITUENTS", "M2") == "all"
const CONSTITUENTS = MULTI ? ALL_CONSTITUENTS : (:M2,)
# Quadratic bottom drag. The stress is −Cᴰ·u·|u| (see
# NumericalEarth.Oceans.u_quadratic_bottom_drag), so Cᴰ is the DIMENSIONLESS drag
# coefficient, not a rate, and the flux has the m² s⁻² a kinematic flux BC wants.
#
# Is 0.003 right? Cᴰ is calibrated against the velocity at a log-layer reference
# height, so the answer depends on which cell it multiplies — and on THIS vertical
# grid the shelf is well resolved rather than a single column: a 50 m column has
# 10 wet levels with the lowest 11.6 m thick, so z_ref ≈ 5.8 m. Then
#
#     Cᴰ = (κ / ln(z_ref/z₀))²  =  2.2e-3 (z₀ = 1 mm)
#                                  2.9e-3 (z₀ = 3 mm)
#                                  4.2e-3 (z₀ = 10 mm)
#
# so 3e-3 sits in the middle of the plausible range for shelf roughness, and also
# matches ocean_simulation's default. Known approximation: strictly Cᴰ should then
# vary with the local bottom-cell thickness — in 20 m of water the lowest cell is
# 1.07 m and the log law wants ≈6.5e-3 — but a single value is the standard
# compromise and the shelf-wide error it implies is well inside the uncertainty in
# z₀ itself.
const Cᴰ       = parse(Float64, get(ENV, "MAB_CD", "0.003"))
# "const"    — one Cᴰ everywhere (the classical choice)
# "variable" — depth-dependent log-law Cᴰ, applied SEMI-IMPLICITLY
# "variable_explicit" — same field, explicit; unstable in 384 inner-shelf cells
const DRAG = get(ENV, "MAB_DRAG", "const")
const z₀   = parse(Float64, get(ENV, "MAB_Z0", "0.003"))   # bottom roughness [m]
# Multiplier on the log-law Cᴰ field. Only meaningful with DRAG = "variable",
# because the semi-implicit treatment has no stability ceiling — which is what
# makes it possible to ask whether dissipation can fix the shelf AT ALL, by
# pushing Cᴰ far past anything the explicit scheme could run.
const CD_SCALE = parse(Float64, get(ENV, "MAB_CD_SCALE", "1.0"))
# Vertical viscosity. This is NOT a detail: bottom drag removes momentum from the
# lowest cell only, and ν sets how far that stress reaches up the column. The
# oscillatory boundary-layer thickness is δ = √(2ν/ω), which at ν = 1e-4 is 1.2 m —
# so in 33 m of shelf water the drag damps one thin cell, |u| there goes to zero,
# the quadratic drag goes to zero with it, and the rest of the column never feels
# it. A real tidal bottom boundary layer has ν ~ κ u★ z ≈ 1e-2–1e-1 m²/s.
#
# MEASURED, because this turned out to be the single most important number in the
# configuration. Bottom drag removes momentum from the lowest cell only, and ν sets
# how far that stress reaches up the column; the oscillatory boundary layer is
# δ = √(2ν/ω). At ν = 1e-4 (δ = 1.2 m) the drag damps one thin cell in 33 m of
# shelf water, |u| there goes to zero, the quadratic drag goes to zero with it, and
# the column above never feels it — bottom drag is effectively INERT. That is why
# Cᴰ sweeps from 0 to 20× the log-law value changed nothing, and why more drag
# appeared to make things worse.
#
#   ν_z      δ       domain cRMS   shelf<50 cRMS/A   |η|max d10   amphidrome error
#   1e-4    1.2 m      0.1326           0.90           2.84 m         79 km
#   1e-2     12 m      0.1075           0.73           1.48 m
#   1e-1     38 m      0.0888           0.60           1.23 m         12 km
#
# THIS IS A PROBE, NOT A PHYSICAL SETTING. A real tidal bottom boundary layer has
# ν ≈ κ u★ z — zero at the bed, growing upward, capped by the layer thickness, and
# varying through the tidal cycle. A uniform 0.1 m²/s gives the abyss a 38 m
# boundary layer where the tidal velocity is ~1 cm/s. It works here only because
# the value is roughly right ON THE SHELF, which is where it matters. The correct
# fix is a closure that computes it — which is exactly what CATKE does and what
# `ocean_simulation` would have supplied. Stripping the model down to a barotropic
# configuration removed the closure, and with it the thing that makes bottom drag
# function at all. That is the lesson for the realistic script, not this number.
const ν_z = parse(Float64, get(ENV, "MAB_NUZ", "1e-1"))

# Vertical closure, the fix the ν_z note above asks for:
#   "nuz"      — the uniform ν_z probe (kept for comparison; reproduces the earlier results)
#   "catke"    — CATKEVerticalDiffusivity (DEFAULT), which COMPUTES the bottom boundary layer:
#                its mixing length is limited by height above the bottom (Cᵇ) and an
#                explicit bottom-cell TKE sink (Cᵂϵ) damps TKE at solid bottoms, so the
#                drag-induced shear builds a law-of-the-wall viscosity instead of a uniform guess.
#   "kepsilon" — TKEDissipationVerticalDiffusivity, the two-equation k-ϵ closure
#                (Burchard & Bolding 2001; Umlauf & Burchard 2003, 2005). It has NEITHER of
#                CATKE's bottom terms above — no bottom-distance mixing-length limiter and no
#                bottom-cell TKE sink, only a wind/wave-forced flux at the TOP — so any BBL
#                behavior here comes entirely from the dynamic e-ϵ shear-production coupling.
# CATKE and k-ϵ both require a buoyancy model — with buoyancy = nothing they have no
# top_buoyancy_flux method — so both carry SeawaterBuoyancy(TEOS10) with UNIFORM T and S.
# N² = g(α∂zT − β∂zS) is then exactly zero: this step isolates the closure from
# stratification, and the next step changes only the initial T/S profile.
#
# WHY CATKE IS THE DEFAULT (measured 2026-09-11; 10-day M2 runs, constant Cᴰ):
#
#   vertical mixing       stratification   domain cRMS   shelf <50 m cRMS/A
#   uniform ν_z = 0.1         none            0.0860            0.59
#   CATKE                     none            0.1070            0.73
#   CATKE                     deep/slope      0.1054            0.72
#   CATKE                     shelf           0.1068            0.73
#
# The uniform 0.1 probe scores best, but it is not a candidate: a constant vertical
# viscosity that large would also erode the SURFACE stratification, and it gives the
# abyss a 38 m boundary layer. CATKE computes the mixing, and scores like
# ν_z ≈ 1e-2. Stratification in either regime moves shelf skill by ≤ 1%.
#
# OPEN — the bottom boundary layer. When shelf stratification cut CATKE's shelf
# viscosity ~5× (end-of-run column max 1.2e-2 → 2.5e-3 m²/s) the skill did not move,
# so the probe's advantage may be blunt damping of all vertical shear rather than
# BBL physics; unresolved. In very shallow water the surface and bottom boundary
# layers overlap, which is where this gets genuinely complicated. To be revisited in
# the more realistic, higher-resolution configurations.
const CLOSURE = get(ENV, "MAB_CLOSURE", "catke")
const T_uniform = 15.0    # °C
const S_uniform = 35.0    # g/kg

# Stratification (requires MAB_CLOSURE=catke, which is what carries buoyancy):
#   "none"  — uniform T/S, N² = 0 (step 1: the closure alone)
#   "deep"  — horizontally uniform mean GLORYS profile of columns deeper than
#             1000 m: realistic at the shelf break and slope, where barotropic-to-
#             internal tide conversion happens
#   "shelf" — the same for columns shallower than 200 m: realistic on the shelf,
#             where the error is; held constant below ~190 m
# Horizontally uniform by design — see glorys_profiles.jl for why there is no
# single representative April profile, and why a horizontal gradient is ruled out.
const STRAT = get(ENV, "MAB_STRAT", "none")
STRAT in ("none", "deep", "shelf") || error("MAB_STRAT must be none, deep or shelf, got $STRAT")
STRAT == "none" || get(ENV, "MAB_CLOSURE", "catke") == "catke" ||
    error("MAB_STRAT=$STRAT needs MAB_CLOSURE=catke (the only path that carries buoyancy)")
const ν_h      = 100.0          # m²/s lateral viscosity; damping time L²/ν ≈ 9 days,
                                # long against the 12.42 h period so the tide is not damped
const ramp     = 1days          # shared by the body force and the boundary tide
const n_discard = 4             # periods of ramp + adjustment to discard

@info "M2 smoke test: $(Nλ)×$(Nφ)×$(Nz), $(sim_days) days, body force = $(BODY_FORCE)"

# ---------------- grid + bathymetry ----------------
# VERTICAL GRID, chosen by measurement rather than inherited from script 02.
# Two things had to change, and both are about representing the water DEPTH,
# because for a barotropic tide the wave speed is √(gH) and the transport Flather
# imposes scales with H — depth error is the dominant error source.
#
#  (a) The BOTTOM: script 01/02 use -4000 m, but this box reaches 5365 m, so every
#      deeper cell was truncated. Measured depth error in the deep ocean: 12.9%,
#      i.e. ~6% in √(gH). Extending to -5500 m brings it to 1.65%.
#  (b) The TOP CELL: with Nz=12 and scale=1400 it is 65.5 m thick, which on a
#      shelf of median depth 50 m is hopeless. It also DESTABILISED the run: a
#      65.5 m cell next to a 13 m partial cell is an ill-conditioned thickness
#      ratio, and the blow-up was completely insensitive to Δt (tested at 120, 30
#      and 10 s — all blew up), confirming a thin-cell pathology rather than CFL.
#      Nz=24 with scale=700 gives a 0.82 m top cell, which both stabilises the run
#      and cuts the shelf depth error from 6.6% to 2.2%.
#
# A cautionary measurement worth keeping: raising minimum_fractional_cell_height
# to 0.5 also made the Nz=12 grid stable, and was nearly adopted on that basis —
# but it floors shallow water at ~33 m and gives a shelf depth error of 49%, WORSE
# than GridFittedBottom. A stable configuration is not automatically a better one.
z = ExponentialDiscretization(Nz, -5500, 0; scale = 700)
underlying_grid = LatitudeLongitudeGrid(CPU(); size = (Nλ, Nφ, Nz),
                                        longitude = λ_bounds, latitude = φ_bounds, z,
                                        halo = (7, 7, 7))
bottom_height = regrid_bathymetry(underlying_grid; dataset = ETOPO2022(),
                                  height_above_water = 1, minimum_depth = 10,
                                  major_basins = 1, interpolation_passes = 10)

# PartialCellBottom, NOT GridFittedBottom. GridFittedBottom snaps the sea floor to
# the nearest z FACE, which on this shelf turns 658 of 2377 shelf cells into DRY
# LAND and leaves a 39% rms depth error on the survivors. PartialCellBottom lets
# the bottom cell take a partial thickness matching the real depth: zero shelf
# cells lost, 2.2% rms depth error. Script 9a used PartialCellBottom for the same
# reason; script 02 uses GridFittedBottom because its signal is the baroclinic
# circulation rather than the tide.
grid = ImmersedBoundaryGrid(underlying_grid,
                            PartialCellBottom(bottom_height;
                                              minimum_fractional_cell_height = 0.2))

# SHALLOW WATER IS NOT THE PROBLEM, AND MUST NOT BE MASKED AWAY. Shallow columns
# run fine given enough vertical resolution, and on this grid they have it: the
# 0.82 m top cell means 5 m of water still gets 4 wet levels, 10 m gets 6, 20 m
# gets 8, and the PartialCellBottom floor is 0.165 m rather than the 13 m it was on
# the old grid. The Nz=12 failure was a RESOLUTION failure wearing a shallow-water
# disguise: a 65.5 m top cell cannot represent a 14 m column at all.
#
# This was tested rather than assumed. Masking every open-boundary cell shallower
# than 30 m to land changed the peak spurious velocity from 4.42 to 4.39 m/s, and a
# 50 m threshold only reached 2.60 — so shallow depth was not the cause, and
# masking would have discarded the MAB shelf, which is the physically interesting
# part of this domain. The real cause was a boundary-condition conflict (see the
# inflow_timescale note below).
#
# The binding limit on shallowness here is `minimum_depth = 10` in
# regrid_bathymetry above — a judgement about what a ~9 km cell can meaningfully
# represent, not a numerical constraint. The vertical grid would support less.

# ---------------- the shared astronomy ----------------
# The forcing reads the package astronomy; the analysis below keeps the script version, unchanged,
# so the skill numbers stay comparable with the earlier runs. The two agree to 0.02° in phase.
harmonics = TidalHarmonics(start_date; constituents = CONSTITUENTS, ramp_time = ramp)
analysis_harmonics = LocalTides.TidalHarmonics(start_date; constituents = CONSTITUENTS)
println(harmonics)

# ---------------- TPXO, and the boundary conditions ----------------
@info "sampling the atlas at the model open boundaries:"
tide_bcs = tidal_boundary_conditions(grid, harmonics, TPXO10Atlas(); dir = TPXO_DIR)

@info "loading TPXO10 over the MAB box for the skill analysis…"
window = load_tpxo(CONSTITUENTS; λ_bounds = (data_λ[1] - 1, data_λ[2] + 1),
                                 φ_bounds = (data_φ[1] - 1, data_φ[2] + 1))

# ---------------- bottom drag ----------------
if DRAG == "const"
    Cᴰ_field = nothing
    # NumericalEarth 0.8 passes these drag functions (μ, ub); 0.7 passes Cᴰ alone. ub = 0 is the same drag.
    drag_parameters = pkgversion(NumericalEarth) >= v"0.8" ? (μ = Cᴰ, ub = 0.0) : Cᴰ
    drag = (u_bottom = FluxBoundaryCondition(NumericalEarth.Oceans.u_quadratic_bottom_drag;
                                             discrete_form = true, parameters = drag_parameters),
            u_immersed = ImmersedBoundaryCondition(bottom =
                FluxBoundaryCondition(NumericalEarth.Oceans.u_immersed_bottom_drag;
                                      discrete_form = true, parameters = drag_parameters)),
            v_bottom = FluxBoundaryCondition(NumericalEarth.Oceans.v_quadratic_bottom_drag;
                                             discrete_form = true, parameters = drag_parameters),
            v_immersed = ImmersedBoundaryCondition(bottom =
                FluxBoundaryCondition(NumericalEarth.Oceans.v_immersed_bottom_drag;
                                      discrete_form = true, parameters = drag_parameters)))
    @info "bottom drag: constant Cᴰ = $Cᴰ"
else
    Cᴰ_field = log_law_drag_coefficient(grid; z₀, Cᴰ_bounds = (1e-3, 1e-2 * max(1, CD_SCALE)))
    CD_SCALE == 1 || (interior(Cᴰ_field, :, :, 1) .*= CD_SCALE)
    @info "bottom drag: depth-dependent log-law Cᴰ (z₀ = $z₀ m, ×$CD_SCALE), $(DRAG == "variable" ? "SEMI-IMPLICIT" : "explicit")"
    report_drag_field(Cᴰ_field, grid; Δt = 240.0, u_scale = 1.0)
    drag = variable_drag_boundary_conditions(Cᴰ_field; implicit = (DRAG == "variable"))
end

# u, v radiate toward a zero exterior: the tide arrives barotropically (see header).
#
# inflow_timescale = Inf, NOT 0, and this was measured the hard way. With
# inflow_timescale = 0 the scheme CLAMPS the boundary value to the exterior on
# inflow, which fights the split-explicit corrector: the clamp resets u while the
# corrector restores the barotropic transport, and the result is a boundary-trapped
# jet that grew to |v| = 4.4 m/s. It is the same failure script 9a documented --
# "the boundary's instruction on inflow is arrive at rest" -- and it is the CLAMP
# that does it, not the zero: prescribing the true TPXO tidal velocity as the
# exterior and keeping tau_in = 0 gave an identical 4.4 m/s.
#
# Measured alternatives, 18 h each, peak |v|:
#     zero exterior,  tau_in = 0      4.42 m/s   (jet)
#     zero exterior,  tau_in = Inf    0.94 m/s
#     tidal exterior, tau_in = 0      4.40 m/s   (jet)
#     tidal exterior, tau_in = 1 day  0.94 m/s
#
# The tidal exterior makes NO measurable difference here -- 0.9393 m/s either way,
# and |η|max is 1.0238 m in all four -- which is the empirical confirmation of the
# header's claim that the tide arrives barotropically and this condition only sets
# the baroclinic deviation. With no stratification there is none, so the zero
# exterior is not a shortcut but the honest choice. Once stratification is added
# (a baroclinic tide to radiate) the tidal exterior will start to matter, and
# tpxo_boundaries.jl already has what it needs to supply it.
open_scheme = NormalRadiation(inflow_timescale = Inf, outflow_timescale = Inf)

# The TANGENTIAL component at each boundary (u at south/north, v at west/east) has used
# NormalRadiation as a stand-in -- a scheme built for the component crossing the boundary,
# applied here to the one running along it. ObliqueRadiation is the scheme actually built for
# that role (Raymond & Kuo 1984): it radiates using both the normal and along-boundary
# derivatives, which a corner needs. Same inflow/outflow timescales as above, for the same
# reason: Inf on both, not the clamp-on-inflow default, so it does not fight the split-explicit
# corrector.
# "normal"  -- NormalRadiation on the tangential component too (the existing baseline)
# "oblique" -- ObliqueRadiation on the tangential component
const TANGENTIAL = get(ENV, "MAB_TANGENTIAL", "normal")
tangential_scheme = if TANGENTIAL == "normal"
    open_scheme
elseif TANGENTIAL == "oblique"
    ObliqueRadiation(inflow_timescale = Inf, outflow_timescale = Inf)
else
    error("MAB_TANGENTIAL must be \"normal\" or \"oblique\", got \"$TANGENTIAL\"")
end

u_bcs = FieldBoundaryConditions(
    west  = NormalFlowBoundaryCondition(0; scheme = open_scheme),
    east  = NormalFlowBoundaryCondition(0; scheme = open_scheme),
    south = ValueBoundaryCondition(0; scheme = tangential_scheme),
    north = ValueBoundaryCondition(0; scheme = tangential_scheme),
    bottom = drag.u_bottom,
    immersed = drag.u_immersed)
v_bcs = FieldBoundaryConditions(
    south = NormalFlowBoundaryCondition(0; scheme = open_scheme),
    north = NormalFlowBoundaryCondition(0; scheme = open_scheme),
    west  = ValueBoundaryCondition(0; scheme = tangential_scheme),
    east  = ValueBoundaryCondition(0; scheme = tangential_scheme),
    bottom = drag.v_bottom,
    immersed = drag.v_immersed)

boundary_conditions = (u = u_bcs, v = v_bcs,
                       U = tide_bcs.U, V = tide_bcs.V, η = tide_bcs.η)

# ---------------- stratification ----------------
if STRAT != "none"
    prof = glorys_mean_profile(joinpath(@__DIR__, "..", "data"), start_date; region = Symbol(STRAT))
    T_of_z = profile_function(prof.depth, prof.T)
    S_of_z = profile_function(prof.depth, prof.S)
    @info "stratification: horizontally uniform GLORYS '$STRAT' profile" prof.ncolumns deepest_level = last(prof.depth)
    for d in (1, 10, 30, 50, 100, 200, 500, 1000)
        d <= last(prof.depth) + 1 || continue
        @printf("    %5d m   T %6.2f °C   S %6.3f\n", d, T_of_z(-d), S_of_z(-d))
    end

    # Open-boundary tracers: the exterior is the same horizontally uniform profile.
    # NormalRadiation with 1-day inflow nudging (Marchesiello et al. 2001). The
    # inflow CLAMP that produced the 4.4 m/s velocity jet is not a concern here:
    # that was the clamp fighting the split-explicit corrector, and tracers have
    # no corrector to fight.
    tracer_scheme = NormalRadiation(inflow_timescale = 1days, outflow_timescale = Inf)
    exterior(f) = (ξ, z, t) -> f(z)          # a function of the two tangential coordinates + t
    tracer_bcs(f) = FieldBoundaryConditions(
        west  = ValueBoundaryCondition(exterior(f); scheme = tracer_scheme),
        east  = ValueBoundaryCondition(exterior(f); scheme = tracer_scheme),
        south = ValueBoundaryCondition(exterior(f); scheme = tracer_scheme),
        north = ValueBoundaryCondition(exterior(f); scheme = tracer_scheme))
    boundary_conditions = merge(boundary_conditions, (T = tracer_bcs(T_of_z), S = tracer_bcs(S_of_z)))
end

# ---------------- the model ----------------
forcing = BODY_FORCE ? tidal_forcing(harmonics) : NamedTuple()

# A vertically implicit vertical closure is present in EVERY configuration: the
# semi-implicit drag term is solved by the vertical tridiagonal solver and is
# silently inert without one. CATKE is vertically implicit by default.
if CLOSURE == "catke"
    # Cᵇ scales the bottom-boundary-layer mixing length by height above the bottom
    # (Cˢ plays the same role for depth below the surface, default 1.131); the CATKE
    # default 0.28 is Greg Wagner's suggested lever for the inner-shelf skill gap.
    Cᵇ = parse(Float64, get(ENV, "MAB_CB", "0.28"))
    vertical_closure = CATKEVerticalDiffusivity(; mixing_length = CATKEMixingLength(; Cᵇ))
    buoyancy = SeawaterBuoyancy(equation_of_state = NumericalEarth.Oceans.TEOS10EquationOfState())
    tracers  = (:T, :S)                       # :e is added by CATKE itself
    closure  = (HorizontalScalarDiffusivity(ν = ν_h), vertical_closure)
elseif CLOSURE == "kepsilon"
    vertical_closure = TKEDissipationVerticalDiffusivity()
    buoyancy = SeawaterBuoyancy(equation_of_state = NumericalEarth.Oceans.TEOS10EquationOfState())
    tracers  = (:T, :S)                       # :e and :ϵ are added by the closure itself
    # NO horizontal diffusivity paired here: TKEDissipationVerticalDiffusivity does not yet
    # support a Tuple closure (time_step_tke_dissipation_equations! has an explicit "TODO:
    # properly handle closure tuples" — kepsilon_tuple_closure_patch.jl fixes the one part of
    # that gap that IS patchable, top_dissipation_flux, but not this deeper one). So this run
    # carries no ν_h grid-scale damping and is not apples-to-apples with catke/nuz.
    closure  = vertical_closure
elseif CLOSURE == "nuz"
    vertical_closure = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν = ν_z)
    buoyancy = nothing
    tracers  = ()
    closure  = (HorizontalScalarDiffusivity(ν = ν_h), vertical_closure)
else
    error("MAB_CLOSURE must be \"nuz\", \"catke\" or \"kepsilon\", got \"$CLOSURE\"")
end

model = HydrostaticFreeSurfaceModel(grid;
    free_surface = SplitExplicitFreeSurface(grid; substeps = 70),
    coriolis     = HydrostaticSphericalCoriolis(),
    momentum_advection = WENOVectorInvariant(),
    buoyancy, tracers,
    closure,
    forcing, boundary_conditions)

if CLOSURE != "nuz"
    if STRAT == "none"
        set!(model, T = T_uniform, S = S_uniform)
    else
        set!(model, T = (λ, φ, z) -> T_of_z(z), S = (λ, φ, z) -> S_of_z(z))
    end
end

# CATKE's and k-ϵ's diagnostic fields. With a TUPLE closure, closure_fields is a tuple
# aligned with model.closure, so the TKE-based entry has to be picked out by TYPE —
# reaching for `model.closure_fields.κu` directly is a FieldError that, the first time,
# threw away a completed 10-day run's skill analysis.
turbulence_closure_fields(m) = m.closure isa Tuple ?
    m.closure_fields[findfirst(c -> c isa Union{CATKEVerticalDiffusivity, TKEDissipationVerticalDiffusivity},
                               m.closure)] : m.closure_fields

@info "model built" BODY_FORCE CLOSURE

# A fast path for checking the configuration without paying for the full run:
# report what the boundaries impose, take a few steps, and stop. Catching a
# setup error here costs seconds instead of ten minutes.
if get(ENV, "MAB_SETUP_ONLY", "0") == "1"
    println("\n=== what the atlas holds along the boundaries (constituent $(first(CONSTITUENTS))) ===")
    λb = λnodes(grid, Center()) ; φb = φnodes(grid, Center())
    boundary_nodes = (west  = [(first(λb), φ) for φ in φb], east  = [(last(λb), φ) for φ in φb],
                      south = [(λ, first(φb)) for λ in λb], north = [(λ, last(φb)) for λ in λb])
    for (side, nodes) in pairs(boundary_nodes)
        amplitudes = Float64[]
        for (λ, φ) in nodes
            constant, n = tpxo_elevation(window, first(CONSTITUENTS), λ, φ)
            isfinite(real(constant)) && push!(amplitudes, first(amplitude_and_phase(constant)))
        end
        @printf("  %-6s n=%3d wet=%3d   |η| mean %.3f max %.3f m\n",
                side, length(nodes), length(amplitudes), mean(amplitudes), maximum(amplitudes))
    end
    println("\n=== 30 minutes of stepping ===")
    probe = Simulation(model; Δt = 60.0, stop_time = 1800.0)
    run!(probe)
    ηp = model.free_surface.displacement
    up, vp = model.velocities
    @printf("  |η|max %.5f m   |u|max %.4f   |v|max %.4f m/s   all finite: %s\n",
            maximum(abs, interior(ηp)), maximum(abs, up), maximum(abs, vp),
            all(isfinite, interior(ηp)) && all(isfinite, interior(up)))
    if haskey(model.tracers, :T)
        # confirm the initial T/S actually landed (immersed cells excluded via wet columns)
        Tw = filter(x -> x != 0, Array(interior(model.tracers.T)))
        Sw = filter(x -> x != 0, Array(interior(model.tracers.S)))
        @printf("  T range %.2f – %.2f °C   S range %.3f – %.3f   max %s κu %.3e m²/s\n",
                minimum(Tw), maximum(Tw), minimum(Sw), maximum(Sw), CLOSURE,
                maximum(turbulence_closure_fields(model).κu))
    end
    exit(0)
end

# ---------------- run, sampling η hourly for harmonic analysis ----------------
simulation = Simulation(model; Δt = 120.0, stop_time = sim_days * days)
conjure_time_step_wizard!(simulation, cfl = 0.3, max_Δt = 4minutes)

times   = Float64[]
η_store = Matrix{Float64}[]

function sample_eta!(sim)
    push!(times, sim.model.clock.time)
    push!(η_store, Array(interior(sim.model.free_surface.displacement, :, :, 1)))
    return nothing
end
add_callback!(simulation, sample_eta!, TimeInterval(1hours))

# Locate any velocity spike well above the smooth background, hourly -- coarser sampling (the
# 12-hour progress line) only reports the domain max, not where it is. u lives at (Face, Center),
# v at (Center, Face); using each field's own node locations, not the cell-center grid, avoids a
# half-cell misattribution right at a boundary, which is exactly where a spike is expected to be.
const SPIKE_THRESHOLD = 2.0  # m/s, comfortably above this run's smooth background of ~1.3-2.2
λf, φf = λnodes(grid, Face()), φnodes(grid, Face())
λc, φc = λnodes(grid, Center()), φnodes(grid, Center())
zc = znodes(grid, Center())

function report_velocity_spike!(sim)
    u, v, w = sim.model.velocities
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

# Low-pass filtered output (NumericalEarth's LowPassFilter), plus hourly η to check it against.
# MAB_LOWPASS=1 needs a NumericalEarth that has LowPassFilter (the ~/dev/NumericalEarth.jl branch).
if get(ENV, "MAB_LOWPASS", "0") == "1"
    η_out = (; η = model.free_surface.displacement)
    simulation.output_writers[:hourly] = JLD2Writer(model, η_out; dir = OUT_DIR, filename = TAG * "_hourly",
                                                    schedule = TimeInterval(1hours), overwrite_existing = true)
    simulation.output_writers[:daily] = JLD2Writer(model, η_out; dir = OUT_DIR, filename = TAG * "_daily",
                                                   schedule = LowPassFilter(1days), overwrite_existing = true)
    simulation.output_writers[:pentad] = JLD2Writer(model, η_out; dir = OUT_DIR, filename = TAG * "_pentad",
                                                    schedule = LowPassFilter(5days; window = 10days, cutoff = 10days),
                                                    overwrite_existing = true)
end

function progress(sim)
    η = sim.model.free_surface.displacement
    u, v = sim.model.velocities
    @printf("  %s  Δt=%s  |η|max=%.3f m  |u|max=%.2f  |v|max=%.2f m/s\n",
            prettytime(sim), prettytime(sim.Δt),
            maximum(abs, interior(η)), maximum(abs, u), maximum(abs, v))
    return nothing
end
add_callback!(simulation, progress, TimeInterval(12hours))

@info "running $(sim_days) days…"
run!(simulation)
@printf("collected %d hourly samples over %.2f days\n", length(times), last(times)/days)

if CLOSURE != "nuz"
  try
    # The point of this step: compare the closure's COMPUTED viscosity with the uniform
    # 0.1 m²/s probe it replaces. Expect it large in the shelf bottom boundary layer
    # and small in the abyss — the uniform probe gave the abyss a 38 m boundary layer.
    κu = Array(interior(turbulence_closure_fields(model).κu))
    Hc = [static_column_depthᶜᶜᵃ(i, j, grid) for i in 1:Nλ, j in 1:Nφ]
    println("\n$CLOSURE vertical viscosity κu at the end of the run [m²/s]")
    println("  class            columns   median(column max)   median(column mean)")
    for (nm, lo, hi) in (("shelf <50 m", 10.0, 50.0), ("shelf 50–200", 50.0, 200.0),
                         ("slope 200–1000", 200.0, 1000.0), ("deep >1000", 1000.0, 1e9))
        cols = [(i, j) for i in 1:Nλ, j in 1:Nφ if lo < Hc[i, j] <= hi]
        isempty(cols) && continue
        colmax  = [maximum(κu[i, j, :]) for (i, j) in cols]
        colmean = [sum(κu[i, j, :]) / count(>(0), κu[i, j, :] .+ 1e-30) for (i, j) in cols]
        @printf("  %-15s %7d %20.3e %21.3e\n", nm, length(cols), median(colmax), median(colmean))
    end
  catch err
    # A diagnostic must never cost the run its skill analysis.
    @warn "closure viscosity diagnostic failed; continuing to the skill analysis" exception = err
  end
end

# ---------------- MULTI-CONSTITUENT MECHANICS CHECK (MAB_CONSTITUENTS=all) ----------------
# Purpose: confirm that all ten constituents are in and working — NOT a skill
# assessment. A 10-day record cannot separate them by harmonic analysis (S2–K2,
# K1–P1 need 183 days; M2–S2 15), so the M2 analysis below is skipped entirely:
# an M2-only fit to a record that also holds S2 would be biased, not just vague.
# Instead compare model η(t) directly with the full ten-constituent predictions,
# which needs no separation because the forcing is known. The spring–neap envelope
# and the diurnal inequality showing up in the model is what "all in" looks like.
using CairoMakie
if MULTI
    p_rec = reconstruction_parameters(analysis_harmonics)
    λc = λnodes(grid, Center()) ; φc = φnodes(grid, Center())
    Hc = [static_column_depthᶜᶜᵃ(i, j, grid) for i in 1:Nλ, j in 1:Nφ]
    keep = findall(t -> t > n_discard * constituent_period(:M2), times)
    tk = times[keep]
    tpxo_constants(i, j) = begin
        A = Float64[] ; G = Float64[]
        for name in CONSTITUENTS
            z, n = tpxo_elevation(window, name, λc[i], φc[j])
            isfinite(real(z)) || return nothing
            a, g = amplitude_and_phase(z) ; push!(A, a) ; push!(G, g)
        end
        (A, G)
    end
    anomaly(x) = x .- mean(x)

    # one domain-wide guard number, over wet interior cells away from the open boundaries
    buf = 6 ; ratios = Float64[]
    for i in (buf+1):(Nλ-buf), j in (buf+1):(Nφ-buf)
        Hc[i, j] > 10 || continue
        c = tpxo_constants(i, j) ; c === nothing && continue
        m = anomaly([η_store[n][i, j] for n in keep])
        all(isfinite, m) || continue
        pr = [reconstruct(c[1], c[2], t, p_rec) for t in tk]
        push!(ratios, sqrt(mean(abs2, m .- pr)) / sqrt(mean(abs2, pr)))
    end
    @printf("\nALL-10 MECHANICS CHECK — %d wet interior cells, %.2f-day record\n", length(ratios), (tk[end]-tk[1])/days)
    @printf("  rms(model − TPXO 10-constituent prediction) / rms(prediction):  median %.3f   90th pct %.3f\n",
            median(ratios), quantile(ratios, 0.9))

    # gauges: model vs the full ten-constituent TPXO and NOAA predictions
    function nearest_wet_cell(lon, lat; maxr = 4)
        i0 = argmin(abs.(λc .- lon)) ; j0 = argmin(abs.(φc .- lat))
        for r in 0:maxr, di in -r:r, dj in -r:r
            max(abs(di), abs(dj)) == r || continue
            i, j = i0 + di, j0 + dj
            (1 <= i <= Nλ && 1 <= j <= Nφ && Hc[i, j] > 10) && return (i, j)
        end
        return nothing
    end
    gauges = (("8651370", "Duck, NC"), ("8534720", "Atlantic City, NJ"))
    fig = Figure(size = (1400, 560))
    for (row, (sid, name)) in enumerate(gauges)
        hit = nearest_wet_cell(NOAA_HARCON[sid].lon, NOAA_HARCON[sid].lat) ; hit === nothing && continue
        i, j = hit ; c = tpxo_constants(i, j)
        An, Gn = noaa_constants(sid, CONSTITUENTS)
        model = anomaly([η_store[n][i, j] for n in keep])
        tpxo  = [reconstruct(c[1], c[2], t, p_rec) for t in tk]
        noaa  = [reconstruct(An, Gn, t, p_rec) for t in tk]
        @printf("  %-18s  corr(model,TPXO) %.3f  corr(model,NOAA) %.3f   rms diff vs TPXO %.3f m, vs NOAA %.3f m   (signal rms %.3f m)\n",
                name, cor(model, tpxo), cor(model, noaa), sqrt(mean(abs2, model .- tpxo)),
                sqrt(mean(abs2, model .- noaa)), sqrt(mean(abs2, noaa)))
        ax = Axis(fig[row, 1], title = "$name — all 10 constituents", ylabel = "η (m)",
                  xlabel = row == length(gauges) ? "days since $(Dates.format(start_date, "yyyy-mm-dd"))" : "")
        lines!(ax, tk ./ days, model, color = :steelblue, linewidth = 2.5, label = "model")
        lines!(ax, tk ./ days, noaa,  color = :black, linestyle = :dash, linewidth = 2, label = "NOAA (10)")
        lines!(ax, tk ./ days, tpxo,  color = :firebrick, linestyle = :dot, linewidth = 2, label = "TPXO (10)")
        row == 1 && axislegend(ax, position = :rt, framevisible = false)
    end
    Label(fig[0, 1], "All-constituent mechanics check — $(TAG)", fontsize = 18)
    fpath = joinpath(OUT_DIR, "$(TAG)_timeseries.png") ; save(fpath, fig)
    @info "wrote $fpath"
    println("ALL10 CHECK DONE")
    exit(0)
end

# ---------------- diagnostics ----------------
λc = λnodes(grid, Center()) ; φc = φnodes(grid, Center())
H  = Array(interior(bottom_height, :, :, 1))
wet = H .< -10                      # bottom_height is negative in water

T_M2 = constituent_period(:M2)
keep = findall(t -> t > n_discard * T_M2, times)
if length(keep) < 8
    @warn """too few samples after discarding spin-up -- the run ended early.
             Collected $(length(times)) samples spanning $(isempty(times) ? 0 : last(times)/days) days;
             the analysis needs more than $(n_discard*T_M2/days) days. Skipping diagnostics."""
    exit(1)
end
t_fit = times[keep]
@printf("\nharmonic analysis over %.2f days (%d samples), after discarding %.2f days\n",
        (t_fit[end]-t_fit[1])/days, length(keep), n_discard*T_M2/days)

A_mod = fill(NaN, Nλ, Nφ) ; G_mod = fill(NaN, Nλ, Nφ)
A_tpx = fill(NaN, Nλ, Nφ) ; G_tpx = fill(NaN, Nλ, Nφ)

for i in 1:Nλ, j in 1:Nφ
    wet[i, j] || continue
    series = [η_store[n][i, j] for n in keep]
    all(isfinite, series) || continue
    fit = harmonic_analysis(t_fit, series, analysis_harmonics)
    A_mod[i, j] = fit.amplitude[1] ; G_mod[i, j] = fit.phase_lag[1]

    z, n = tpxo_elevation(window, :M2, λc[i], φc[j])
    if isfinite(real(z)) && n > 0
        a, g = amplitude_and_phase(z)
        A_tpx[i, j] = a ; G_tpx[i, j] = rad2deg(g)
    end
end

# Interior only: exclude a buffer next to the open boundaries, where the
# condition is imposed rather than solved, so the comparison measures the
# model's solution and not its own boundary data echoed back.
const buf = 6
interior_mask = falses(Nλ, Nφ)
interior_mask[(buf+1):(Nλ-buf), (buf+1):(Nφ-buf)] .= true
valid = @. interior_mask & isfinite(A_mod) & isfinite(A_tpx)

@printf("comparison over %d wet interior cells (%.0f%% of the domain)\n",
        count(valid), 100*count(valid)/length(valid))

Zm = @. A_mod * cis(-deg2rad(G_mod))
Zt = @. A_tpx * cis(-deg2rad(G_tpx))

dA  = A_mod[valid] .- A_tpx[valid]
dG  = @. mod(G_mod[valid] - G_tpx[valid] + 180, 360) - 180
crms = sqrt(mean(abs2, Zm[valid] .- Zt[valid]))

println()
println("="^74)
println("M2 SKILL vs TPXO10  (interior, open-boundary buffer of $buf cells excluded)")
println("="^74)
@printf("  mean TPXO M2 amplitude          %.4f m   (range %.3f–%.3f)\n",
        mean(A_tpx[valid]), minimum(A_tpx[valid]), maximum(A_tpx[valid]))
@printf("  mean model M2 amplitude         %.4f m   (range %.3f–%.3f)\n",
        mean(A_mod[valid]), minimum(A_mod[valid]), maximum(A_mod[valid]))
@printf("  amplitude bias                  %+.4f m\n", mean(dA))
@printf("  amplitude RMS difference        %.4f m\n", sqrt(mean(abs2, dA)))
@printf("  phase bias                      %+.2f °\n", mean(dG))
@printf("  phase RMS difference            %.2f °\n", sqrt(mean(abs2, dG)))
@printf("  COMPLEX RMS DIFFERENCE          %.4f m   ← the standard skill metric\n", crms)
@printf("  as a fraction of mean amplitude %.3f\n", crms / mean(A_tpx[valid]))

# ---- skill by depth class, and against the bathymetry difference ----
# The domain-wide number hides the structure: if the error is a shelf phenomenon,
# it should separate cleanly by depth. And if it is driven by the model and TPXO
# disagreeing about H, the amplitude ratio should track the depth ratio — Green's
# law for a shoaling wave gives A ∝ H^(-1/4), so a model shelf that is too shallow
# amplifies too much.
H_tpx = fill(NaN, Nλ, Nφ)
for i in 1:Nλ, j in 1:Nφ
    wet[i, j] || continue
    d, n = tpxo_depth(window, λc[i], φc[j])
    n > 0 && d > 0 && (H_tpx[i, j] = d)
end
H_mod = [static_column_depthᶜᶜᵃ(i, j, grid) for i in 1:Nλ, j in 1:Nφ]

println()
println("skill by depth class, and the bathymetry difference in the same cells")
println("  class            n     |  A_mod   A_tpx   cRMS(m)  cRMS/A  |  H_mod   H_tpx  H_mod/H_tpx")
for (name, lo, hi) in (("shelf <50 m", 0.0, 50.0), ("shelf 50–200 m", 50.0, 200.0),
                       ("slope 200–1000", 200.0, 1000.0), ("deep >1000 m", 1000.0, 1e9))
    sel = @. valid & isfinite(H_tpx) & (H_mod > lo) & (H_mod <= hi)
    count(sel) < 20 && continue
    c = sqrt(mean(abs2, Zm[sel] .- Zt[sel]))
    @printf("  %-15s %5d  |  %6.3f  %6.3f  %7.4f  %6.2f  |  %6.0f  %6.0f  %11.3f\n",
            name, count(sel), mean(A_mod[sel]), mean(A_tpx[sel]), c, c/mean(A_tpx[sel]),
            median(H_mod[sel]), median(H_tpx[sel]), median(H_mod[sel] ./ H_tpx[sel]))
end

# Does the amplitude error track the depth disagreement, as Green's law predicts?
shelfsel = @. valid & isfinite(H_tpx) & (H_mod <= 200)
if count(shelfsel) > 50
    ratio_A = A_mod[shelfsel] ./ A_tpx[shelfsel]
    ratio_H = H_mod[shelfsel] ./ H_tpx[shelfsel]
    green   = ratio_H .^ (-0.25)
    @printf("\n  on the shelf: corr(A_mod/A_tpx, (H_mod/H_tpx)^-1/4) = %+.3f over %d cells\n",
            cor(ratio_A, green), count(shelfsel))
    @printf("  median A ratio %.3f   median Green's-law prediction from depth alone %.3f\n",
            median(ratio_A), median(green))
    println("  (a high correlation would implicate the bathymetry difference;")
    println("   a low one means the over-amplification is dynamical, not bathymetric)")
end

# NOAA gauges, where the model has a wet cell
# Nearest WET model cell to a gauge, searching outward — a coastal gauge often
# falls in a land cell at 1/12 deg, and silently comparing against a land cell (or
# skipping the station) would misrepresent the model.
function nearest_wet(lon, lat; maxr = 4)
    i0 = argmin(abs.(λc .- lon)) ; j0 = argmin(abs.(φc .- lat))
    best = nothing ; bestd = Inf
    for r in 0:maxr, di in -r:r, dj in -r:r
        max(abs(di), abs(dj)) == r || continue
        i, j = i0 + di, j0 + dj
        (1 <= i <= Nλ && 1 <= j <= Nφ) || continue
        isfinite(A_mod[i, j]) || continue
        d = di^2 + dj^2
        if d < bestd; bestd = d; best = (i, j, r); end
    end
    return best
end

println()
println("gauge comparison (nearest wet model cell)")
println("  station                model A   NOAA A      TPXO A    model G   NOAA G       ΔG  cells")
gauge_hits = Tuple{String, Int, Int, Float64, Float64}[]
for sid in sort(collect(keys(NOAA_HARCON)))
    r = NOAA_HARCON[sid]
    hit = nearest_wet(r.lon, r.lat)
    hit === nothing && continue
    i, j, off = hit
    A_n, G_n = r.con[:M2]
    ΔG = mod(G_mod[i,j] - G_n + 180, 360) - 180
    @printf("  %-22s %8.3f %9.3f %10.3f %9.1f %8.1f %+8.1f %5d  %s\n",
            r.name, A_mod[i,j], A_n, A_tpx[i,j], G_mod[i,j], G_n, ΔG, off, r.kind)
    push!(gauge_hits, (r.name, i, j, A_n, G_n))
end

# ---------------- cotidal figure ----------------
using CairoMakie
maskland(A) = replace(x -> isfinite(x) ? x : NaN, A)
fig = Figure(size = (1500, 900))
amax = maximum(filter(isfinite, A_tpx))

for (col, (A, G, ttl)) in enumerate(((A_mod, G_mod, "model"), (A_tpx, G_tpx, "TPXO10")))
    ax = Axis(fig[1, col], title = "M2 amplitude, $ttl (m)", xlabel = "°E", ylabel = "°N",
              aspect = DataAspect())
    hm = heatmap!(ax, λc, φc, maskland(A), colormap = :viridis, colorrange = (0, amax))
    Colorbar(fig[2, col], hm, vertical = false)
    ax2 = Axis(fig[3, col], title = "M2 Greenwich phase, $ttl (°)", xlabel = "°E", ylabel = "°N",
               aspect = DataAspect())
    hm2 = heatmap!(ax2, λc, φc, maskland(G), colormap = :twilight, colorrange = (0, 360))
    Colorbar(fig[4, col], hm2, vertical = false)
end

axd = Axis(fig[1, 3], title = "amplitude difference, model − TPXO (m)",
           xlabel = "°E", ylabel = "°N", aspect = DataAspect())
dAf = fill(NaN, Nλ, Nφ) ; dAf[valid] .= A_mod[valid] .- A_tpx[valid]
hmd = heatmap!(axd, λc, φc, dAf, colormap = :balance, colorrange = (-0.15, 0.15))
Colorbar(fig[2, 3], hmd, vertical = false)

axs = Axis(fig[3, 3], title = "M2 amplitude: model vs TPXO", xlabel = "TPXO (m)", ylabel = "model (m)")
scatter!(axs, A_tpx[valid], A_mod[valid], markersize = 3, color = (:steelblue, 0.3))
lines!(axs, [0, amax], [0, amax], color = :black, linestyle = :dash)

Label(fig[0, 1:3], "MAB M2 tide smoke test — $(TAG), $(sim_days) d, body force = $(BODY_FORCE)",
      fontsize = 20)

figpath = joinpath(OUT_DIR, "$(TAG)_cotidal.png")
save(figpath, fig)
@info "wrote $figpath"


# ---------------- proper cotidal chart, and amphidrome comparison ----------------
# A cotidal chart is co-amplitude with CO-PHASE LINES drawn on it; amphidromes are
# where the co-phase lines all converge and the amplitude goes to zero. Contouring
# phase directly is wrong because phase is cyclic — the 0/360 branch cut produces a
# spurious contour straight across the map. Instead, for each phase θ we contour the
# SIGNED angular difference d = mod(G − θ + 180, 360) − 180 at level zero, having
# masked |d| > 90 so the contouring never sees the wrap at all.
signed_phase_diff(G, θ) = @. mod(G - θ + 180, 360) - 180
function cophase_field(G, θ)
    d = signed_phase_diff(G, θ)
    return map(x -> (isfinite(x) && abs(x) < 90) ? x : NaN, d)
end

"Candidate amphidromes: local minima of amplitude below `frac` of the median."
function amphidromes(A, λ, φ; frac = 0.25, halfwidth = 3)
    Nx, Ny = size(A)
    thresh = frac * median(filter(isfinite, A))
    out = Tuple{Float64,Float64,Float64}[]
    for i in (halfwidth+1):(Nx-halfwidth), j in (halfwidth+1):(Ny-halfwidth)
        a = A[i, j]
        (isfinite(a) && a < thresh) || continue
        nb = @view A[(i-halfwidth):(i+halfwidth), (j-halfwidth):(j+halfwidth)]
        all(x -> !isfinite(x) || x >= a, nb) || continue
        push!(out, (λ[i], φ[j], a))
    end
    return out
end

θs = 0:30:330
figc = Figure(size = (1700, 700))
for (col, (A, G, ttl)) in enumerate(((A_mod, G_mod, "model"), (A_tpx, G_tpx, "TPXO10")))
    ax = Axis(figc[1, col], title = "M2 cotidal chart — $ttl", xlabel = "°E", ylabel = "°N",
              aspect = DataAspect())
    hm = heatmap!(ax, λc, φc, maskland(A), colormap = :viridis, colorrange = (0, amax))
    for θ in θs
        contour!(ax, λc, φc, cophase_field(G, θ); levels = [0.0],
                 color = (:white, 0.85), linewidth = 1.2)
    end
    amp = amphidromes(A, λc, φc)
    isempty(amp) || scatter!(ax, [p[1] for p in amp], [p[2] for p in amp];
                             color = :red, marker = :xcross, markersize = 16, strokewidth = 2)
    Colorbar(figc[2, col], hm, vertical = false, label = "amplitude (m)")
end

# Both sets of co-phase lines on one axis: the direct amphidrome comparison.
axo = Axis(figc[1, 3], title = "co-phase every 30°: model (blue) vs TPXO (red)",
           xlabel = "°E", ylabel = "°N", aspect = DataAspect())
heatmap!(axo, λc, φc, maskland(A_tpx), colormap = (:grays, 0.35), colorrange = (0, amax))
for θ in θs
    contour!(axo, λc, φc, cophase_field(G_mod, θ); levels = [0.0], color = (:steelblue, 0.9), linewidth = 1.4)
    contour!(axo, λc, φc, cophase_field(G_tpx, θ); levels = [0.0], color = (:firebrick, 0.9),
             linewidth = 1.4, linestyle = :dash)
end
amp_m = amphidromes(A_mod, λc, φc) ; amp_t = amphidromes(A_tpx, λc, φc)
isempty(amp_m) || scatter!(axo, [p[1] for p in amp_m], [p[2] for p in amp_m];
                           color = :steelblue, marker = :xcross, markersize = 18, strokewidth = 2)
isempty(amp_t) || scatter!(axo, [p[1] for p in amp_t], [p[2] for p in amp_t];
                           color = :firebrick, marker = :cross, markersize = 18, strokewidth = 2)
Label(figc[0, 1:3], "M2 cotidal chart — $(TAG)", fontsize = 20)
cpath = joinpath(OUT_DIR, "$(TAG)_cotidal_chart.png")
save(cpath, figc)
@info "wrote $cpath"

println()
println("amphidrome candidates (amplitude minima below 25% of the median)")
@printf("  model: %d    TPXO: %d\n", length(amp_m), length(amp_t))
for (nm, lst) in (("model", amp_m), ("TPXO ", amp_t))
    for p in lst
        @printf("    %s  (%.2f°E, %.2f°N)  A = %.4f m\n", nm, p[1], p[2], p[3])
    end
end
if !isempty(amp_m) && !isempty(amp_t)
    for pm in amp_m
        d, k = findmin([sqrt(((pm[1]-pt[1])*cosd(pm[2])*111)^2 + ((pm[2]-pt[2])*111)^2) for pt in amp_t])
        @printf("  nearest TPXO amphidrome to model (%.2f, %.2f): %.0f km away\n", pm[1], pm[2], d)
    end
end

# ---------------- save the analysis fields (NetCDF; NCDatasets is already a dep) ----------------

# ---------------- SSH time series at the gauges ----------------
# All three curves are M2 ONLY, which is the only fair comparison: the model is
# forced with M2 alone, so plotting it against a full 37-constituent NOAA
# prediction would show a mismatch that is entirely the missing constituents.
# NOAA and TPXO here are reconstructions from their own published M2 constants
# through the same `reconstruct` and the same astronomy the model was forced with.
p_rec = reconstruction_parameters(analysis_harmonics)
show_gauges = filter(g -> g[1] in ("Duck, NC", "Atlantic City, NJ", "Montauk, NY",
                                   "Sandy Hook, NJ", "Sewells Point, VA"), gauge_hits)
if !isempty(show_gauges)
    ng = length(show_gauges)
    figts = Figure(size = (1400, 260 * ng))
    # last 4 days, so the curves are readable
    tshow = findall(t -> t > last(times) - 4days, times)
    for (row, (name, i, j, A_n, G_n)) in enumerate(show_gauges)
        tt = times[tshow] ./ days
        model = [η_store[n][i, j] for n in tshow]
        noaa  = [reconstruct([A_n], [deg2rad(G_n)], times[n], p_rec) for n in tshow]
        tpxo  = [reconstruct([A_tpx[i,j]], [deg2rad(G_tpx[i,j])], times[n], p_rec) for n in tshow]

        ax = Axis(figts[row, 1],
                  title = @sprintf("%s  —  model %.3f m @ %.0f°, NOAA %.3f m @ %.0f°, TPXO %.3f m @ %.0f°",
                                   name, A_mod[i,j], G_mod[i,j], A_n, G_n, A_tpx[i,j], G_tpx[i,j]),
                  xlabel = row == ng ? "days since $(Dates.format(start_date, "yyyy-mm-dd"))" : "",
                  ylabel = "η (m)")
        lines!(ax, tt, model, color = :steelblue, linewidth = 2.5, label = "model")
        lines!(ax, tt, noaa,  color = :black, linestyle = :dash, linewidth = 2, label = "NOAA M2")
        lines!(ax, tt, tpxo,  color = :firebrick, linestyle = :dot, linewidth = 2, label = "TPXO M2")
        row == 1 && axislegend(ax, position = :rt, framevisible = false)
    end
    Label(figts[0, 1], "M2 sea level at MAB tide gauges — $(TAG)  (M2 only, all three curves)",
          fontsize = 18)
    tspath = joinpath(OUT_DIR, "$(TAG)_timeseries.png")
    save(tspath, figts)
    @info "wrote $tspath"
end
