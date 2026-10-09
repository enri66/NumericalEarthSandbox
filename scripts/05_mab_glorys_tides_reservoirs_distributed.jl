# ==================================================================
# 05 — script 04, run as a domain-decomposed MPI job (one rank per node, x-partition by default)
#
# Physics, boundary conditions, forcing and output are those of 04_mab_glorys_tides_reservoirs.jl. What differs is
# how the grid is held:
#   - Every rank builds the same WHOLE-DOMAIN quantities serially, exactly as 04 does (bathymetry, the GLORYS
#     boundary slabs, tidal constants, sponge and surface-salinity tables), on a whole-domain grid `ggrid`.
#   - The model runs on a Distributed grid; each rank holds a slab of it. Only code that runs inside the model
#     kernels sees LOCAL indices, so it adds this rank's offset (`I_OFF`, `J_OFF`) to reach the whole-domain tables.
#   - The split-explicit free surface gets an explicit substep number from a barotropic CFL estimate, because
#     NumericalEarth's default for distributed grids (fixed from an estimated maximum Δt) is too small over deep
#     trenches at Δt = 5 min and blows up.
# Output and checkpoints are per rank (`_rank{r}` suffix); a restart needs the same number of ranks.
#
# Launch (triton, bundled MPICH):
#   julia --project=<env with MPI> mpiexec.jl -launcher slurm -n R -ppn 1 julia -t 16 --project=<env> 05_mab_....jl
# MAB_STAGE = "bathy" | "model" stops early for testing (see the checks below); default runs the simulation.
#
# What follows is 04's own header.
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

using MPI
using NumericalEarth
using NumericalEarth.DataWrangling: Metadata, MetadataSet, BoundingBox
using NumericalEarth.NestedModels: Interpolated
using CopernicusMarine              # activates the GLORYS download backend
using CopernicusClimateDataStore    # activates the ERA5 download backend
if get(ENV, "MAB_RIVERS", "false") == "true"
    using CDSAPI                    # activates the GloFAS backend (its `download` also returns files already on disk)
end
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: ExponentialDiscretization, MutableVerticalDiscretization, znodes, λnodes, φnodes
using Profile    # MAB_STAGE=profile
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEMixingLength, CATKEEquation
using Oceananigans.BoundaryConditions: PerturbationAdvection, NormalRadiation, ObliqueRadiation, TracerReservoir,
                                       GravityWaveRadiationBoundaryCondition,
                                       SurfaceWaveRadiationBoundaryCondition
using Dates, Printf, Statistics
using NumericalEarth.DataWrangling: jldopen

include(joinpath(@__DIR__, "glorys_bathymetry.jl"))

# ---------------- MPI: ranks, partition, and this rank's offset into the whole domain ----------------
MPI.Init()
const comm   = MPI.COMM_WORLD
const rank   = MPI.Comm_rank(comm)
const nranks = MPI.Comm_size(comm)
say(args...) = rank == 0 && (println(args...); flush(stdout))
allmax(x) = MPI.Allreduce(x, max, comm)
allmin(x) = MPI.Allreduce(x, min, comm)
allsum(x) = MPI.Allreduce(x, +, comm)


const PARTITION_X = parse(Int, get(ENV, "MAB_PARTITION_X", string(nranks)))
const PARTITION_Y = parse(Int, get(ENV, "MAB_PARTITION_Y", "1"))
PARTITION_X * PARTITION_Y == nranks ||
    error("MAB_PARTITION_X × MAB_PARTITION_Y = $(PARTITION_X * PARTITION_Y) must equal the number of MPI ranks ($nranks)")
const STAGE = get(ENV, "MAB_STAGE", "run")

# See 02_mab_glorys_obc.jl for why this defaults to the shared out-of-Dropbox cache.
const DATA_DIR   = get(ENV, "MAB_DATA_DIR", joinpath(homedir(), "Data", "NumericalEarth"))
# Horizontal resolution: MAB_CELLS_PER_DEGREE (default 12, GLORYS's own) model cells per degree, a multiple of 12 so
# that each GLORYS cell holds a whole number of model cells. The domain is the same at every resolution.
const src_resolution  = 1 / 12                                       # GLORYS
const CELLS_PER_DEGREE = parse(Int, get(ENV, "MAB_CELLS_PER_DEGREE", "12"))
CELLS_PER_DEGREE % 12 == 0 || error("MAB_CELLS_PER_DEGREE must be a multiple of 12, got $CELLS_PER_DEGREE")
const resolution = 1 / CELLS_PER_DEGREE
const refinement = CELLS_PER_DEGREE ÷ 12                              # model cells per GLORYS cell, in each direction
# vertical grid: Nz levels to MAB_ZBOTTOM metres (default 5500, below the deepest water, about 5400 m, so the basin keeps
# its true depth; 4000 before 2026-10-05), exponentially stretched so that the top layer is MAB_DZ_TOP metres thick
const Nz     = parse(Int, get(ENV, "MAB_NZ", "50"))
const Δz_top = parse(Float64, get(ENV, "MAB_DZ_TOP", "2"))
const Z_BOTTOM = parse(Float64, get(ENV, "MAB_ZBOTTOM", "5500"))
# MAB_LAND_FRACTION (default 0 = off): a cell is land when more than this fraction of ETOPO's 1/60° points in it are above
# sea level (0.5 is typical), instead of regrid_bathymetry's smoothed-mean-height test, which keeps narrow land such as
# Long Island wet (see land_fraction_mask.jl)
const LAND_FRACTION = parse(Float64, get(ENV, "MAB_LAND_FRACTION", "0"))

const n_pad      = 2                       # GLORYS cells between the GLORYS box and the model's open boundaries
const data_λ     = (-76.0, -64.0)          # GLORYS box (what is on disk)
const data_φ     = ( 34.0,  42.0)
const λ_bounds   = (data_λ[1] + n_pad*src_resolution, data_λ[2] - n_pad*src_resolution)   # model
const φ_bounds   = (data_φ[1] + n_pad*src_resolution, data_φ[2] - n_pad*src_resolution)
const Nλ = round(Int, (λ_bounds[2] - λ_bounds[1]) / resolution)
const Nφ = round(Int, (φ_bounds[2] - φ_bounds[1]) / resolution)
const Nλ_src = round(Int, (data_λ[2] - data_λ[1]) / src_resolution)
const Nφ_src = round(Int, (data_φ[2] - data_φ[1]) / src_resolution)

const start_date = DateTime(get(ENV, "MAB_START_DATE", "2019-04-01"))
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
# MAB_TIDES=false switches the tides off: no TPXO tidal transport and elevation at the open boundaries and no
# equilibrium tidal body force (the GLORYS subtidal exterior is unchanged)
const TIDES = get(ENV, "MAB_TIDES", "true") == "true"
# MAB_TIDE_TRANSPORT_SCALING (default true since 2026-10-05) scales the TPXO tidal transports at the open boundaries by
# H_model / H_TPXO, so the tide enters with TPXO's depth-mean velocity, as the GLORYS subtidal transport does; false
# passes TPXO's transport as is
const TIDE_TRANSPORT_SCALING = get(ENV, "MAB_TIDE_TRANSPORT_SCALING", "true") == "true"
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
const SPONGE_CORNER = get(ENV, "MAB_SPONGE_CORNER", "product")  # "product" (default) | "max" — see sponge_masks
# Surface salinity restoring toward GLORYS as a salt flux with this piston velocity (m/day; 0 turns it off). Over a
# mixed layer of depth h it damps salinity differences in h / piston: ~30 days for a 15 m summer mixed layer,
# ~120 days for 60 m, slow enough to keep the model's eddies and fast enough to hold the seasonal cycle.
const SSS_PISTON   = parse(Float64, get(ENV, "MAB_SSS_PISTON", "0.5")) / days
# Flood ERA5 from the ocean over land before interpolating it to the model (see era5_land_flooding.jl), so coastal
# cells do not take land values: "true" (default) or "false" (ERA5 as delivered)
const ERA5_FLOOD = get(ENV, "MAB_ERA5_FLOOD", "true") == "true"
ERA5_FLOOD && include(joinpath(@__DIR__, "era5_land_flooding.jl"))
# Filling one ERA5 file takes well under a millisecond, so it is redone on every read rather than cached on disk
const ERA5_KW = ERA5_FLOOD ? (; inpainting = ERA5_FLOODING, cache_inpainted_data = false) : (;)

# CATKE parameter changes for sensitivity runs, e.g. MAB_CATKE="Cᵉc=0,Cˢ=0.967": fields of CATKEMixingLength or
# CATKEEquation; the rest keep NumericalEarth's defaults (Oceananigans' calibrated values, with Cᵂϵ = 1)
const CATKE_CHANGES = Dict(Symbol(strip(first(kv))) => parse(Float64, last(kv))
                           for kv in split.(filter(!isempty, split(get(ENV, "MAB_CATKE", ""), ",")), "="))

# checkpoint/restart — PICKUP itself is parsed further down, right before it's used, since it
# can be a Bool, an iteration number, or a filepath (see the comment there)
const CHECKPOINT_EVERY = parse(Float64, get(ENV, "MAB_CHECKPOINT_EVERY", "5")) * days
# For exact-restart debugging: MAB_STOP_ITERATION stops the run after that iteration (default: run to MAB_DAYS), and
# MAB_CHECKPOINT_ITERATIONS writes a checkpoint every that many iterations instead of every MAB_CHECKPOINT_EVERY days.
const STOP_ITERATION        = parse(Float64, get(ENV, "MAB_STOP_ITERATION", "Inf"))
const CHECKPOINT_ITERATIONS = parse(Int, get(ENV, "MAB_CHECKPOINT_ITERATIONS", "0"))

mkpath(DATA_DIR)

# ---------------- grid + bathymetry (as script 02, with a finer vertical grid) ----------------
# The top layer of a right-biased ExponentialDiscretization over depth H is H expm1(H / Nz / h) / expm1(H / h), which
# grows with the scale h; find the h that makes it Δz_top.
function exponential_scale(N, H, Δtop)
    top(h) = H * expm1(H / N / h) / expm1(H / h)
    lo, hi = H / 50, 100H
    for _ in 1:200
        mid = (lo + hi) / 2
        top(mid) > Δtop ? (hi = mid) : (lo = mid)
    end
    return (lo + hi) / 2
end

z = ExponentialDiscretization(Nz, -Z_BOTTOM, 0; scale = exponential_scale(Nz, Z_BOTTOM, Δz_top))

# MAB_BOTTOM=partial: partial bottom cells (PartialCellBottom: the bottom cell of each column is shortened to the real
# depth, keeping at least MAB_PARTIAL_MIN of its full height) instead of whole-cell steps (GridFittedBottom, the default)
const BOTTOM = get(ENV, "MAB_BOTTOM", "gridfitted")
const PARTIAL_MIN = parse(Float64, get(ENV, "MAB_PARTIAL_MIN", "0.2"))
BOTTOM in ("gridfitted", "partial") || error("MAB_BOTTOM must be gridfitted or partial, got $BOTTOM")
model_bottom(h) = BOTTOM == "partial" ? PartialCellBottom(h; minimum_fractional_cell_height = PARTIAL_MIN) : GridFittedBottom(h)
# MAB_ZSTAR=true: a z-star vertical coordinate on the model grid (each column's cells stretch with the free surface by
# (H + η) / H, and the surface freshwater flux adds volume); the whole-domain grid used for the tables stays static
const ZSTAR = get(ENV, "MAB_ZSTAR", "false") == "true"

# The whole-domain grid and bathymetry: the same on every rank, computed serially exactly as in 04.
whole_grid = LatitudeLongitudeGrid(CPU(); size = (Nλ, Nφ, Nz),
                                   longitude = λ_bounds, latitude = φ_bounds, z, halo = (7, 7, 7))
# Every rank regrids the bathymetry itself, with the disk cache off. `regrid_bathymetry` calls `download`, whose `@root` has
# an MPI barrier, only when the cache misses, so with the cache on the number of barriers depends on whether each rank
# sees a cache file: ranks drift out of step (hang) or all regrid and write the file at once (truncated file). With the
# cache off every rank makes exactly one matching barrier call, reads nothing and writes nothing. The cost is a regrid
# at every start (seconds at 1/12 degree, minutes at 1/60 degree).
make_bathymetry() = regrid_bathymetry(whole_grid; dataset = ETOPO2022(), height_above_water = 1, cache = false,
                                      minimum_depth = 10, major_basins = 1, interpolation_passes = 10)
bottom_height = make_bathymetry()
if LAND_FRACTION > 0
    include(joinpath(@__DIR__, "land_fraction_mask.jl"))
    apply_land_fraction!(bottom_height, whole_grid, LAND_FRACTION; minimum_depth = 10, say)
end

# MAB_MASK_CLEANUP=true (default false): clean the coastline (mask_cleanup.jl: water cut off from the open ocean, one-cell
# dead-end bays) and apply the cells you flagged (MAB_MASK_OVERRIDES, a CSV); every rank does the same computation
const MASK_CLEANUP = get(ENV, "MAB_MASK_CLEANUP", "false") == "true"
if MASK_CLEANUP
    include(joinpath(@__DIR__, "mask_cleanup.jl"))
    h = Array(interior(bottom_height))[:, :, 1]
    wet0 = h .< 0
    overrides = isempty(get(ENV, "MAB_MASK_OVERRIDES", "")) ? [] : read_overrides(ENV["MAB_MASK_OVERRIDES"])
    wet = clean_mask(wet0; overrides, λ = collect(λnodes(whole_grid, Center())), φ = collect(φnodes(whole_grid, Center())), say)
    h[wet0 .& .!wet] .= 1.0                       # cells that became land (the height_above_water used above)
    h[.!wet0 .& wet] .= -10.0                     # cells an override made wet: the minimum depth
    set!(bottom_height, h)
    Oceananigans.BoundaryConditions.fill_halo_regions!(bottom_height)
end

if MATCH_BATHY
    match_boundary_bathymetry!(bottom_height, whole_grid, DATA_DIR;
                               region = BoundingBox(longitude = data_λ, latitude = data_φ),
                               n_match = N_MATCH, minimum_depth = 10)
end

ggrid = ImmersedBoundaryGrid(whole_grid, model_bottom(bottom_height))    # whole-domain immersed grid

# This rank's slab of the domain. `I_OFF`, `J_OFF` turn a local index into a whole-domain one.
(Nλ % PARTITION_X == 0 && Nφ % PARTITION_Y == 0) ||
    error("the grid ($Nλ × $Nφ) must divide evenly by the partition ($PARTITION_X × $PARTITION_Y)")
arch = Distributed(CPU(); partition = Partition(x = PARTITION_X, y = PARTITION_Y))
const I_OFF = (arch.local_index[1] - 1) * (Nλ ÷ PARTITION_X)
const J_OFF = (arch.local_index[2] - 1) * (Nφ ÷ PARTITION_Y)

model_z = ZSTAR ? MutableVerticalDiscretization(collect(znodes(whole_grid, Face()))) : z
dist_grid = LatitudeLongitudeGrid(arch; size = (Nλ, Nφ, Nz),
                                  longitude = λ_bounds, latitude = φ_bounds, z = model_z, halo = (7, 7, 7))
local_bottom = Field{Center, Center, Nothing}(dist_grid)
set!(local_bottom, reshape(Array(interior(bottom_height))[I_OFF+1:I_OFF+dist_grid.Nx, J_OFF+1:J_OFF+dist_grid.Ny, 1],
                           dist_grid.Nx, dist_grid.Ny, 1))
Oceananigans.BoundaryConditions.fill_halo_regions!(local_bottom)
grid = ImmersedBoundaryGrid(dist_grid, model_bottom(local_bottom))
say("bottom: $(BOTTOM == "partial" ? "partial cells (minimum fraction $PARTIAL_MIN)" : "whole-cell steps"), vertical coordinate: $(ZSTAR ? "z-star" : "z")")
say("ranks = $nranks (partition $PARTITION_X × $PARTITION_Y), local grid on rank 0: $(dist_grid.Nx) × $(dist_grid.Ny) × $Nz, threads/rank = $(Threads.nthreads())")

# Check that each rank holds the right slab: a position-weighted sum of the bathymetry over all ranks must
# equal the same sum over the whole-domain array.
let whole = Array(interior(bottom_height))[:, :, 1], slab = Array(interior(local_bottom))[:, :, 1]
    reference = sum(i * whole[i, j] for i in 1:Nλ, j in 1:Nφ) + sum(1000 * j * whole[i, j] for i in 1:Nλ, j in 1:Nφ)
    distributed = allsum(sum((I_OFF + i) * slab[i, j] for i in axes(slab, 1), j in axes(slab, 2)) +
                         sum(1000 * (J_OFF + j) * slab[i, j] for i in axes(slab, 1), j in axes(slab, 2)))
    say(@sprintf("bathymetry slab check: whole-domain %.10e, distributed %.10e, %s", reference, distributed,
                 isapprox(reference, distributed; rtol = 1e-12) ? "MATCH" : "MISMATCH"))
end
STAGE == "bathy" && (say("stopping after the grid and bathymetry (MAB_STAGE=bathy)"); exit(0))

# ---------------- GLORYS as FieldTimeSeries on the model grid ----------------
region = BoundingBox(longitude = data_λ, latitude = data_φ)

src_grid = LatitudeLongitudeGrid(CPU(); size = (Nλ_src, Nφ_src, Nz),
                                 longitude = data_λ, latitude = data_φ, z, halo = (7, 7, 7))
dates  = start_date : Day(1) : stop_date
glorys = GLORYSDaily()
meta(name) = Metadata(name; dataset = glorys, dates, dir = DATA_DIR, region)

say("building GLORYS FieldTimeSeries on the whole-domain source grid (downloads/inpaints as needed)…")
build_glorys_series() = (FieldTimeSeries(meta(:u_velocity),  src_grid), FieldTimeSeries(meta(:v_velocity), src_grid),
                         FieldTimeSeries(meta(:temperature), src_grid), FieldTimeSeries(meta(:salinity),   src_grid),
                         FieldTimeSeries(meta(:free_surface), src_grid))
# Every rank builds the series at the same time. NumericalEarth downloads files and writes the inpainted caches on rank 0
# only, inside `@root`, which ends in a barrier on every rank, so the ranks must make the same calls in the same order:
# running the builder on rank 0 first and then on the others paired those barriers wrongly whenever a cache was missing,
# and hung the run.
fts_u, fts_v, fts_T, fts_S, fts_η = build_glorys_series()
say("  done: $(length(fts_u.times)) times, $(Dates.format(start_date,"yyyy-mm-dd")) → $(Dates.format(stop_date,"yyyy-mm-dd"))")

# ---------------- the SUBTIDAL barotropic exterior, from GLORYS (see 02_mab_glorys_obc.jl) ----------------
const zc_src = collect(znodes(src_grid, Center()))
const Δz_src = diff(collect(znodes(src_grid, Face())))
const bh     = Array(interior(bottom_height))[:, :, 1]
const floors = ((bh[1, :], bh[end, :]), (bh[:, 1], bh[:, end]))

# Along a boundary, GLORYS column s covers the model cells (s - n_pad - 1) refinement + 1 … (s - n_pad) refinement;
# its floor is the shallowest of theirs (with one model cell per GLORYS cell, cell s - n_pad)
function source_column_floor(floor_line, s)
    cells = clamp((s - n_pad - 1) * refinement + 1, 1, length(floor_line)):clamp((s - n_pad) * refinement, 1, length(floor_line))
    return maximum(floor_line[cells])
end

function wet_transport(col, floor_line)
    Ns = size(col, 1)
    T, Hwet = zeros(Ns), zeros(Ns)
    for s in 1:Ns
        b = source_column_floor(floor_line, s)
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
say("building the native-resolution GLORYS u/v for the true depth mean… ($Nz_native levels vs the model's $Nz)")
build_native_series() = (FieldTimeSeries(meta(:u_velocity), src_grid_native_z), FieldTimeSeries(meta(:v_velocity), src_grid_native_z))
fts_u_native, fts_v_native = build_native_series()
Hg_src = glorys_deptho_on_grid(src_grid, DATA_DIR; region)
true_floors = ((-Hg_src[1, :], -Hg_src[end, :]), (-Hg_src[:, 1], -Hg_src[:, end]))

function true_depth_mean(col, zc, Δz, floor_line)
    Ns = size(col, 1)
    ubar = zeros(Ns)
    for s in 1:Ns
        b = source_column_floor(floor_line, s)
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
    # model cell i along the boundary sits at GLORYS index n_pad + 1/2 + (i - 1/2) / refinement (n_pad + i at 1/12°)
    σ = clamp(n_pad + 0.5 + (i - 0.5) / refinement, 1, length(a))
    s₁ = clamp(floor(Int, σ), 1, length(a) - 1); r = σ - s₁
    return (1 - w) * ((1 - r) * a[s₁] + r * a[s₁+1]) + w * ((1 - r) * b[s₁] + r * b[s₁+1])
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
say("tides: $(harmonics)")

λᶜ, λᶠ = λnodes(ggrid, Center()), λnodes(ggrid, Face())     # whole-domain node positions
φᶜ, φᶠ = φnodes(ggrid, Center()), φnodes(ggrid, Face())

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
if TIDES && TIDE_TRANSPORT_SCALING
    # TPXO's own depth at the boundary nodes, from its grid file (tpxo.jl, read in a module of its own)
    @eval module TPXOFiles
        include(joinpath($(@__DIR__), "tpxo.jl"))
    end
    tpxo_window = TPXOFiles.load_tpxo((:m2,); λ_bounds = (λ_bounds[1] - 1, λ_bounds[2] + 1),
                                      φ_bounds = (φ_bounds[1] - 1, φ_bounds[2] + 1), dir = TPXO_DIR, verbose = false)
    function scale_to_model_depth!(C, nodes, model_depths)
        ratios = Float64[]
        for (k, (λ, φ)) in enumerate(nodes)
            Ht = first(TPXOFiles.tpxo_depth(tpxo_window, λ, φ)); Hm = model_depths[k]
            (isfinite(Ht) && Ht > 0 && Hm > 0) || continue
            C[k, :] .*= Hm / Ht; push!(ratios, Hm / Ht)
        end
        return ratios
    end
    r = vcat(scale_to_model_depth!(west_U_const,  [(first(λᶠ), φ) for φ in φᶜ], [min(-bh[1, j], Z_BOTTOM) for j in 1:Nφ]),
             scale_to_model_depth!(east_U_const,  [(last(λᶠ), φ) for φ in φᶜ],  [min(-bh[end, j], Z_BOTTOM) for j in 1:Nφ]),
             scale_to_model_depth!(south_V_const, [(λ, first(φᶠ)) for λ in λᶜ], [min(-bh[i, 1], Z_BOTTOM) for i in 1:Nλ]),
             scale_to_model_depth!(north_V_const, [(λ, last(φᶠ)) for λ in λᶜ],  [min(-bh[i, end], Z_BOTTOM) for i in 1:Nλ]))
    say(@sprintf("tidal transports scaled to the model depth at %d boundary nodes: H_model / H_TPXO median %.2f, range %.2f-%.2f",
                 length(r), median(r), minimum(r), maximum(r)))
end
if !TIDES
    for c in (west_U_const, west_η_const, east_U_const, east_η_const, south_V_const, south_η_const, north_V_const, north_η_const)
        c .= 0
    end
    say("tides switched off (MAB_TIDES=false)")
end

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

# The `*_g` functions take the WHOLE-DOMAIN index s along the boundary (j for west/east, i for south/north); the
# functions the model calls receive this rank's LOCAL index, so they add the rank's offset first.
function U_west_g(j, k, grid, clock, f)
    Us, ηsub = lerp_slab(tU, U_lo, j, clock.time), lerp_slab(tU, η_w, j, clock.time)
    Ut, ηt = tidal_UV_eta(west_U_const, west_η_const, harmonics, clock.time, j)
    return (Us + Ut, ηsub + ηt)
end
function U_east_g(j, k, grid, clock, f)
    Us, ηsub = lerp_slab(tU, U_hi, j, clock.time), lerp_slab(tU, η_e, j, clock.time)
    Ut, ηt = tidal_UV_eta(east_U_const, east_η_const, harmonics, clock.time, j)
    return (Us + Ut, ηsub + ηt)
end
function V_south_g(i, k, grid, clock, f)
    Vs, ηsub = lerp_slab(tV, V_lo, i, clock.time), lerp_slab(tV, η_s, i, clock.time)
    Vt, ηt = tidal_UV_eta(south_V_const, south_η_const, harmonics, clock.time, i)
    return (Vs + Vt, ηsub + ηt)
end
function V_north_g(i, k, grid, clock, f)
    Vs, ηsub = lerp_slab(tV, V_hi, i, clock.time), lerp_slab(tV, η_n, i, clock.time)
    Vt, ηt = tidal_UV_eta(north_V_const, north_η_const, harmonics, clock.time, i)
    return (Vs + Vt, ηsub + ηt)
end
U_west(j, k, grid, clock, f)  = U_west_g(j + J_OFF, k, grid, clock, f)
U_east(j, k, grid, clock, f)  = U_east_g(j + J_OFF, k, grid, clock, f)
V_south(i, k, grid, clock, f) = V_south_g(i + I_OFF, k, grid, clock, f)
V_north(i, k, grid, clock, f) = V_north_g(i + I_OFF, k, grid, clock, f)

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
    Utotal :: U   # (s, k, grid, clock, fields) -> (barotropic transport, η), whole-domain index s (the `*_g` functions)
    table  :: T   # times, depth integral of the interpolated profile per time and boundary point, wet depth, floor
    off    :: Int # this rank's offset along the boundary: local index + off = whole-domain index
end

regularize_boundary_condition(c::ConsistentNormalFlow, grid, loc, dim, Side, args...) =
    ConsistentNormalFlow(regularize_boundary_condition(c.base, grid, loc, dim, Side, args...), c.Utotal, c.table, c.off)

@inline function getbc(c::ConsistentNormalFlow, s::Integer, k::Integer, grid::Oceananigans.Grids.AbstractGrid, clock = nothing, args...)
    u = getbc(c.base, s, k, grid, clock, args...)
    tb = c.table
    sg = s + c.off
    1 <= sg <= length(tb.floor) || return u
    Oceananigans.Grids.znode(1, 1, k, grid, Center(), Center(), Center()) > tb.floor[sg] || return u
    t = isnothing(clock) ? 0.0 : clock.time
    n1, n2, w = frame(tb.times, t)
    Uint = (1 - w) * tb.Uint[n1, sg] + w * tb.Uint[n2, sg]
    return u + (c.Utotal(sg, k, grid, (; time = t), nothing)[1] - Uint) / tb.Hwet[sg]
end

function make_consistency_table(floor_line)
    Hwet = [sum(Δz_src[k] for k in 1:Nz if zc_src[k] > b; init = 0.0) for b in floor_line]
    return (; times = collect(fts_u.times), Uint = zeros(length(fts_u.times), length(floor_line)),
            Hwet, floor = collect(floor_line))
end

consistency_tables = (u_west  = make_consistency_table(floors[1][1]), u_east  = make_consistency_table(floors[1][2]),
                      v_south = make_consistency_table(floors[2][1]), v_north = make_consistency_table(floors[2][2]))

normal_condition(fts, Ufun, table, off) = CONSISTENT_UBC ? ConsistentNormalFlow(Interpolated(fts), Ufun, table, off) : Interpolated(fts)

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
    west  = NormalFlowBoundaryCondition(normal_condition(fts_u, U_west_g, consistency_tables.u_west, J_OFF); scheme = normal_scheme),
    east  = NormalFlowBoundaryCondition(normal_condition(fts_u, U_east_g, consistency_tables.u_east, J_OFF); scheme = normal_scheme),
    south = tangential(fts_u),
    north = tangential(fts_u))

v_bcs = FieldBoundaryConditions(
    south = NormalFlowBoundaryCondition(normal_condition(fts_v, V_south_g, consistency_tables.v_south, I_OFF); scheme = normal_scheme),
    north = NormalFlowBoundaryCondition(normal_condition(fts_v, V_north_g, consistency_tables.v_north, I_OFF); scheme = normal_scheme),
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
# MAB_ETA_BC: "chapman" (default, SurfaceWaveRadiation on η at the open boundaries) or "default" (Oceananigans' default
# halo fill for η), for testing
const ETA_BC = get(ENV, "MAB_ETA_BC", "chapman")
η_bcs = ETA_BC == "chapman" ?
    FieldBoundaryConditions(grid, (Center(), Center(), Face());
        west = SurfaceWaveRadiationBoundaryCondition(), east = SurfaceWaveRadiationBoundaryCondition(),
        south = SurfaceWaveRadiationBoundaryCondition(), north = SurfaceWaveRadiationBoundaryCondition()) :
    FieldBoundaryConditions(grid, (Center(), Center(), Face()))

boundary_conditions = (u = u_bcs, v = v_bcs, T = tracer_bcs(fts_T), S = tracer_bcs(fts_S),
                       U = U_bcs, V = V_bcs, η = η_bcs)

# ---------------- GLORYS restoring sponge (MAB_SPONGE_VARS) ----------------
# The mask falls from ≈1 at an open boundary to 0 at SPONGE_WIDTH cells from it as cos²(π d / 2W), where d is
# the distance in cells from the boundary face. A cell belongs to a side's band only if it is connected to that
# side's wet boundary cell along its row (west/east) or column (south/north), so bays and sounds behind land
# are left alone. MAB_SPONGE_CORNER sets how two sides' masks combine near a corner: "product" (default) takes
# 1 - (1 - m₁)(1 - m₂), which is smooth there and a little larger than "max" near the corner; "max" takes the
# larger, i.e. follows the smaller distance, so the mask has a crease along the diagonal into the corner (the
# default before 2026-09-29; "product" was better at both open-open corners in a 90-day test, corner_compare.jl).
function sponge_masks(wet, W)
    Nx, Ny = size(wet)
    # d[i, j, side]: distance in cells from the west, east, south and north boundary faces
    d = fill(Inf, Nx, Ny, 4)
    for j in 1:Ny
        for i in 1:min(W, Nx)
            wet[i, j] || break
            d[i, j, 1] = i - 0.5
        end
        for i in Nx:-1:max(Nx - W + 1, 1)
            wet[i, j] || break
            d[i, j, 2] = Nx - i + 0.5
        end
    end
    for i in 1:Nx
        for j in 1:min(W, Ny)
            wet[i, j] || break
            d[i, j, 3] = j - 0.5
        end
        for j in Ny:-1:max(Ny - W + 1, 1)
            wet[i, j] || break
            d[i, j, 4] = Ny - j + 0.5
        end
    end
    m = [d[i, j, s] < W ? cos(π * d[i, j, s] / 2W)^2 : 0.0 for i in 1:Nx, j in 1:Ny, s in 1:4]
    μᶜ = SPONGE_CORNER == "max" ? dropdims(maximum(m; dims = 3); dims = 3) :
                                  [1 - prod(1 .- m[i, j, :]) for i in 1:Nx, j in 1:Ny]
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
        wet = [wet_node(i, j, k, ggrid, LX, LY) for k in 1:Nz]
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
                    Oceananigans.Fields.interpolate(Oceananigans.Grids.node(i, j, k, ggrid, LX, LY, Center()), frame_field, loc, fts.grid) : 0.0
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
    ci = clamp(i + I_OFF, 1, size(p.column_index, 1))    # whole-domain index of this rank's local (i, j)
    cj = clamp(j + J_OFF, 1, size(p.column_index, 2))
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
forcing = !TIDES ? sponge_forcing :
          merge(sponge_forcing,
                (u = :u in SPONGE_VARS ? (tides.u, sponge_forcing.u) : tides.u,
                 v = :v in SPONGE_VARS ? (tides.v, sponge_forcing.v) : tides.v))
isempty(SPONGE_VARS) || say(@sprintf("GLORYS sponge on %s: %d cells, τ = %.2f days at the boundary, %d wet cells with μ > 0.01",
                                     join(SPONGE_VARS, ","), SPONGE_WIDTH, SPONGE_TAU / days, count(>(0.01), sponge.μᶜ)))

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
            X = Oceananigans.Grids.node(i, j, Nz, ggrid, Center(), Center(), Center())
            target[i, j, n] = Oceananigans.Fields.interpolate(X, frame_field, loc, fts.grid)
        end
    end
    return SurfaceSalinityRestoring(target, times, piston)
end

@inline function getbc(r::SurfaceSalinityRestoring, i::Integer, j::Integer, grid::Oceananigans.Grids.AbstractGrid, clock, fields)
    ci = clamp(i + I_OFF, 1, size(r.target, 1))    # whole-domain index of this rank's local (i, j)
    cj = clamp(j + J_OFF, 1, size(r.target, 2))
    n₁, n₂, w = frame(r.times, clock.time)
    Sᴳ = @inbounds (1 - w) * r.target[ci, cj, n₁] + w * r.target[ci, cj, n₂]
    return @inbounds - r.piston * (Sᴳ - fields.S[i, j, grid.Nz])
end

sss_restoring = SSS_PISTON > 0 ? SurfaceSalinityRestoring(fts_S, SSS_PISTON) : nothing
additional_surface_fluxes = isnothing(sss_restoring) ? NamedTuple() : (; S = sss_restoring)
isnothing(sss_restoring) || say(@sprintf("surface salinity restoring toward GLORYS: piston velocity %.2f m/day", SSS_PISTON * days))

# ---------------- ocean simulation: GLORYS boundaries + equilibrium tidal body force ----------------
# Barotropic substeps for the split-explicit free surface. NumericalEarth's default for a distributed grid fixes
# them from an ESTIMATED maximum Δt (14 at 1/12°), which is too few over 5 km trenches at Δt = 5 min: the free
# surface blows up. Take them from the fastest surface gravity wave (deepest water) and the smallest cell instead.
const Δt_baroclinic = parse(Float64, get(ENV, "MAB_DT", "300"))
function barotropic_substeps(Δt; cfl = parse(Float64, get(ENV, "MAB_BAROTROPIC_CFL", "0.5")))
    R  = 6.371e6
    Δy = R * deg2rad((φ_bounds[2] - φ_bounds[1]) / Nφ)
    Δx = R * deg2rad((λ_bounds[2] - λ_bounds[1]) / Nλ) * cosd(maximum(abs, φ_bounds))
    c  = sqrt(9.80665 * -minimum(bh))
    dτ = cfl / (c * sqrt(1 / Δx^2 + 1 / Δy^2))
    return ceil(Int, Δt / dτ)
end
const SUBSTEPS = parse(Int, get(ENV, "MAB_SUBSTEPS", string(barotropic_substeps(Δt_baroclinic))))
say("split-explicit free surface: $SUBSTEPS barotropic substeps per Δt = $(Δt_baroclinic) s (deepest water $(round(Int, -minimum(bh))) m)")
# MAB_EXTEND_HALOS (default true): with open boundaries Oceananigans substeps the barotropic mode into extended halos and
# refills only the physical boundary halos at each substep; "false" fills all halos (with communication) at each substep
const EXTEND_HALOS = get(ENV, "MAB_EXTEND_HALOS", "true") == "true"
free_surface = SplitExplicitFreeSurface(grid; substeps = SUBSTEPS, extend_halos = EXTEND_HALOS)

function catke_closure(changes)
    ml  = filter(p -> first(p) in fieldnames(CATKEMixingLength), changes)
    tke = filter(p -> first(p) in fieldnames(CATKEEquation), changes)
    unknown = setdiff(keys(changes), keys(ml), keys(tke))
    isempty(unknown) || error("MAB_CATKE: not CATKE parameters: $(join(unknown, ", "))")
    return CATKEVerticalDiffusivity(VerticallyImplicitTimeDiscretization();
                                    mixing_length = CATKEMixingLength(; ml...),
                                    turbulent_kinetic_energy_equation = CATKEEquation(; Cᵂϵ = 1.0, tke...))
end

closure_kw = isempty(CATKE_CHANGES) ? (;) : (; closure = catke_closure(CATKE_CHANGES))
isempty(CATKE_CHANGES) || say("CATKE parameter changes: " * join(["$k = $v" for (k, v) in CATKE_CHANGES], ", "))
# Momentum advection: "default" (NumericalEarth's WENOVectorInvariant: WENO order 9 for the vorticity flux, order 5
# for vertical advection, divergence and the kinetic-energy gradient) or "weno9" (order 9 for all four, less
# dissipative; needs halos of at least 5, the grid has 7)
const MOMENTUM_ADVECTION = get(ENV, "MAB_MOMENTUM_ADVECTION", "default")
MOMENTUM_ADVECTION in ("default", "weno9") || error("MAB_MOMENTUM_ADVECTION must be default or weno9, got $MOMENTUM_ADVECTION")
advection_kw = MOMENTUM_ADVECTION == "weno9" ?
    (; momentum_advection = WENOVectorInvariant(order = 9, time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl = 0.5))) : (;)
MOMENTUM_ADVECTION == "default" || say("momentum advection: WENOVectorInvariant, order 9 throughout")
# Quadratic bottom drag coefficient (NumericalEarth's default 0.003, semi-implicit)
const BOTTOM_DRAG = parse(Float64, get(ENV, "MAB_BOTTOM_DRAG", "0.003"))
BOTTOM_DRAG == 0.003 || say("bottom drag coefficient Cᴰ = $BOTTOM_DRAG")
# MAB_IMPLICIT_DRAG=false: the bottom drag computed explicitly in the tendencies instead of in the vertically implicit step
const IMPLICIT_DRAG = get(ENV, "MAB_IMPLICIT_DRAG", "true") == "true"
IMPLICIT_DRAG || say("bottom drag: explicit")
# River discharge (MAB_RIVERS=true): GloFAS daily discharge at the river mouths inside the domain, deposited on the coastal
# wet cells as a freshwater flux (download_glofas.jl fetches the files into DATA_DIR/glofas). The ocean also gets extra
# vertical mixing in the top MAB_RIVER_MIXING_DEPTH m of the cells receiving a river (river_mouth_vertical_diffusivity),
# so a plume held in one surface cell cannot drive the salinity to zero. The surface salinity restoring (MAB_SSS_PISTON)
# pulls salinity back toward GLORYS and so works against the rivers: turn it down or off together with them.
const RIVERS = get(ENV, "MAB_RIVERS", "false") == "true"
const RIVER_MIXING_DEPTH = parse(Float64, get(ENV, "MAB_RIVER_MIXING_DEPTH", "10"))
# MAB_RIVER_EXTRA (default true) adds the Hudson and the Delaware by hand: GloFAS's automatic mouth detection misses them
# (see glofas_land.jl)
const RIVER_EXTRA = get(ENV, "MAB_RIVER_EXTRA", "true") == "true"
RIVERS && include(joinpath(@__DIR__, "glofas_land.jl"))
# How far, in model cells, a mouth may be from the wet cell that receives it (NumericalEarth's 5 cells is 0.4 degree at
# 1/12 degree), and over how many wet cells nearest the mouth each river's discharge is split equally (NumericalEarth's 8
# is tuned for 1/12 degree). Both scale with CELLS_PER_DEGREE / 12 by default, so the footprint keeps its size in km.
const RIVER_SEARCH_CELLS = parse(Int, get(ENV, "MAB_RIVER_SEARCH_CELLS", string(round(Int, 5 * CELLS_PER_DEGREE / 12))))
const RIVER_SPREAD_CELLS = parse(Int, get(ENV, "MAB_RIVER_SPREAD_CELLS", string(round(Int, 8 * CELLS_PER_DEGREE / 12))))
land = RIVERS ? glofas_land_with_mouths(grid; extra_mouths = RIVER_EXTRA ? MAB_EXTRA_MOUTHS : [], start_date, end_date = stop_date,
                                        dir = joinpath(DATA_DIR, "glofas"), region = BoundingBox(longitude = data_λ, latitude = data_φ),
                                        maximum_spread_cells = RIVER_SPREAD_CELLS, maximum_search_radius = RIVER_SEARCH_CELLS, say,
                                        routing_grid = ggrid, block = (I_OFF, J_OFF, dist_grid.Nx, dist_grid.Ny)) : nothing
RIVERS && say("rivers: GloFAS discharge, each river split over $RIVER_SPREAD_CELLS cells, routed onto the coast, river-mouth mixing over the top $(RIVER_MIXING_DEPTH) m")
river_kw = RIVERS ? (; river_routing = land.river_routing, river_mouth_mixing_depth = RIVER_MIXING_DEPTH) : (;)
ocean = ocean_simulation(grid; free_surface, boundary_conditions, forcing, additional_surface_fluxes, closure_kw..., advection_kw...,
                         river_kw..., bottom_drag_coefficient = BOTTOM_DRAG, implicit_bottom_drag = IMPLICIT_DRAG)

if CONSISTENT_UBC
    # Built on the whole-domain grid, on every rank, so it does not depend on which boundaries a rank owns.
    # (name, whole-domain boundary function, source series, normal dimension, boundary face index)
    faces = ((:u_west,  U_west_g,  fts_u, 1, 1),
             (:u_east,  U_east_g,  fts_u, 1, Nλ + 1),
             (:v_south, V_south_g, fts_v, 2, 1),
             (:v_north, V_north_g, fts_v, 2, Nφ + 1))
    for (name, Utotal, fts, dim, iface) in faces
        tb = consistency_tables[name]
        # Integrate each source frame exactly as the boundary condition samples it (the same node and
        # spatial interpolation), from the frame fields themselves: querying the time series at
        # arbitrary times here would read slices that are not in memory yet.
        loc = Oceananigans.instantiated_location(fts)
        for n in eachindex(tb.times), s in eachindex(tb.Hwet)
            acc = 0.0
            for k in 1:Nz
                zc_src[k] > tb.floor[s] || continue
                X = dim == 1 ? Oceananigans.Grids.node(iface, s, k, ggrid, Face(), Center(), Center()) :
                               Oceananigans.Grids.node(s, iface, k, ggrid, Center(), Face(), Center())
                acc += Δz_src[k] * Oceananigans.Fields.interpolate(X, fts[n], loc, fts.grid)
            end
            tb.Uint[n, s] = acc
        end
        δ = [abs(Utotal(s, 1, ggrid, (; time = 0.0), nothing)[1] - tb.Uint[1, s]) / tb.Hwet[s]
             for s in eachindex(tb.Hwet) if tb.Hwet[s] > 0]
        say(@sprintf("consistent normal velocity, %s: depth-uniform correction at t=0 has max %.2e m/s, mean %.2e m/s",
                     name, maximum(δ), mean(δ)))
    end
end

# Initial conditions from the first GLORYS frame (the one for `start_date`). 04 uses NumericalEarth's
# `set!(model, MetadataSet(...))`, which on a distributed grid reloads the inpainted cache for each rank's LOCAL slab,
# finds the shapes do not match, and OVERWRITES the shared cache files with slab-shaped data (and inpaints each slab
# on its own, so the result depends on the partition). Each rank already holds the whole-domain frames on `src_grid`
# (the same series the boundary conditions use), so interpolate those onto this rank's nodes instead.
function set_from_frame!(field, frame)
    ℓ = Oceananigans.instantiated_location(field)
    source_ℓ = Oceananigans.instantiated_location(frame)
    A = interior(field)
    for k in axes(A, 3), j in axes(A, 2), i in axes(A, 1)
        X = Oceananigans.Grids.node(i, j, k, field.grid, ℓ...)
        A[i, j, k] = Oceananigans.Fields.interpolate(X, frame, source_ℓ, frame.grid)
    end
    Oceananigans.BoundaryConditions.fill_halo_regions!(field)
    return field
end

set_from_frame!(ocean.model.velocities.u, fts_u[1])
set_from_frame!(ocean.model.velocities.v, fts_v[1])
set_from_frame!(ocean.model.tracers.T,    fts_T[1])
set_from_frame!(ocean.model.tracers.S,    fts_S[1])
# Start at GLORYS's own sea level: the Flather exterior values are GLORYS's zos, so removing the domain mean
# here would make the boundaries import or export that mean in the first few hours.
set_from_frame!(ocean.model.free_surface.displacement, fts_η[1])

atmosphere = ERA5PrescribedAtmosphere(; start_date, end_date = stop_date, region,
                                      dir = joinpath(DATA_DIR, "era5"), ERA5_KW...)
radiation  = ERA5PrescribedRadiation(;  start_date, end_date = stop_date, region,
                                      dir = joinpath(DATA_DIR, "era5"), ERA5_KW...)
model = OceanOnlyModel(ocean; atmosphere, radiation, land)
simulation = Simulation(model; Δt = Δt_baroclinic, stop_time = sim_days * days, stop_iteration = STOP_ITERATION)

# MAB_STAGE=model: everything is built and initialised; report the initial state and stop before time stepping.
if STAGE == "model"
    oc0 = ocean.model
    say(@sprintf("initial state:  max|u| = %.8e  max|v| = %.8e  max|η| = %.8e  sum(T) = %.10e  sum(S) = %.10e",
                 allmax(maximum(abs, interior(oc0.velocities.u))), allmax(maximum(abs, interior(oc0.velocities.v))),
                 allmax(maximum(abs, interior(oc0.free_surface.displacement))),
                 allsum(sum(interior(oc0.tracers.T))), allsum(sum(interior(oc0.tracers.S)))))
    say("stopping before time stepping (MAB_STAGE=model)")
    exit(0)
end

# MAB_STAGE=steps: take MAB_STEPS time steps of Δt and write the prognostic fields of this rank after each step listed in
# MAB_DUMP_STEPS (comma-separated, 0 = before the first step) to <tag>_state<n>_rank<r>.jld2, then stop. For comparing layouts
# step by step (compare_layouts.jl): the same case on one rank and split should agree to rounding.
if STAGE == "steps"
    oc = ocean.model
    nsteps = parse(Int, get(ENV, "MAB_STEPS", "10"))
    dumps = parse.(Int, split(get(ENV, "MAB_DUMP_STEPS", "0,1,2,5,10"), ","))
    # MAB_DUMP_SYNC=true: wait for the tracer halo exchanges, which are in flight between steps, before writing. Only the
    # tracers: synchronizing a field whose exchange has completed unpacks its receive buffers again, and for the
    # velocities that would undo the barotropic correction of their halos.
    dump_sync = get(ENV, "MAB_DUMP_SYNC", "false") == "true"
    function dump_state(n)
        if dump_sync
            for fld in (oc.tracers.T, oc.tracers.S, oc.tracers.e)
                Oceananigans.DistributedComputations.synchronize_communication!(fld)
            end
        end
        jldopen("$(TAG)_state$(n)_rank$(rank).jld2", "w") do f
            f["I_OFF"] = I_OFF; f["J_OFF"] = J_OFF
            f["u"] = Array(interior(oc.velocities.u)); f["v"] = Array(interior(oc.velocities.v)); f["w"] = Array(interior(oc.velocities.w))
            f["T"] = Array(interior(oc.tracers.T)); f["S"] = Array(interior(oc.tracers.S))
            f["e"] = Array(interior(oc.tracers.e)); f["eta"] = Array(interior(oc.free_surface.displacement))
            # whole arrays with their halos, and the halo widths, to compare halo cells across layouts
            fs = oc.free_surface
            for (name, fld) in (("eta", fs.displacement), ("U", fs.barotropic_velocities.U), ("V", fs.barotropic_velocities.V),
                                ("u", oc.velocities.u), ("v", oc.velocities.v), ("T", oc.tracers.T))
                f[name * "_p"] = Array(parent(fld)); f[name * "_halo"] = collect(Oceananigans.Grids.halo_size(fld.grid))
            end
            # the net surface fluxes handed to the ocean, and the air-sea momentum flux with its halo
            nf = model.interfaces.net_fluxes.ocean
            f["tau_x"] = Array(interior(nf.u)); f["tau_y"] = Array(interior(nf.v)); f["J_T"] = Array(interior(nf.T))
            # every 2D input of the air-sea flux kernel, with halos: the exchanger's ocean and atmosphere states
            for (side, st) in (("ocn", model.interfaces.exchanger.ocean.state), ("atm", model.interfaces.exchanger.atmosphere.state))
                for (n, fld) in pairs(st)
                    fld isa Oceananigans.Fields.AbstractField || continue
                    a = Array(parent(fld)); ndims(a) == 3 && size(a, 3) > 1 && (a = a[:, :, end:end])
                    f[side * "_" * string(n) * "_p"] = a; f[side * "_" * string(n) * "_halo"] = collect(Oceananigans.Grids.halo_size(fld.grid))
                end
            end
            ao = model.interfaces.atmosphere_ocean_interface.fluxes
            f["ao_x_p"] = Array(parent(ao.x_momentum)); f["ao_x_halo"] = collect(Oceananigans.Grids.halo_size(ao.x_momentum.grid))
            # the surface flux each prognostic field receives (top boundary condition), where it is a field
            for (name, fld) in (("u", oc.velocities.u), ("v", oc.velocities.v), ("T", oc.tracers.T), ("S", oc.tracers.S))
                q = fld.boundary_conditions.top.condition
                q isa Oceananigans.Fields.AbstractField && (f["Jtop_" * name] = Array(interior(q)))
            end
        end
    end
    # MAB_SYNC_HALOS=true: fill the ocean fields' halos and wait for the exchange to finish before each step
    sync_halos = get(ENV, "MAB_SYNC_HALOS", "false") == "true"
    function sync_ocean_halos!()
        for f in (oc.velocities.u, oc.velocities.v, oc.tracers.T, oc.tracers.S, oc.tracers.e)
            Oceananigans.BoundaryConditions.fill_halo_regions!(f)
            Oceananigans.DistributedComputations.synchronize_communication!(f)
        end
    end
    0 in dumps && dump_state(0)
    for n in 1:nsteps
        sync_halos && sync_ocean_halos!()
        time_step!(model, Δt_baroclinic)
        n in dumps && dump_state(n)
    end
    say("took $nsteps steps (MAB_STAGE=steps)")
    exit(0)
end

# MAB_STAGE=profile: after MAB_WARMUP steps, time MAB_STEPS steps, then sample MAB_PROFILE_STEPS more with the profiler;
# rank 0 writes the samples as a flat list, by count, to <tag>_profile_rank0.txt
if STAGE == "profile"
    warmup = parse(Int, get(ENV, "MAB_WARMUP", "3"))
    nsteps = parse(Int, get(ENV, "MAB_STEPS", "10"))
    nprofile = parse(Int, get(ENV, "MAB_PROFILE_STEPS", "5"))
    for _ in 1:warmup
        time_step!(model, Δt_baroclinic)
    end
    t₀ = time_ns()
    for _ in 1:nsteps
        time_step!(model, Δt_baroclinic)
    end
    say(@sprintf("wall time per step: %.3f s (mean of %d steps after %d warm-up steps)", (time_ns() - t₀) / 1e9 / nsteps, nsteps, warmup))
    Profile.init(n = 10^8, delay = 0.005)
    Profile.@profile for _ in 1:nprofile
        time_step!(model, Δt_baroclinic)
    end
    if rank == 0
        open("$(TAG)_profile_rank0.txt", "w") do io
            Profile.print(IOContext(io, :displaysize => (100000, 400)); format = :flat, sortedby = :count, mincount = 10)
        end
    end
    say("profiled $nprofile steps (MAB_STAGE=profile)")
    exit(0)
end

wall = Ref(time())
function progress(sim)
    u, v, w = sim.model.ocean.model.velocities
    η = sim.model.ocean.model.free_surface.displacement
    date = start_date + Second(round(Int, sim.model.clock.time))
    ηwet = filter(isfinite, interior(η))
    umax, vmax = allmax(maximum(abs, interior(u))), allmax(maximum(abs, interior(v)))
    ηmin, ηmax = allmin(minimum(ηwet)), allmax(maximum(ηwet))
    S = sim.model.ocean.model.tracers.S
    Smin, Smax = allmin(minimum(filter(isfinite, interior(S)))), allmax(maximum(filter(isfinite, interior(S))))
    say(@sprintf("%s  |u|=%.2f |v|=%.2f  η∈[%+.2f,%+.2f]  S∈[%.2f,%.2f]  (%.0f s)",
                 Dates.format(date, "yyyy-mm-dd HH:MM"), umax, vmax, ηmin, ηmax, Smin, Smax, time() - wall[]))
    wall[] = time()
    return nothing
end
# MAB_PROGRESS_MINUTES (default 360) sets how often the progress line is printed; shorten it to watch a run that may blow up
const PROGRESS_MINUTES = parse(Float64, get(ENV, "MAB_PROGRESS_MINUTES", "360"))
add_callback!(simulation, progress, TimeInterval(PROGRESS_MINUTES * 60))

const λf = λnodes(ggrid, Face());   const λc = λnodes(ggrid, Center())     # whole-domain node positions
const φf = φnodes(ggrid, Face());   const φc = φnodes(ggrid, Center())
const zc = znodes(ggrid, Center())
function report_velocity_spike!(sim)
    u, v, w = sim.model.ocean.model.velocities
    ua, va = Array(interior(u)), Array(interior(v))

    # the rank holding the global maximum reports it, with whole-domain indices
    iu = argmax(abs.(ua))
    if allmax(abs(ua[iu])) > SPIKE_THRESHOLD && abs(ua[iu]) == allmax(abs(ua[iu]))
        i, j, k = iu.I; i += I_OFF; j += J_OFF
        @printf("  SPIKE u=%+.2f m/s at (i=%d,j=%d,k=%d) λ=%.3f φ=%.3f z=%.1f m  t=%s\n",
                ua[iu], i, j, k, λf[i], φc[j], zc[k], prettytime(sim))
        # the whole column there: u and v at a few levels (k = 1 is the bottom cell, Nz the surface) and the depth mean
        # of u over the wet cells, to tell a depth-uniform (barotropic) flow from a bottom-intensified one
        col, vcol = ua[iu.I[1], iu.I[2], :], va[iu.I[1], iu.I[2], :]
        dz = diff(znodes(ggrid, Face())); wetk = findall(!iszero, col)
        levels = unique(clamp.([1, 2, 3, 5, 10, 20, 40, 60, 80, Nz], 1, Nz))
        @printf("    column at that cell: depth-mean u = %+.3f m/s over %d wet levels; u(k) %s; v(k) %s\n",
                sum(col[wetk] .* dz[wetk]) / sum(dz[wetk]), length(wetk),
                join([@sprintf("%d:%+.2f", k, col[k]) for k in levels], " "), join([@sprintf("%d:%+.2f", k, vcol[k]) for k in levels], " "))
    end

    iv = argmax(abs.(va))
    if allmax(abs(va[iv])) > SPIKE_THRESHOLD && abs(va[iv]) == allmax(abs(va[iv]))
        i, j, k = iv.I; i += I_OFF; j += J_OFF
        @printf("  SPIKE v=%+.2f m/s at (i=%d,j=%d,k=%d) λ=%.3f φ=%.3f z=%.1f m  t=%s\n",
                va[iv], i, j, k, λc[i], φf[j], zc[k], prettytime(sim))
    end
    return nothing
end
add_callback!(simulation, report_velocity_spike!, TimeInterval(1hours))

# MAB_MOORINGS: hourly model columns at fixed positions, for comparison with moored current profilers. "pioneer" (the
# OOI Coastal Pioneer New England Shelf moorings, 2019 positions) or "name:λ:φ,...". Each position takes the nearest
# wet cell. The rank that owns the cell keeps u and v (averaged to the cell centre), T, S, CATKE's κc and κu (at cell
# faces) and TKE e, and at the surface the atmospheric wind the fluxes are computed from (ua, va), the turbulent momentum
# fluxes (τx, τy: NumericalEarth's x_momentum, y_momentum, N/m², positive upward) and the friction velocity u★, and
# rewrites <tag>_moorings_rank<r>.jld2 once a simulated day and at the end.
const PIONEER_MOORINGS = "OSSM:-70.8869:39.9375,PMUO:-70.7702:39.9393,PMCO:-70.8792:40.0968,CNSM:-70.7783:40.1333"
const MOORINGS = let s = get(ENV, "MAB_MOORINGS", "")
    s = s == "pioneer" ? PIONEER_MOORINGS : s
    [(String(a), parse(Float64, b), parse(Float64, c)) for (a, b, c) in split.(filter(!isempty, split(s, ",")), ":")]
end

if !isempty(MOORINGS)
    JLD2m = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "JLD2")]
    wet_columns = [(i, j) for i in 1:Nλ, j in 1:Nφ if bh[i, j] < 0]
    columns = map(MOORINGS) do (name, λ, φ)
        i, j = wet_columns[argmin([hypot((λc[i] - λ) * cosd(φ), φc[j] - φ) for (i, j) in wet_columns])]
        say(@sprintf("mooring %s (%.4f, %.4f): cell (%d, %d) at (%.4f, %.4f), model depth %.0f m", name, λ, φ, i, j, λc[i], φc[j], -bh[i, j]))
        (; name, λ, φ, i, j, depth = -bh[i, j])
    end
    mine = filter(c -> I_OFF < c.i <= I_OFF + dist_grid.Nx && J_OFF < c.j <= J_OFF + dist_grid.Ny, columns)
    mooring_record = Dict(c.name => Dict(:time => Float64[], :u => Vector{Float64}[], :v => Vector{Float64}[],
                                         :T => Vector{Float64}[], :S => Vector{Float64}[], :κc => Vector{Float64}[],
                                         :κu => Vector{Float64}[], :e => Vector{Float64}[],
                                         :ua => Float64[], :va => Float64[], :τx => Float64[], :τy => Float64[], :u★ => Float64[]) for c in mine)
    mooring_file = joinpath(@__DIR__, "..", TAG * "_moorings_rank$(rank).jld2")

    function save_moorings()
        isempty(mine) && return nothing
        JLD2m.jldopen(mooring_file, "w") do file
            file["z_center"] = collect(zc); file["z_face"] = collect(znodes(ggrid, Face()))
            file["start_date"] = string(start_date)
            for c in mine
                r = mooring_record[c.name]
                file["$(c.name)/position"] = (c.λ, c.φ); file["$(c.name)/cell"] = (c.i, c.j, λc[c.i], φc[c.j], c.depth)
                file["$(c.name)/time"] = r[:time]
                for q in (:u, :v, :T, :S, :κc, :κu, :e)
                    file["$(c.name)/$q"] = isempty(r[q]) ? zeros(0, 0) : reduce(hcat, r[q])
                end
                for q in (:ua, :va, :τx, :τy, :u★)
                    file["$(c.name)/$q"] = r[q]
                end
            end
        end
        return nothing
    end

    function record_moorings!(sim)
        oc = sim.model.ocean.model
        u, v = oc.velocities.u, oc.velocities.v
        T, S = oc.tracers.T, oc.tracers.S
        # with rivers the closure is a tuple (CATKE plus the river-mouth diffusivity): take CATKE's fields
        catke_fields = oc.closure_fields isa Tuple ? first(filter(f -> hasproperty(f, :κc), oc.closure_fields)) : oc.closure_fields
        κc = catke_fields.κc; κu = catke_fields.κu; e = oc.tracers.e
        atm = sim.model.interfaces.exchanger.atmosphere.state
        ao  = sim.model.interfaces.atmosphere_ocean_interface.fluxes
        for c in mine
            i, j = c.i - I_OFF, c.j - J_OFF
            r = mooring_record[c.name]
            push!(r[:time], oc.clock.time)
            push!(r[:u], [(u[i, j, k] + u[i+1, j, k]) / 2 for k in 1:Nz])
            push!(r[:v], [(v[i, j, k] + v[i, j+1, k]) / 2 for k in 1:Nz])
            push!(r[:T], [T[i, j, k] for k in 1:Nz]); push!(r[:S], [S[i, j, k] for k in 1:Nz])
            push!(r[:κc], [κc[i, j, k] for k in 1:Nz+1]); push!(r[:κu], [κu[i, j, k] for k in 1:Nz+1])
            push!(r[:e], [e[i, j, k] for k in 1:Nz])
            push!(r[:ua], atm.u[i, j, 1]); push!(r[:va], atm.v[i, j, 1])
            push!(r[:τx], ao.x_momentum[i, j, 1]); push!(r[:τy], ao.y_momentum[i, j, 1]); push!(r[:u★], ao.friction_velocity[i, j, 1])
        end
        return nothing
    end
    add_callback!(simulation, record_moorings!, TimeInterval(1hours))
    add_callback!(simulation, sim -> save_moorings(), TimeInterval(1days))
end

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
# MAB_BAROTROPIC_OUTPUT=true: hourly barotropic transports U, V (m²/s) of the split-explicit free surface, for
# tidal-current maps
if get(ENV, "MAB_BAROTROPIC_OUTPUT", "false") == "true"
    simulation.output_writers[:barotropic] = JLD2Writer(oc, oc.free_surface.barotropic_velocities;
        filename = joinpath(@__DIR__, "..", TAG * "_barotropic.jld2"),
        schedule = TimeInterval(1hours), overwrite_files = fresh_start)
end
simulation.output_writers[:eta_daily] = JLD2Writer(oc, η_out;
    filename = joinpath(@__DIR__, "..", TAG * "_eta_daily.jld2"),
    schedule = FilteredTimeInterval(LanczosKernel(6days; cutoff = 40hours); interval = 1days), overwrite_files = fresh_start)
simulation.output_writers[:eta_pentad] = JLD2Writer(oc, η_out;
    filename = joinpath(@__DIR__, "..", TAG * "_eta_pentad.jld2"),
    schedule = FilteredTimeInterval(LanczosKernel(10days; cutoff = 10days); interval = 5days), overwrite_files = fresh_start)
# MAB_CATKE_OUTPUT=true: daily means of CATKE's tracer diffusivity and its shear and convective parts, with N², S²,
# the TKE and the surface buoyancy flux (catke_diagnostics.jl), sampled every 4 time steps
if get(ENV, "MAB_CATKE_OUTPUT", "false") == "true"
    include(joinpath(@__DIR__, "catke_diagnostics.jl"))
    simulation.output_writers[:catke_daily] = JLD2Writer(oc, catke_diagnostics(oc);
        filename = joinpath(@__DIR__, "..", TAG * "_catke_daily.jld2"),
        schedule = AveragedTimeInterval(1days; stride = 4), overwrite_files = fresh_start)
end

simulation.output_writers[:checkpointer] = Checkpointer(model;
    schedule = CHECKPOINT_ITERATIONS > 0 ? IterationInterval(CHECKPOINT_ITERATIONS) : TimeInterval(CHECKPOINT_EVERY), dir = checkpoint_dir,
    prefix = basename(TAG) * "_checkpoint", overwrite_files = false, cleanup = false)

say("running on $nranks ranks: tangential=$TANGENTIAL  Uᵉˣᵗ=$UEXT_MODE  τ_in=$(TAU_IN/86400) day(s)  " *
    "tides=$(join(TIDE_CONSTITUENTS, ",")) velocity=$VELOCITY_SCHEME$(VELOCITY_SCHEME == "oblique" ? "(τ_in=$(OBLIQUE_TAU_IN/days)d, τ_out=$(OBLIQUE_TAU_OUT/days)d, w=$PHASE_SPEED_WEIGHT)" : "") consistent_ubc=$CONSISTENT_UBC tracers=$TRACER_SCHEME reservoir(L_in=$RESERVOIR_L_IN, L_out=$RESERVOIR_L_OUT)  " *
    "$(sim_days) days  pickup=$PICKUP")
# MAB_PROFILE_RUN=true: as MAB_STAGE=profile, but stepping the whole simulation (callbacks and output writers included)
if get(ENV, "MAB_PROFILE_RUN", "false") == "true"
    warmup = parse(Int, get(ENV, "MAB_WARMUP", "3"))
    nsteps = parse(Int, get(ENV, "MAB_STEPS", "10"))
    nprofile = parse(Int, get(ENV, "MAB_PROFILE_STEPS", "5"))
    for _ in 1:warmup
        time_step!(simulation)
    end
    t₀ = time_ns()
    block = parse(Int, get(ENV, "MAB_PROFILE_BLOCK", string(nsteps)))
    for b in 1:cld(nsteps, block)
        stats = @timed for _ in 1:block
            time_step!(simulation)
        end
        say(@sprintf("steps %4d-%4d: %.3f s per step, %.1f MB allocated per step, %.0f%% of the time in GC",
                     (b - 1) * block + 1, b * block, stats.time / block, stats.bytes / block / 1e6, 100 * stats.gctime / stats.time))
    end
    say(@sprintf("wall time per simulation step: %.3f s (mean of %d steps after %d warm-up steps)", (time_ns() - t₀) / 1e9 / (cld(nsteps, block) * block), cld(nsteps, block) * block, warmup))
    Profile.init(n = 10^8, delay = 0.005)
    Profile.@profile for _ in 1:nprofile
        time_step!(simulation)
    end
    if rank == 0
        open("$(TAG)_profile_run_rank0.txt", "w") do io
            Profile.print(IOContext(io, :displaysize => (100000, 400)); format = :flat, sortedby = :count, mincount = 10)
        end
    end
    say("profiled $nprofile simulation steps (MAB_PROFILE_RUN)")
    exit(0)
end

run!(simulation; pickup = PICKUP, checkpoint_at_end = true)
isempty(MOORINGS) || save_moorings()
say("\n✅ done — $(TAG)")
