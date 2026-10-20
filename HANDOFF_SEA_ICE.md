# Handoff: sea ice in the regional NumericalEarth setup (session "NE-Ice")

Written 2026-10-14 from the main MAB session. Goal: add a coupled sea-ice model to a regional simulation, starting from
the production MAB configuration, and validate it against observations.

## Where things stand

- The regional ocean setup is `scripts/05_mab_glorys_tides_reservoirs_distributed.jl` (MPI, per-rank output), run on
  triton. Agreed production settings and open items: `PRODUCTION_CONFIG.md` (z*, partial cells, fixed open boundaries,
  bounds-preserving WENO-7 for T/S, GloFAS rivers, SSS restoring with a river mask).
- Open boundaries: oblique radiation for baroclinic velocities, Flather-type barotropic conditions, tracer reservoirs,
  GLORYS T/S/u/v sponge, tides. These are working; the related Oceananigans PRs are #6091, #6092, #6182–#6185.
- Mixed-layer depth bias work with Greg Wagner is ongoing in the main session (not needed here).

## What NumericalEarth provides

- `sea_ice_simulation(grid, ocean; ...)` (`src/SeaIces/sea_ice_simulation.jl` in NumericalEarth.jl) builds a
  ClimaSeaIce `SeaIceModel` on the ocean grid:
  - slab ice thermodynamics with a conductive internal heat flux (conductivity 2 W/m/K), a snow slab
    (conductivity 0.31), ice salinity 4 psu, ice density 900 kg/m³, snow 330 kg/m³, consolidation thickness 0.05 m;
  - bottom boundary at the freezing point of the ocean surface salinity (`IceWaterThermalEquilibrium`);
  - dynamics from `sea_ice_dynamics`: elastic-viscous-plastic rheology, semi-implicit ice-ocean stress
    (drag 3.24×10⁻³), wind stress, Coriolis, ocean surface tilt, split-explicit solver; `SplitRungeKutta3` stepping.
- The Earth-system model couples ocean, sea ice, atmosphere and radiation: ice-ocean heat and salt fluxes, the
  freezing-limited ocean temperature, and ice-concentration-weighted air-sea fluxes.
- Data: GLORYS has `siconc`, `sithick`, `usi`, `vsi` (`:sea_ice_concentration`, `:sea_ice_thickness`,
  `:sea_ice_u_velocity`, `:sea_ice_v_velocity`); ECCO has `SIarea`, `SIheff`. `set!` can initialise the ice from
  metadata, and the restoring machinery accepts `:sea_ice_thickness` and `:sea_ice_area_fraction`.
- Reference uses: NumericalEarth tests `test_ocean_sea_ice_model.jl`, `test_sea_ice_ocean_heat_fluxes.jl`,
  `test_snow_model_integration.jl`, `test_ocean_sea_ice_oceanbiome.jl`; examples `one_degree_simulation.jl`,
  `global_climate_simulation.jl`.

## The main decision: the region

The MAB box (34–42°N, 76–64°W) has essentially no sea ice, so it cannot test the coupling. A new or extended domain
is needed. Candidates to discuss with Enri:

- Gulf of St. Lawrence (seasonal first-year ice, strong river input, well observed).
- Scotian Shelf / Newfoundland–Labrador shelf (ice edge, Labrador Current).
- An Arctic shelf sea.

A new region means new bathymetry (`regrid_bathymetry` from ETOPO), a land mask clean-up, the open-boundary set-up
(the OBC code is region-independent), and new GLORYS/ERA5/GloFAS downloads for the region and period.

## Decision (Enri, 2026-10-14)

- Region: the Gulf of St. Lawrence and the Newfoundland–Labrador shelf in one box, about 70.5–48°W, 45–56°N.
- Period: winter 2013–14, starting ice-free on 2013-12-01 from GLORYS ocean fields, to mid/late May 2014.
- Development resolution: 1/12°, 50 levels; finer once the coupling works.

GLORYS sea ice in that box (`siconc`, `sithick`, `usi`, `vsi`; daily files on antares A3):

| date | ice area (10³ km²) | ice volume (km³) | southward ice volume flux across 56°N (km³/day) |
|---|---|---|---|
| 2013-12-01 | 0.8 | 0.0 | 0.0 |
| 2013-12-31 | 61 | 33 | 0.7 |
| 2014-01-30 | 258 | 132 | 2.3 |
| 2014-03-01 | 381 | 219 | 3.6 |
| 2014-03-16 | 403 | 249 | 0.7 |
| 2014-04-15 | 177 | 142 | 0.8 |
| 2014-05-15 | 80 | 77 | 1.1 |
| 2014-05-30 | 28 | 20 | 2.3 |

The ice entering across the northern boundary (Labrador shelf) is of the order of the ice volume in the box over the
winter, so ice open boundaries (sponge restoring of thickness and concentration at least, possibly ice velocity) are
needed from the start, not as a refinement. The eastern boundary carries almost none.

Other winters checked (Gulf of St. Lawrence part, area 10³ km² on Feb 15 / Mar 10 / Apr 10): 2003 177/212/86,
2008 143/177/30, 2014 171/159/99, 2015 162/152/30, 2019 115/138/7, 2022 53/115/6, 2024 18/8/0.6 (a contrasting
low-ice year for later).

The data box is 71.5–48°W, 45–56°N (west edge moved from 70.5°W so the St. Lawrence estuary closes near Québec
instead of crossing an open boundary). The 45°N southern edge crosses the Bay of Fundy; plan to mask Fundy/Gulf of Maine
as land (e.g. water west of 64.2°W and south of 45.95°N) rather than open a boundary through its tides.

## Data (staged 2026-10-14)

Cut from the antares archives with `scripts/regional_cache_from_archives.jl` (checked against NumericalEarth's own
downloads for the MAB box: GLORYS and GloFAS identical, ERA5 identical up to float32 rounding), 2013-11-25 to
2014-06-05, and synced to triton's `/t0/workdir/enrique/Data/NumericalEarth` (region suffix `_-71.5_-48.0_45.0_56.0`):
GLORYS daily thetao, so, uo, vo, zos, siconc, sithick, usi, vsi; ERA5 hourly u10, v10, t2m, d2m, sp, tp, sf, ssrd, strd;
GloFAS daily discharge; GLORYS static deptho and the ERA5 land-sea mask (`scripts/stage_static_region_data.jl`, and
on triton `fetch_lsm_gsl.jl`: antares' CDS file fetch fails with an SSL certificate error).

## Distributed sea ice (2026-10-14, `restart_test_harness/distributed_sea_ice.jl`, triton `runs/seaice`)

Coupled ocean_simulation + sea_ice_simulation + prescribed cold, windy atmosphere, 80 × 40 × 10 at 1/12° near 50°N,
land strip + seamount, an ice patch across the rank edges, 4 steps of 300 s, distributed vs serial:

| case | ocean u, v, η, T, S | ice h, ℵ | ice u, v |
|---|---|---|---|
| 1 × 1 (Distributed arch, one rank) | exact | exact | exact |
| ocean only (no sea ice), 2 × 1 | exact | – | – |
| ice thermodynamics only, 2 × 1 | T 2.5e-2 °C, u 8e-4 m/s at the rank edge | 1.5e-9 | – |
| full (EVP, 8 substeps), 2 × 1 | as above | h 4.7e-2 m | 0.31 of 0.41 m/s, at the rank edge |
| full (EVP, 8 substeps), 1 × 2 | T 1.1e-2 °C | h 5.0e-2 m | 0.26 of 0.35 m/s, at the rank edge |

Findings:
- The ocean differs at rank edges as soon as sea ice is coupled, even without ice dynamics; ocean-only is exact. Likely an
  un-exchanged halo in NumericalEarth's ocean–sea-ice flux path (the same family as NumericalEarth #746 and the CATKE
  Jᵇ halo of Oceananigans #6184).
- ClimaSeaIce's split-explicit EVP on a distributed grid gives wrong ice velocities at rank edges.
- Its distributed design widens the velocity halos to 2·substeps + 3 cells and does not communicate during substeps;
  with NumericalEarth's default 240 substeps (RK3 ocean) that is 483 cells, which must fit inside each rank's slab, so
  it cannot run with any practical regional partition. A distributed EVP needs halo communication every substep (or a
  much smaller halo budget).
- Harness note: NumericalEarth's default free surface on a distributed grid fixes the barotropic substeps from an
  estimated maximum Δt while the serial default adapts them; give both builds an explicit `SplitExplicitFreeSurface`
  (script 05 already does).

### Rank-edge fix (2026-10-15, side session; triton `runs/seaice/icefix`, envs `mpi05_icefix`, `mpi05_icefix2`)

Cause of the ocean differences: `set!(sea_ice.model, h = …, ℵ = …)` fills the interior only, and the coupled model
computes its first interface fluxes when it is built. `_assemble_net_ocean_fluxes!` averages ℵ to the u/v faces
(`ℑxᶠᵃᵃ(ℵ)`, `ℑyᵃᶠᵃ(ℵ)`) for the (1 − ℵ)-weighted wind stress, so the first face of each rank read a zero ℵ halo
(serially, the neighbour's ℵ). Only the first coupled step's stress was wrong (τˣ 5.6e-9 at face 41), but it grew to
T 2.5e-2 °C in 4 steps. After the first ice step, ClimaSeaIce's own `update_state!` fills the halos.
Harness control `FILLICE=true` (fill h, ℵ halos after `set!`) made everything exact before any NumericalEarth change.

NumericalEarth fork branches (from `everything` d9dfcbf6, not PR'd):
- `sea-ice-halos` (8ced8eab), the fix:
  - ef05b7b0 the sea-ice exchanger's `initialize!` calls `update_state!(sea_ice.model)` (fills the ice halos);
  - 8ced8eab `above_freezing_ocean_temperature!` refills the T halos after clipping the initial T to the freezing
    point (control without it, with a supercooled patch: initial τˣ 5.4e-8 on 2 × 1, τʸ 4.9e-8 on 1 × 2).
- `sea-ice-halos-frazil` (200414de) = the above + `update_state!(ocean.model)` after the frazil adjustment, for both
  ClimaSeaIce and `FreezingLimitedOceanTemperature`. The frazil clip changes T after the ocean's `update_state!`;
  the distributed RK step recomputes the diffusivities/pressure next to rank edges after its halo exchange (from the
  clipped T), while the rest of the domain and the serial run keep those of the unclipped T. A halo fill alone is not
  enough (1e-8 left). Cost: a second ocean `update_state!` per coupled step. It also runs the boundary-condition
  updates a second time per step, so check the oblique-radiation stored state before using it in script 05. Enri to decide.
- `distributed-sea-ice-halos` (d4bd4eae) is superseded (halo-only frazil fill); delete it.

Results (4 steps of 300 s, immersed, `distributed_sea_ice.jl` at 94aef96; "exact" = 0.00e+00 in u, v, η, T, S, h, ℵ,
uᵢ, vᵢ and the net ocean fluxes τˣ, τʸ, Jᵀ, Jˢ, Jʷ):

| case | `everything` | `sea-ice-halos` | `sea-ice-halos-frazil` |
|---|---|---|---|
| nodynamics, 2 × 1 | T 2.5e-2, u 8e-4 | exact | exact |
| nodynamics, 1 × 2 | T 1.8e-2, v 5.9e-4 | exact | – |
| nodynamics + supercooled, at build (NSTEPS=0), 2 × 1 / 1 × 2 | – | exact / exact | – |
| nodynamics + supercooled, 2 × 1 | – | T 6.5e-4, u 2.6e-6 (frazil every step) | exact |
| nodynamics + supercooled, 1 × 2 | – | – | exact |
| noice + supercooled, 2 × 1 | – | T 8.7e-9 (run with ef05b7b0 alone) | exact |
| full EVP (8 substeps), 1 × 2 | T 1.1e-2, uᵢ 0.26 | exact (ice velocities too) | exact |
| full EVP, 2 × 1 | T 2.5e-2, uᵢ 0.31 | ocean T 7.0e-6, h 2.9e-6, uᵢ 0.19 of 0.41 at (43, 40) | – |
| full EVP, no ice advection, 2 × 1 | – | ocean T 7.0e-6, uᵢ 1.5e-8, vᵢ 1.2e-7 | – |
| noice, 2 × 1 | exact | exact | – |

EVP notes (cause found, see "EVP 2 × 1 fix" below): with the ℵ halo fixed, the 1 × 2 EVP case is exact. The 2 × 1 one still differs: small ice
velocity differences at the x rank edge from the first step (1-step: uᵢ 1e-8, vᵢ 3e-8, ocean v 2e-6), and, with
WENO ice advection, a large one (0.19 m/s) on the north-wall row next to the rank edge, (43, 40). Likely candidates are
the EVP substep halo handling in x (`only_local_halos` fills with the 2·substeps + 3 halo) and the ice advection at
the wall/rank-edge corner. The 483-cell halo issue above still stands.

### EVP 2 × 1 fix (2026-10-15, side session; triton `runs/seaice/evp2x1`, env `mpi05_icefix4`)

Cause: the ice–ocean drag coefficient of the coupled model's first step. On a distributed grid ClimaSeaIce's
split-explicit solver puts the ice velocities on a halo-extended grid (2·substeps + 3), so `materialize_stress` turns the
`SemiImplicitStress` ocean velocities (`uₑ`, `vₑ` = NumericalEarth's surface `SSU`, `SSV`) into extended *copies*,
refreshed (interior + halos) only inside `time_step_momentum!`. Serially they are the live ocean fields.
NumericalEarth's `_compute_sea_ice_ocean_stress!` evaluates the drag with that same object, `λ = ρ Cᴰ |Δu|` with
`ℑxyᶜᶠᵃ(uₑ)` / `ℑxyᶠᶜᵃ(vₑ)`, so at build time and at the first step (before any ice momentum step) it read the copies'
unfilled halos at the rank edge. τ = −λ uᵢ is still exact there (uᵢ = 0 at start), which is why NSTEPS=0 compared
exact; the ocean uses λ implicitly, so its first step went wrong at i = 40/41 (η 3.8e-7, v 2e-6) and everything else
(ice 1e-8 after one step, uᵢ 0.19 at (43, 40) after 4 with WENO ice advection) followed from that. Found with the
per-substep probe (`EVPPROBE=true`, `restart_test_harness/evp_substep_probe.jl`): the ice fields entered the first
momentum step exact to 1e-16 while η already differed in the interior at the rank edge.

Fix, NumericalEarth fork branch `sea-ice-ocean-drag-halos` (95c5e374, on top of `sea-ice-halos`, not PR'd): one call,
`update_external_stress!(τₛ, grid)` before `_compute_sea_ice_ocean_stress!` in `compute_sea_ice_ocean_fluxes!`
(a no-op serially). It also covers a script that sets the ocean velocities after building the sea ice (the copies'
interiors are taken at construction). It applies to `sea-ice-halos-frazil` as well (cherry-pick).

Results (`immersed`, DT 300, 8 EVP substeps unless noted):

| case | `sea-ice-halos` | `sea-ice-ocean-drag-halos` |
|---|---|---|
| EVP 2 × 1, 1 step | ocean v 2e-6, uᵢ 1e-8 | exact |
| EVP 2 × 1, 4 steps (WENO ice advection) | uᵢ 0.19 at (43, 40) | exact |
| EVP 1 × 2, 4 steps | exact | exact |
| EVP 2 × 1, no ice advection, 4 steps | uᵢ 1.5e-8, vᵢ 1.2e-7 | exact |
| nodynamics 2 × 1, 4 steps | exact | exact |
| EVP 2 × 1, 16 substeps, 8 steps | – | exact |

("exact" = 0.00e+00 in all 14 compared fields; jobs: probe 11808–11812, verification `verify_11814.log` on the fork
branch with no probe. Probe with `EVPFIX=stress` alone, on `sea-ice-halos`, was already exact.)

### Wall/rank-edge corner fix (2026-10-15, same side session; envs `mpi05_icefix5`, `mpi05_icefix6`)

2 × 2 with the centred patch was exact on `sea-ice-ocean-drag-halos`, but a new harness feature `wallice` (the ice
patch centred on the south wall at the x rank edge, pressed against it by the southward wind) was not: 2 × 1 and
2 × 2 both uᵢ 6.2e-3, h 2.9e-4 at (40, 1) after 4 steps (identical, so the 2 × 2 interior corner adds nothing).

Cause (probe, `wall_11817.log`): `fill_halo_regions!` fills a wall's (Bounded side's) halo along the interior of the
other direction only (`fill_halo_size` → `size(grid, loc)`); neither `only_local_halos` nor the corner exchange (no
corner neighbour beyond a wall) fills the corner blocks where a wall halo meets the extended rank-edge halo. The
substeps update u, v in the halo columns, but those corner cells keep the start-of-step values: after substep 1,
u(41, 0) on rank 0 is stale, so at substep 2 σ₁₁ at the wall-halo centre (40, 0) is off by 196 and ζᶠᶠ at the wall node
(40, 1); the wall-row shear stress then moves u(40, 1) (7e-7 at substep 2, uᵢ 9e-4 after one step). With the ice at the
walls marginal (centred patch) the free drift takes over and reads no stresses, which is why it did not show before.

Fix, ClimaSeaIce fork enri66/ClimaSeaIce.jl branch `extended-halo-corners` (81b196b, from v0.5.10, not PR'd):
`time_step_momentum!` fills u, v with boundary-condition kernels rebuilt (`construct_boundary_conditions_kernels`)
with indices spanning the connected sides' halos, so the wall fills also cover those columns (true domain corners
are left as serially). u, v only; h, ℵ, σ and the external fields needed nothing.

Results, all exact (0.00e+00 in the 14 fields), `mpi05_icefix5` = NumericalEarth `sea-ice-ocean-drag-halos` +
ClimaSeaIce `extended-halo-corners` (`csi_11820.log`, `csi_11821.log`):
2 × 1 wallice 1 and 4 steps; 2 × 2 wallice 4 steps; 2 × 1, 1 × 2, 2 × 2 centred patch 4 steps; 2 × 1 wallice with 16
substeps, 8 steps.

Frazil: NumericalEarth branch `sea-ice-halos-frazil-drag` (39047c81) = `sea-ice-halos-frazil` + the drag fix;
`mpi05_icefix6` = that + the ClimaSeaIce branch (`frz_11823.log`, `frz_11824.log`):
EVP + supercooled 2 × 1 and 2 × 2, EVP + supercooled + wallice 2 × 1, nodynamics + supercooled 2 × 2, 4 steps: all
exact. This stack is the full fix set.

Still open:
- The 483-cell halo issue (default 240 substeps) is unchanged: distributed EVP at production settings needs
  per-substep halo exchange in ClimaSeaIce or far fewer substeps.
- `sea-ice-halos-frazil`'s second ocean `update_state!` per step (cost; oblique-radiation stored state unchecked).
- Script 05's GLORYS ice sponge uses global indices (one rank only); it fills the halos after restoring, so it only
  needs porting to local indices for a distributed GSL run.
- Not chased: why 1 × 2 was exact before the drag fix (λʸ at rank 1's first v face reads `uₑ` at j = 0).

### EVP at production substep counts: per-substep halo exchange (2026-10-16; env `mpi05_icefix7`)

ClimaSeaIce fork branch `substep-halo-exchange` (bf463e9, on `extended-halo-corners`, not PR'd): `SplitExplicitSolver(grid;
substeps, extend_halos = false)`, mirroring Oceananigans' `SplitExplicitFreeSurface(; extend_halos)`. The velocities
stay on the grid's own halos (no 2·substeps + 3 extension), the velocity kernels run on the interior, and after every u
and v update `fill_halo_regions!` exchanges the halos and the corner pass of `extended-halo-corners` refills the
walls' halos across the connected halos (needed on 2-D partitions: no corner neighbour beyond a wall). Default
`extend_halos = true` is unchanged. The stress kernels already cover the halos, so nothing else changes.

Harness (`EVP_EXTEND=false`, `xch_11826/7/8.log`), all exact: 2 × 1, 1 × 2, 2 × 2 (centred patch); 2 × 1 and 2 × 2
wallice; 2 × 2 supercooled; **240 substeps 2 × 1 and 2 × 2 wallice** (impossible with extended halos: 483 cells);
the extended mode (default) is still exact (2 × 1 wallice, 2 × 2 supercooled).

Script 05 now runs sea ice on several ranks (ice sponge on local indices; `extend_halos = false` when nranks > 1; needs
`mpi05_icefix7` or equivalent) — `scripts/gsl_mpi.sbatch`, sandbox copy on triton `Models/NESandbox_icempi`, runs in
`runs/seaice/gslmpi`. One GSL day from 2013-12-01, s per simulated hour: 1 rank × 48 threads 40.7; 2 ranks × 24
threads on one node 38.7; 2 nodes × 48 threads 38.2 (no speed-up: likely the 3 × 240 EVP substeps per step, small
latency-bound kernels now each with a halo exchange; not profiled).

**New blocker: the GSL ocean is not layout-exact on its own.** `MAB_STAGE=steps`, 1 vs 2 ranks, `compare_layouts.jl`
(ice fields added to the dumps): with `MAB_SEA_ICE=false` the ocean differs from step 1 — η 1.2e-6 at (139, 120), u
2.9e-5 at the rank edge (139/140) and v 3.8e-5 at (203, 1, 32), on the southern open boundary away from the rank
edge; 4 % of cells after one step. With sea ice the ocean differences are identical and the ice is still exact (no ice
in the first 12 steps on Dec 1), so it is not the sea ice. After one day (`compare_rank_outputs.jl`, sea ice on): ocean
u 1.9e-3, T 2.2e-3 near (141, 9); ice velocities 0.83 m/s (the ice follows the ocean). The MAB ocean was exact; what
differs here (southern boundary across the Bay of Fundy mask, the closed estuary, the largest-basin filter, tides on
444 boundary nodes, …) is the next thing to chase.

Regression test idea (NumericalEarth `test/test_distributed_utils.jl`, which already has `rank_window`): build
`OceanSeaIceModel(ocean_simulation(grid), sea_ice_simulation(grid, ocean; dynamics = nothing); atmosphere)` with a
constant-wind `PrescribedAtmosphere` on a small LatitudeLongitudeGrid, serially and on `Partition(2, 1)` and
`Partition(1, 2)`; `set!(sea_ice.model, ℵ = patch)` with the patch crossing the rank edge; compare
`interior(coupled.interfaces.net_fluxes.ocean.u)` (and `.v`) with the serial window, no time step needed. It fails on
`everything` (first face of the second rank) and passes on `sea-ice-halos`. Same test with an initial T below freezing
across the edge covers the clip refill.

## Suggested plan

1. Pick the region and a winter period with ice growth (check GLORYS `siconc`).
2. Download GLORYS (including the sea-ice variables) and ERA5 for it. The global ERA5 and GLORYS12 archives are on
   antares (see the global-forcing-archives memory note); sync only what is needed to triton.
3. Build the regional grid and get an ocean-only run working with the production settings.
4. Add `sea_ice_simulation(grid, ocean)`, initialise ice from GLORYS, decide how ice is handled at the open boundaries
   (restoring ice thickness/concentration in the sponge, or not).
5. Validate: ice extent and concentration against GLORYS and satellite products (e.g. NSIDC / OSI SAF), ice thickness
   where available, and the ocean's surface temperature and salinity under the ice.

## Open questions

- Which region and period?
- Do ClimaSeaIce and its EVP dynamics run with the distributed (MPI) layout used by script 05? Test early on a small
  2-rank case.
- Open boundaries for ice: is there an ice equivalent of the tracer reservoir, or is sponge restoring enough?
- Rivers: GloFAS runoff should be freshwater with zero salinity entering under ice; check how it interacts with the
  ice model.

## Practical notes

- Sandbox: `~/Dropbox/Models/NumericalEarthSandbox` (Mac) and `/t0/workdir/enrique/Models/NumericalEarthSandbox`
  (triton). Run on triton via antares (`ssh antares`, then `ssh triton`), Slurm; exclude triton24 node[38,41]
  (`sbatch --exclude=node[38,41]`). Environments: `/t0/workdir/enrique/mpi05_ib` (MPI over InfiniBand) and
  `mpi05_rivers`; `JULIA_PKG_GC_AUTO=false`.
- NumericalEarth and Oceananigans come from the fork `everything` branches (enri66), which carry the open-boundary and
  distributed fixes.
- Upstream code (Oceananigans, NumericalEarth, ClimaSeaIce) must follow NumericalEarth's restraint rules; PRs need
  Enri's OK before opening.

## One-rank GSL runs (2026-10-15)

Script 05 now takes `MAB_REGION=gsl` and `MAB_SEA_ICE=true` (one rank only); runs use an isolated sandbox copy on triton,
`/t0/workdir/enrique/Models/NESandbox_ice`, with `scripts/gsl_one_rank.sbatch` (copy of `runs/seaice/gsl/gsl1.sbatch`:
one triton48 node, mpirun -np 1, production settings z*, partial cells, bounded WENO-7 T/S, SSS restoring 0.5 m/day,
T,S,u,v sponge; rivers off so far).

- Grid 278 × 128 × 50 at 1/12°. The bathymetry's largest-basin filter removes the Bay of Fundy and Gulf of Maine (they
  connect to the domain only through the 45°N edge); the estuary closes at 70.9°W; the west side has no wet cells.
- Sea ice: ClimaSeaIce via `sea_ice_simulation`, EVP with 240 substeps, WENO-7 ice advection, no ice at the start.
  ClimaSeaIce applies `forcing` only to the ice velocities, so ice enters the open boundaries by restoring thickness and
  concentration toward GLORYS inside the ocean sponge band after every step (`MAB_ICE_SPONGE_TAU`, default 1 day).
  GLORYS stores ice-free ocean as missing (not 0) in `siconc`/`sithick`, so the script reads its files directly with
  missing as no ice; NumericalEarth's Metadata path (inpainting) would spread ice over open water — this also affects
  `set!` from GLORYS ice and `DatasetRestoring` of ice in NumericalEarth.
- Out-of-bounds read: GLORYS daily series are Cyclical with 2 frames in memory; before the first (noon) frame they wrap
  to the last frame, and the `Interpolated` boundary halo fill at t = 0 (in `set_from_frame!`) reads it outside the
  memory window. With the GSL arrays that segfaulted at random, without a backtrace; `--check-bounds=yes` shows it
  (`NestedModels/interpolated_fts_boundary.jl:112` → Oceananigans `field_time_series_indexing.jl:340`). The MAB runs very
  likely read garbage there too (smaller arrays). Script 05: `MAB_GLORYS_TIME_INDEXING=clamp` (the default for every region
  since 2026-10-15, PRODUCTION_CONFIG.md); it also stops the first 12 hours of boundary data from leaning toward the last day of the record.
- Cost on one triton48 node (48 threads): ~30 s per simulated hour ocean-only, ~41 s with sea ice (~50 h per winter).
- 10-day test (Dec 1–11 2013, `scripts/gsl_ice_vs_glorys.jl`): model ice 18.8k km² on the Labrador side (Lake
  Melville and the coast) vs GLORYS 0.4k; Gulf 1.8k (upper estuary) vs 0.1k; ice up to 3 m in Lake Melville by Dec 6
  (dynamic piling?). GLORYS may miss Lake Melville; needs satellite products (OSI SAF / NSIDC, Canadian Ice Service).
- Full winter running: triton job 11793, tag `ice12_winter` (Dec 1 2013 – Jun 1 2014, 6-hourly ice output).

## GSL ocean layout-exact (2026-10-18, side session; triton `runs/gslx`, env `mpi05_gslx`)

> Status 2026-10-20: the two Oceananigans fixes below were opened as CliMA/Oceananigans.jl#6291 and #6292 (see
> `restart_test_harness/pr_drafts_2026-10-16_sea_ice.md`).

The "new blocker" above is solved: GSL production (z*, partial cells, bounded WENO-7, SSS restoring, T,S,u,v sponge,
tides, open boundaries) is bit-exact between 1 and 2 ranks. It was not GSL-specific: the MAB region with the same
settings also differed (`m1`, rank edge 70/71). The MAB exactness of 2026-10-08 predates z*, partial cells and the
bounded WENO. Three causes, found by bisecting the settings (`g1`–`g7`: only z* and the bounds reproduced) and then
with the reduced harness:

1. **z* pressure gradient (Oceananigans, upstream bug).** On a z* grid `∂xᶠᶜᶜ` is the chain rule
   `∂x|r − ∂x z · ∂z ϕ`, and `∂z pHY′` averaged to the face reads the pressure's top and bottom halos of the
   neighbouring column. At a rank edge that column is in the x-halo: the pressure is computed there (surface kernel
   parameters, recomputed in the buffer pass) but `fill_halo_regions!(pHY′; only_local_halos = true)` fills the
   vertical halos over the interior columns only. Serially that column is interior, so its halos are filled. Effect:
   u tendency at the first/last face of each rank (η 1e-6 after one step). Fix: `update_hydrostatic_pressure!`
   also launches the field's own bottom/top fill over the same columns.
2. **Bounds-preserving WENO limiter (Oceananigans, upstream bug).** `update_state!` starts the tracer halo
   exchange asynchronously and computes the limiters (`update_advection!`) before it completes; their
   reconstructions reach across the rank edge (4.7e-11 after one step). Fix: `synchronize_communication!(tracers)`
   before `update_advection!` (a later synchronization re-copies the same buffers, so it is harmless).
3. **Script 05 `ConsistentNormalFlow` (sandbox).** The test "is level k above the boundary floor" used
   `znode(1, 1, k, grid, …)`, the local cell (1, 1); with z* that node moves with η and is a different column on
   each rank (v 3.8e-5 at the south boundary (203, 1, 32), u at the east face). Now `zc_src[k]`, the resting level
   the table is built from (also a correctness fix on one rank: the threshold no longer follows η at one corner).

Fork branches (enri66/Oceananigans.jl, not PR'd, need Enri's OK):
- `zstar-pressure-vertical-halos` (aa08329d6, from main 5309b1c66): fix 1 + test "Distributed z-star hydrostatic
  pressure" in `test/mpi/hydrostatic_model.jl` (z* box, buoyancy varying in x and y, x4/y4/2×2 vs serial).
- `bounds-limiter-tracer-halos` (f6ddb06d1, from main): fix 2 + test "Distributed bounds-preserving tracer
  advection" (WENO-5 bounded front advected across the rank edges).
- `everything-gsl-exact` (1d906c244) = `everything` (7956c9727) + both, for the runs; triton checkout
  `dev/Oceananigans-gslx`, env `mpi05_gslx` (NumericalEarth `sea-ice-halos-frazil-drag`, ClimaSeaIce
  `substep-halo-exchange`). `everything` itself and the production envs are unchanged.

Tests on triton (`restart_test_harness/zstar_pressure_mpi_test.jl`, `bounds_limiter_mpi_test.jl` run the test
blocks outside Oceananigans' test environment, envs `zt_main`/`zt_fix`/`zt_bnd`): on main 7–9 of 12 z* checks fail
per rank and 2 of 3 bounds partitions fail; with the fixes 12/12 and 3/3 pass.

GSL results (`MAB_STAGE=steps`, `compare_layouts.jl`, 1 rank × 48 threads vs 2 ranks × 24, every dumped field
0.0e+00): production 48 steps (`j0`), production + sea ice 24 steps (`j2`), the stripped configuration 12 steps
(`j1`; on the old env it had a 3.6e-15 difference in one cell at (138, 8, 42), not attributed, gone with the
fixes). Before: η 1.2e-6, u 2.9e-5, v 3.8e-5 after one step (`g0`). One day with sea ice (`MAB_STAGE=run`, `scripts/gsl_mpi.sbatch`
with `JENV=mpi05_gslx SANDBOX=Models/NESandbox_gslx`, tags `x1`/`x2n1` in `runs/seaice/gslmpi`,
`compare_rank_outputs.jl`): every output exact, ocean and ice (ice up to 7 cm, ice velocities to 1.3 m/s); the same
comparison before had ocean u 1.9e-3 and ice velocities 0.83 m/s. The MAB region with the production settings is exact
too on the fixed stack (`m2`, 12 steps).

Reduced reproducer: `restart_test_harness/distributed_features.jl` feature `zstar` (`FEATURES=latlon,teos10,zstar`,
2 ranks, 2 steps: u 2.4e-8 before, exact after; also with immersed, rk3, open, oblique on 2×1, 1×2, 2×2).
`SHOWDIFF` now also compares internal fields (pressure, tendencies, transport velocities, z* fields) with halos;
`zstar_step_probe.jl` (PROBE=1) takes a step one operation at a time on both models, `PROBE_AT=i:j:k` prints the
u tendency's inputs. Script 05's step dumps include the z* grid fields when `MAB_ZSTAR=true`.

Not done: 2-D partitions of GSL; the cost of the extra synchronization (one wait per `update_state!`, after the
diagnostics that overlap it).

## Observed sea ice (2026-10-15, `scripts/sea_ice_obs_vs_model.jl`)

Data on antares `~/Data/sea_ice_obs/<product>` (README.md there: versions, grids, sources, download script), synced to
triton `/t0/workdir/enrique/Data/sea_ice_obs`. All open access (no account, no licence click-through), 2013-11-25 to
2014-06-05:

| product | version | grid | use |
|---|---|---|---|
| OSI SAF OSI-450-a sea-ice concentration CDR (`osisaf_osi450a_v3`) | v3.0 | 25 km EASE2 north, daily | concentration offshore |
| NOAA/NSIDC CDR G02202 (`nsidc_g02202_v6`) | v6 (v06r00) | 25 km NSIDC polar stereographic, daily | concentration offshore |
| Canadian Ice Service regional charts, SIGRID-3, NSIDC G02171 (`cis_g02171/East_Coast`) | v1 | weekly (Mondays) polygons | concentration everywhere incl. Lake Melville, the estuary, fast ice; thickness from stage of development |
| AWI CryoSat-2/SMOS(+S3) merged thickness (`cs2smos_v300`) | v300 | 12.5 km EASE2, 7-day windows moved daily, Oct-Apr | ice volume |

The CIS Eastern Arctic charts (`cis_g02171/Eastern_Arctic`) lie north of the box; no separate SMOS-only thin-ice
product for 2013-14 is openly available (CS2SMOS includes it).

The script regrids everything to the 278 × 128 model grid and compares daily area (Σ ℵ A) and extent (ℵ ≥ 0.15) for the
Gulf, the Labrador side (≥ 51.7°N), Lake Melville (61-58.3°W, 53-54.35°N, part of the Labrador side) and the shelf,
without the 8-cell sponge band, using the masks of `gsl_ice_vs_glorys.jl`. Model = daily mean of `<tag>_sea_ice_rank0.jld2`.
On triton: `RUN_TAG=ice12_winter sbatch -p triton24 --exclude=node[38,41] -J iceobs -o .../iceobs_%j.log --export=ALL
scripts/sea_ice_obs.sbatch` (a few minutes; SANDBOX= a checkout with the script). Output in the run directory:
`<tag>_ice_obs_vs_model.csv` (daily), `<tag>_ice_obs_vs_model_cis.csv` (chart dates, incl. volumes) and
`<tag>_ice_obs_maps.jld2` (all fields on the model grid). Choices:
- Passive microwave has land spill-over: in early December both products show 0.2-0.6 along every Gulf coastline that
  the charts have ice free. Pixels within ~50 km of land are dropped, so those products cover only 28% (OSI) / 59%
  (G02202) of the Gulf's water, 70-87% of the Labrador side and shelf, none of Lake Melville; model and GLORYS are summed
  over the same cells for those columns.
- CIS: each cell takes the polygon its centre falls in; 9% of the model's water cells (coasts, and much of Lake Melville,
  which the 1/12° mask draws wider) fall on chart land and take the nearest chart cell within 2 cells. Thickness per
  stage = WMO mid-range (new 5 cm, grey 12.5, grey-white 22.5, thin FY 50, medium FY 95, thick FY 150 cm), so CIS
  volumes are rough.

Early December 2013 (question from the 10-day test, 10³ km²; CIS on the chart cells, model = ice12_winter, daily mean):

| | Dec 2 | Dec 9 | Dec 16 | Dec 23 | Dec 30 |
|---|---|---|---|---|---|
| Lake Melville: CIS area / GLORYS / model | 1.7 / 0 / 0.2 | 3.2 / 0 / 5.4 | 7.2 / 0 / – | 7.5 / 0 | 8.0 / 0 |
| Labrador side: CIS / GLORYS / model | 2.6 / 0.3 / 0.4 | 5.1 / 0.4 / 11.8 | 26.0 / 4.0 / – | 41.6 / 15.6 | 65.1 / 26.9 |
| Gulf: CIS / GLORYS / model | 0.1 / 0 / 0 | 3.7 / 0 / 0.7 | 30.6 / 10.2 / – | 52.2 / 13.1 | 78.1 / 15.9 |

- GLORYS has no ice in Lake Melville at any time in the winter (CIS: ice from Dec 2, the whole lake by mid-December,
  grey ice 10-15 cm early on, still covered in late May), and far less than CIS on the Labrador coast and in the Gulf
  through December: GLORYS misses coastal, fast and new ice that the charts map. OSI-450-a and CS2SMOS (unfiltered, see
  the maps) also show 0.6-0.7 in western Lake Melville on Dec 9. So the model's early Lake Melville / Labrador-coast ice
  is real, not a model artefact; the model is ahead of CIS (Dec 9: Melville 5.4 vs 3.2, Labrador 11.8 vs 5.1, all
  of Lake Melville at ℵ ≈ 1 where CIS has 0.7-1 in the west and < 1/10 in the east) and too thick (Melville volume 1.5 km³ vs CIS 0.3, i.e. ~28 cm vs
  grey ice ~10 cm). The ice sponge restoring toward GLORYS likely removes coastal ice near the north edge in
  December.
- Offshore (on the passive-microwave cells) GLORYS agrees well with OSI-450-a and G02202 from late December to April
  (e.g. Labrador Feb 10: GLORYS 75.6 vs OSI 69.3; shelf Mar 10: 114.8 vs 112.1), as expected if it assimilates them;
  on the G02202 cells GLORYS has less Gulf ice in February-March (Mar 10: 114 vs 165).
- Over the CIS cells GLORYS is below the charts all winter (Gulf Feb 10: 134 vs 252; Labrador Mar 17: 99 vs 131),
  and in spring the charts keep much more ice (May 26 Labrador 51.7 vs 13.7; volume 60 vs 10.5 km³ with the stage
  thicknesses). CS2SMOS volume is lower than GLORYS mid-winter on its cells (Gulf Feb 24: 64 vs 83 km³).
- The model columns cover the run so far (ice12_winter had reached ~Dec 13 at the time); rerun the sbatch when job
  11793 finishes.

## NumericalEarth fixes for GLORYS ice and the out-of-bounds read (2026-10-16, fork branches, no PRs yet)

> Status 2026-10-20: opened as NumericalEarth#790 (GLORYS ice) and #791 (clockless Interpolated read), and the Oceananigans
> clocked-halo-fill variant as CliMA/Oceananigans.jl#6293; see the PR drafts file.

Both on enri66/NumericalEarth.jl, branched from upstream `main` (526d401c), each with a test that fails on `main`.
- `glorys-sea-ice-open-water` (7870435d): GLORYS `retrieve_data` reads missing `siconc`/`sithick` as 0 (ice-free ocean
  and land). On `main`, 2013-12-01 GSL box onto a 100×50 grid: plain `Field`/`FieldTimeSeries` NaN in 4674 of 5000 cells,
  with `NearestNeighborInpainting(Inf)` (DatasetRestoring's default) ice in 4996 cells; with the fix no NaN and 346 iced
  cells both ways. Test in `test/test_glorys_surface_variables.jl` (synthetic file, offline). Old `siconc_*_inpainted.jld2`
  / `sithick_*_inpainted.jld2` caches hold the spread ice: delete them. ECCO: Field/FTS are fine (no inpainting for sea
  ice, zeros kept), but ECCO4Monthly's `default_mask_value` is 0, so DatasetRestoring of SIarea/SIheff (inpainting
  Inf) would mask all ice-free cells and spread ice the same way (from the code; no ECCO ice files at hand to run).
- `cyclic-fts-boundary-window` (4482e8f5): an `Interpolated` FTS source filled without a clock is read at the first
  snapshot in memory instead of at t = 0. Cause of the GSL crash (job 11791 log): the consistency loop leaves u's window
  at (10, 11), `fts_u[1]` moves it to (1, 2), and the clockless fill in `set_from_frame!` at t = 0 needs (11, 1). Test in
  `test/test_nested_simulation.jl`: BoundsError with `--check-bounds=yes` on `main`, wrong halo value without.
- Not fixed (Oceananigans, follow-up): a *clocked* halo fill also reads outside the window when host reads (`fts[n]`)
  moved it after the last `update_state!` — plain Oceananigans reproduces it (`ValueBoundaryCondition(fts)`, `fts[2]`,
  `set!(model, T = 0)` → BoundsError). `prepare_halo_fill!` cannot fix it: on Oceananigans `main` the Field halo fill
  adapts the boundary conditions (`convert_to_device`) before the hook runs, so the window must be updated on the host
  series in `fill_halo_regions!(::Field, clock, …)` before that conversion. Script 05 is safe with Clamp as it is.
