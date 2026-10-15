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
  likely read garbage there too (smaller arrays). Script 05: `MAB_GLORYS_TIME_INDEXING=clamp` (default for gsl, the MAB
  default unchanged); it also stops the first 12 hours of boundary data from leaning toward the last day of the record.
- Cost on one triton48 node (48 threads): ~30 s per simulated hour ocean-only, ~41 s with sea ice (~50 h per winter).
- 10-day test (Dec 1–11 2013, `scripts/gsl_ice_vs_glorys.jl`): model ice 18.8k km² on the Labrador side (Lake
  Melville and the coast) vs GLORYS 0.4k; Gulf 1.8k (upper estuary) vs 0.1k; ice up to 3 m in Lake Melville by Dec 6
  (dynamic piling?). GLORYS may miss Lake Melville; needs satellite products (OSI SAF / NSIDC, Canadian Ice Service).
- Full winter running: triton job 11793, tag `ice12_winter` (Dec 1 2013 – Jun 1 2014, 6-hourly ice output).

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
