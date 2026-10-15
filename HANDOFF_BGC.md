# Handoff: biogeochemistry in the MAB setup (session "NE-BGC")

Written 2026-10-14 from the main MAB session. Goal: add biogeochemistry (BGC) to the production Mid-Atlantic Bight
simulation, with realistic initial, boundary and river inputs, and validate it against observations.

## Where things stand

- The regional ocean setup is `scripts/05_mab_glorys_tides_reservoirs_distributed.jl` (MPI, per-rank output), run on
  triton. Agreed production settings and open items: `PRODUCTION_CONFIG.md` (z*, partial cells, fixed open boundaries,
  bounds-preserving WENO-7 for T/S, GloFAS rivers, SSS restoring with a river mask).
- Open boundaries: oblique radiation for baroclinic velocities, Flather-type barotropic conditions, tracer reservoirs,
  GLORYS T/S/u/v sponge, tides. These are working.
- Rivers: GloFAS daily discharge (global v4 archive on antares) routed to coastal cells, with hand-placed Hudson and
  Delaware mouths, entering as a surface freshwater flux (`MAB_RIVERS=true`, `mpi05_rivers` environment).
- The analysis package (`scripts/analysis_package.sh`) compares with GLORYS, Argo and the OOI Pioneer moorings.

## What NumericalEarth provides

- BGC models come from OceanBioME (a weak dependency; NumericalEarth's `NumericalEarthOceanBioMEExt` loads when
  `using OceanBioME`): the nutrient-phytoplankton-zooplankton-detritus family (`NutrientsPlanktonDetritus`, including
  LOBSTER-type configurations), optional `Oxygen()` and a carbonate system (`CarbonateSystem()`: DIC and alkalinity).
  Passed to `ocean_simulation(grid; biogeochemistry, ...)`.
- Coupling in the extension (`ext/NumericalEarthOceanBioMEExt/`):
  - `PARFromShortwave(grid)`: surface photosynthetically active radiation from the atmosphere's shortwave;
  - air-sea gas exchange of O₂ and CO₂ (atmospheric `pCO₂` as an atmosphere tracer, pressure from the atmosphere),
    reduced by the sea-ice concentration.
- River tracer loads: `ocean_simulation(...; freshwater_tracer_content = (; Alk = field, ...))` gives the tracer
  carried by the freshwater flux (the coupled test passes river alkalinity this way). This can put nutrients, DIC and
  alkalinity on the GloFAS rivers.
- Data:
  - WOA: nitrate, phosphate, silicate, dissolved oxygen (climatology);
  - ECCO-Darwin: DIC, alkalinity, phosphate, nitrate, silicate, oxygen, iron, DOP, POP;
  - SeaWiFS: chlorophyll;
  - the restoring machinery accepts the BGC tracers (`src/DataWrangling/restoring.jl`), so they can be nudged in the
    same boundary sponge as T and S.
- Reference use: NumericalEarth test `test/test_ocean_sea_ice_oceanbiome.jl` (ocean + sea ice + BGC with PAR, oxygen,
  carbonate system, prescribed pCO₂ and a river alkalinity source).

## Gaps to consider

- No reader yet for the Copernicus Marine global BGC reanalysis (the biogeochemical counterpart of GLORYS). It would
  give boundary conditions consistent with GLORYS physics and daily variability; adding it would follow the GLORYS
  module (`src/DataWrangling/GLORYS/`) and the same `copernicusmarine` download path.
- River nutrient and carbon concentrations: GloFAS gives discharge only. Concentrations would come from a river
  database or literature values for the Hudson, Delaware, and the smaller rivers.
- Light attenuation choice (fixed versus chlorophyll-dependent) and sinking of detritus in the shallow, tidally mixed
  shelf.

## Suggested plan

1. Choose the model with Enri: start simple (nutrients-phytoplankton-zooplankton-detritus with oxygen), add the
   carbonate system later.
2. Prototype at 1/12° on a short run (the production config at 1/12° runs in a few hours on 2 triton24 nodes):
   initialise from WOA or ECCO-Darwin, PAR from ERA5 shortwave, BGC tracers restored in the boundary sponge.
3. Add river loads through `freshwater_tracer_content`.
4. Validate: surface chlorophyll against satellite ocean colour, nutrients and oxygen against WOA and ECCO-Darwin,
   the OOI Pioneer moorings where they measured BGC in 2019, and the seasonal shelf cycle (spring bloom, summer
   bottom-oxygen minimum on the shelf).
5. Then move to 1/24° or 1/60°.

## Open questions

- Which BGC model, and which boundary data (WOA climatology, ECCO-Darwin, or a new Copernicus BGC reader)?
- Which Pioneer BGC measurements exist for the run period (fluorometers, nitrate, oxygen)?
- Cost: every BGC tracer adds a WENO-advected 3D tracer; check the step cost at 1/12° before going finer.
- Distributed runs: check that the BGC extension's kernels (gas exchange, PAR) work with the MPI layout of script 05
  on a small 2-rank case.

## Decisions and status (NE-BGC session, 2026-10-14)

- Model: LOBSTER + `Oxygen()` + `CarbonateSystem()` from the start (NO₃, NH₄, P, Z, DOM, sPOM, bPOM, DIC, Alk, O₂), two-band
  light under `PARFromShortwave`, and `InstantRemineralisationSediment` returning sinking sPOM/bPOM to NH₄ at the seafloor.
  Without the sediment the sinking detritus leaves through the immersed bottom and its nitrogen is lost (column test).
- Data: WOA first (NO₃ and O₂: monthly in the top 200 m, annual below 300 m; the 1° monthly nutrient analyses are
  noisy at depth), DIC and Alk as salinity placeholders; then a reader for the Copernicus global BGC reanalysis
  (GLOBAL_MULTIYEAR_BGC_001_029, 1/4°, daily, 1993 on) for the initial state, boundaries and sponge.
- Prototype period: from 2019-02-01 (winter-spring bloom). Triton's MAB cache holds 2019-04 to 2020-04 only; Feb–Mar
  2019 is in the antares global archives (ERA5 on A2, GLORYS12 and GloFAS on A3) and needs converting into the
  per-region cache files.
- Code: `bgc/` holds the BGC run (`mab_bgc.jl`, a copy of script 05 at 45d50d8 with sections marked "BGC"),
  `bgc_data.jl` (BGC series on the whole-domain source grid) and `bgc_run.sbatch` (production physics). Physics
  changes to script 05 must be carried over by hand.
- Triton environment: `/t0/workdir/enrique/mpi05_bgc` (mpi05_rivers + OceanBioME 0.19.4 + WorldOceanAtlasTools).
- 2-rank test (job 11727, `runs/bgc/bgc2r`, 1/12°, 3 days from 2019-09-01, production physics with rivers): ran to the
  end, all tracers bounded, phytoplankton growing (P max 0.10 → 0.31 mmol N/m³), no discontinuity across the rank
  boundary at t = 0 or day 3. About 50 s per model hour on 2 triton24 nodes (script 05 alone: ~18 s); the first time step
  compiles for ~57 min. OceanBioME warns that the sediment is stepped outside the split Runge-Kutta scheme (first order,
  not mass conserving): check the nitrogen budget before trusting bottom fluxes.
- The output files hold the halos (7 cells): strip them before joining ranks.
- Feb–Mar 2019 forcing for the MAB is now on triton (built from the antares archives by the converter session).
- NumericalEarth fix: the OceanBioME extension's gas exchange read `sea_ice.model.ice_concentration`, which fails with
  `OceanOnlyModel` (the sea ice is `FreezingLimitedOceanTemperature`); it now calls `sea_ice_concentration`. On the
  fork's `everything` branch as d9dfcbf6; not proposed upstream yet.

## Copernicus BGC reader (2026-10-15)

- Product: GLOBAL_MULTIYEAR_BGC_001_029 (PISCES, 1/4°, 75 levels 0.5–5902 m, version 202406, 1993-01 → 2026-08).
  Datasets: `cmems_mod_glo_bgc_my_0.25deg_P1D-m` (chl, no3, po4, si, o2, nppv), `..._P1M-m` (the same plus fe, phyc,
  ph, spco2) and `..._static` (deptho and mask in part "mask"). **No DIC or alkalinity in this product**: PISCES
  dissic/talk exist only in the analysis-forecast product GLOBAL_ANALYSISFORECAST_BGC_001_028, from 2021-11.
- NumericalEarth: `CopernicusBGCStatic`, `CopernicusBGCDaily`, `CopernicusBGCMonthly`, subtypes of `GLORYSDataset` in the
  GLORYS module (`src/DataWrangling/GLORYS/copernicus_bgc.jl`), downloaded through the CopernicusMarine extension.
  Names: `:chlorophyll`, `:nitrate`, `:phosphate`, `:dissolved_silicate`, `:dissolved_oxygen`, `:net_primary_production`
  (daily and monthly); `:dissolved_iron`, `:phytoplankton_carbon`, `:ph`, `:carbon_dioxide_partial_pressure` (monthly);
  `:depth` (static). Units are the product's, with no conversion: mmol/m³ for the concentrations (LOBSTER's units),
  chl mg/m³, nppv mg C/m³/day, ph total scale (dimensionless, 3-D on all 75 levels), spco2 Pa (surface only). The
  native grid is centred exactly on the file's 1/4° points. Fork branch `copernicus-bgc` (fe421393, on d9dfcbf6 =
  fork `everything`), local clone worktree `~/dev/NumericalEarth-copernicus-bgc`, **not pushed**; test
  `test/test_copernicus_bgc_downloading.jl` passes on the Mac (nitracline, O₂ and pH ranges at 35–37°N 200–202°E, static
  depth, 2-D pCO₂). Restraint checklist run; the diff is +132/−11 including docs and the test.
- Sandbox: `MAB_BGC_DATA=copernicus` (bgc_data.jl) takes NO₃, O₂ and P from `CopernicusBGCMonthly` at the same
  mid-month frames as WOA (P = phyc × 16/106, mmol C → mmol N), read for the GLORYS box ±1° and inpainted
  (NearestNeighborInpainting(Inf), no cache); P then has a series, boundary reservoir and sponge, and leaves
  BGC_UNIFORM. DIC/Alk stay salinity placeholders (Enri chose DIC from PISCES pH with Alk from salinity, to be written
  in the NE-BGC session). `bgc_series()` unchanged; WOA stays the default. `bgc/bgc_copernicus_download.jl` fetches
  NO₃, O₂, phyc and pH for a span of months; bgc_run.sbatch takes `BGC_DATA`, `JENV` and `SANDBOX`.
- Triton: MAB subset 2018-12 → 2020-01 (NO₃, O₂, phyc, pH; 56 files, 46 MB) in `~/Data/NumericalEarth/copernicus_bgc`.
  Env `/t0/workdir/enrique/mpi05_bgc_cmems` = mpi05_bgc with NumericalEarth from `/t0/workdir/enrique/dev/NumericalEarth-cmemsbgc`
  (a copy of the branch). Separate sandbox copy `/t0/workdir/enrique/Models/NESandbox_bgc_cmems` (the main one's bgc/
  was in use by job 11794). The copernicusmarine wrapper in `julia_depot/scratchspaces/cd43e856-.../copernicusmarine-v2.4.1`
  had been replaced by the bundled binary (h5py error again) and was restored.
- 2-rank build (job 11802, 1/12°, MAB_STAGE=model, start 2019-02-01, 30 min): frames 2019-01-15 → 2019-03-15;
  initial NO₃ 0.03–22.0, O₂ 131–366, P 0.002–2.74 mmol/m³ (WOA-based NO₃/O₂ and uniform P = 0.1 before); BGC sponge on
  NO₃, O₂, P, Alk, DIC.
- DIC from pH (Enri's choice, 2026-10-15): with MAB_BGC_DATA=copernicus, Alk = 2300 S / 35 and DIC is the value whose
  total-scale pH at the GLORYS T, S and the cell's pressure matches the hindcast's monthly pH (`dic_from_pH`,
  `carbonate_from_pH_series` in bgc_data.jl). Build test job 11805 (2 ranks, 2019-02-01, env mpi05_bgc_cmems, sandbox
  copy `Models/NESandbox_bgc_ph`): DIC 233–2167 mmol/m³ on the first day, domain DIC/Alk 0.895, no cell without a
  solution. DIC inherits any error of the salinity alkalinity almost one for one.
- The reader (NumericalEarth fe421393) is on the fork's `everything` (pushed 2026-10-15 with Enri's OK) and
  `mpi05_bgc` is updated to it (old manifest kept as `Manifest.toml.bak-2026-10-15`), so MAB_BGC_DATA=copernicus
  runs with the default environment. Not proposed upstream.

## Practical notes

- Sandbox: `~/Dropbox/Models/NumericalEarthSandbox` (Mac) and `/t0/workdir/enrique/Models/NumericalEarthSandbox`
  (triton). Run on triton via antares (`ssh antares`, then `ssh triton`), Slurm; exclude triton24 node[38,41]
  (`sbatch --exclude=node[38,41]`). Environments: `/t0/workdir/enrique/mpi05_ib` (MPI over InfiniBand) and
  `mpi05_rivers`; `JULIA_PKG_GC_AUTO=false`. OceanBioME must be added to the environment used.
- NumericalEarth and Oceananigans come from the fork `everything` branches (enri66), which carry the open-boundary and
  distributed fixes.
- Upstream code (Oceananigans, NumericalEarth, OceanBioME) must follow NumericalEarth's restraint rules; PRs need
  Enri's OK before opening.
