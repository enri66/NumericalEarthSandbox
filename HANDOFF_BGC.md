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
- Prototype `runs/bgc/bgc_feb19` (job 11794; 1/12°, 4 triton24 ranks, 2019-02-01 → 04-01, WOA + carbonate
  placeholders): ran to the end in ~15 h (~36 s per model hour after a 55-min first step). Bloom: P up to 4.0 mmol N/m³
  over Nantucket Shoals (70.5°W 41.5°N, uniform over the shallow column), Z up to 2.0, surface O₂ up to 392. Open
  problem: isolated bottom-cell hot spots, e.g. the bottom cell at 67.12°W 40.63°N has NO₃ 34.7 and O₂ 71 against 9 and
  218 in the cell above (9 of 465 222 wet cells with NO₃ > 25, 3 with O₂ < 150; NH₄ max 1.8, also a bottom cell).
  NO₃ dips to −0.01 in places.
- Bottom-cell hot spots, diagnosed 2026-10-16 (scripts in triton `runs/bgc/`: hotspots_feb19.jl, hotspots_vs_dz.jl,
  hotspot_budget.jl, hotspot_ts.jl; bottom-cell heights from the `<TAG>_bottom_rank<r>.jld2` files mab_bgc.jl now
  writes): 60 bottom cells end with NO₃ > 5 mmol/m³ above the cell over them, mostly on the shelf break and upper slope.
  - Not a rank-edge error: they sit at every column position within the ranks.
  - Not thin partial cells: hot and non-hot bottom cells have the same height ÷ level distribution (median 0.33 vs 0.41;
    0.50 % of the cells thinner than 0.3 of their level are hot, 0.47 % of the thicker ones).
  - Not the sediment or sinking detritus: in the worst column (105, 78; 226 m deep, 5.9 m bottom cell) the column holds
    2–13 mmol N/m² of detritus and the bottom-cell NH₄ stays ~0.3, while the bottom cell gains ~190 mmol N/m² as NO₃
    between days 39 and 48; the column inventory rises by ~750 mmol N/m² (lateral inflow of slope water).
  - Values outside every source: NO₃ 37 and O₂ 44 (initial and boundary data: NO₃ ≤ 25, O₂ ≥ ~130), so advection
    creates new extrema there. The same bottom cell is 2–3 °C colder than the cell 6 m above while the vertical velocity
    at its top face is −5 to −11 m/day. This is the bottom-cell advection problem of PRODUCTION_CONFIG.md, visible
    because the BGC tracers are only bounded to [0, 1e5]; T and S stay inside their loose physical bounds.
  - Data-range bounds (restart `bgc_feb19_bnd`, job 11847: day 35 checkpoint of bgc_feb19 to day 50 with
    MAB_BGC_BOUNDS="NO₃:0:30,O₂:50:450"; triton `runs/bgc/compare_bnd.jl`): on the shared days 39–47 the NO₃ maximum stays
    ≤ 24 (37 unbounded) and the O₂ minimum ≥ 84 (59), with domain totals within 4e-5 of the unbounded run. But the hot
    cells remain (39 vs 43 on day 47) and the worst bottom cell still drifts (NO₃ 5 → 23, O₂ 247 → 84): the bounds keep
    the values possible, they do not stop the bottom cell from decoupling from the water above.
  - Full Copernicus runs (59 days from 2019-02-01, 4 ranks, ~15 h each; `bgc_feb19_cm` unbounded, job 11842, and
    `bgc_feb19_cmb`, job 11917, with MAB_BGC_BOUNDS="NO₃:0:30,O₂:50:450,DIC:0:2300"; `runs/bgc/compare_runs.jl`): the
    unbounded run reaches NO₃ 47, O₂ 36 and DIC 2413 in bottom cells; the bounded one stays within NO₃ ≤ 24.3,
    O₂ ≥ 82, DIC ≤ 2162, and its domain sums of NO₃, O₂ and DIC differ by at most 1.3e-4. The number of hot bottom cells
    (13–35 over days 8–53) is almost the same in both (bounded: 2–4 fewer). So the bounds are a safeguard; the cause
    in the bottom-cell advection next to the immersed slope remains open.
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

## Validation of `bgc_feb19_cmb` (2019-02-04 → 03-29 daily frames, 2026-10-20)

`bgc/validate_bgc.jl` (run on triton with the `mpi05_bgc` environment; data from `bgc/download_validation_data.sh`):
satellite chlorophyll (Copernicus L3 multi-sensor, 4 km, only where the satellite saw the surface; 32 % of pixels) and
the OOI Pioneer sensors at CNSM and ISSM (7 m near-surface frame, seafloor node), with the model's chlorophyll taken as
1.31 mg Chl/mmol N × P. Figures and `summary.txt` in `runs/analysis/bgc_bgc_feb19_cmb/`.

- **Surface chlorophyll (maps, matched cell-days).** Whole-domain log10 correlation 0.64, but the shelf-to-ocean contrast is
  too weak. Satellite / model means (mg/m³): shelf (< 100 m) 2.4 / 1.0, upper slope 0.70 / 1.0, deep (> 1000 m) 0.31 / 0.57.
  Region-mean time series: r = 0.05 (shelf), 0.52 (slope), 0.78 (deep; mostly the shared late-winter rise).
  The offshore excess is there on the first day (about 2×), so it comes from the initial and boundary data (PISCES
  phytoplankton carbon × 16/106 × 1.31), not from the model's own dynamics. Near the coast the satellite may itself be too high (or the
  fluorometer too low): at ISSM the satellite is 1.8× the in-situ fluorometer over 28 matched days, at CNSM they agree
  (0.82 vs 0.88), so the shelf bias of the model is somewhere between 1.3× and 2.4× too low.
- **Chlorophyll at the moorings.** CNSM / ISSM mean fluorometer 0.96 / 0.62 against model 1.12 / 0.95; the daily series are
  anticorrelated (r = −0.60, −0.17): the model bloom comes late and large at CNSM (2.3 mg/m³ at days 37–42, the fluorometer
  peaks at 1.5 on day 18), and rises at ISSM where the fluorometer falls.
- **Nitrate (SUNA).** Far too low and lost too early. CNSM: sensor 4.5–8, model 1.8 on day 3 (data too low from the start)
  falling to 0 after day 40; ISSM: sensor 0.3–2.6, model 1.4 falling to 0 by day 38. The model's nitrate runs out as its
  bloom starts, while the sensors still read nitrate. River water carries NO₃ = 0 (the default `freshwater_tracer_content`),
  so it can only dilute the shelf nitrate; river loads are still to be added.
- **Oxygen.** The model starts 83–87 % saturated at 7 m (its initial PISCES O₂) and air-sea exchange takes about 35 days to bring it to ~100 %, so the first
  month of O₂ is a spin-up (not tested why it is so slow: mixed-layer depth and the gas transfer velocity are the
  candidates). The Pioneer
  near-surface optodes read 123–131 % saturation in Feb–Mar with QC flag 1, which is not physical, so they cannot check the
  absolute value (mean offset −92/−104 mmol/m³); their variability is small. The seafloor optode at ISSM (11 days)
  reads 270 against the model's 230 mmol/m³.
- **Not done:** vertical sections (the wire-following profilers have fluorometers), discrete samples or WOA for an independent
  oxygen and nutrient reference, the Copernicus hindcast's own chlorophyll and nitrate against the same observations (to
  separate its bias from the model's), pCO₂.

## Practical notes

- Sandbox: `~/Dropbox/Models/NumericalEarthSandbox` (Mac) and `/t0/workdir/enrique/Models/NumericalEarthSandbox`
  (triton). Run on triton via antares (`ssh antares`, then `ssh triton`), Slurm; exclude triton24 node[38,41]
  (`sbatch --exclude=node[38,41]`). Environments: `/t0/workdir/enrique/mpi05_ib` (MPI over InfiniBand) and
  `mpi05_rivers`; `JULIA_PKG_GC_AUTO=false`. OceanBioME must be added to the environment used.
- NumericalEarth and Oceananigans come from the fork `everything` branches (enri66), which carry the open-boundary and
  distributed fixes.
- Upstream code (Oceananigans, NumericalEarth, OceanBioME) must follow NumericalEarth's restraint rules; PRs need
  Enri's OK before opening.
