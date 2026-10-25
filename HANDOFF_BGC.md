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

## River loads (2026-10-20)

`MAB_RIVER_LOADS` (default on with rivers) in `bgc/mab_bgc.jl`: the river water carries NO₃ and O₂ through
`freshwater_tracer_content`. That field is the flux Σ cᵢ Jʷᵢ (mmol/m²/s), not a concentration, so a callback rescales a
static concentration map by the runoff flux (`model.interfaces.exchanger.land.state.runoff_freshwater_flux`/ρ) every step.
A river cell takes the concentration of the nearest hand-placed mouth (Hudson, Delaware) within 100 km, "other" elsewhere.
- NO₃ (mmol/m³, `MAB_RIVER_NO3`): Hudson 50 (published as "usually > 50 µM", Roland and Cole 1999, via BioNumbers: a
  secondary summary, lower bound), Delaware 50 and other rivers 20 (both ASSUMED: no values found). Replace them with
  measured values (USGS NWIS nitrate+nitrite at Trenton and the Hudson gauges; Boyer et al. 2002 for the rest).
- O₂: fresh-water saturation at 5 °C (399 mmol/m³, `MAB_RIVER_O2_TEMPERATURE`). NH₄, DOM, particles, DIC and Alk: zero
  (zero Alk and DIC matches the salinity-derived placeholders at S = 0); river alkalinity and DIC would need values.
- Check at the start of `bgc_feb19_riv` (job 11941): Q = 1875 m³/s, NO₃ load 91.8 mol/s (mean 49.0 mmol/m³); the Hudson
  (1013 m³/s) and Delaware (775 m³/s) are 95 % of the discharge. Only 6 mouths exceed 20 m³/s and the others total ~87 m³/s:
  the Connecticut River does not seem to be among the routed rivers in Feb 2019 (not investigated).
- NO₃ advection bound raised to 60 for this run (river water at 50 would otherwise sit above the limiter's bound).
- **RC4USCoast** (Gomez et al. 2023, ESSD 15, 2223, doi:10.5194/essd-15-2223-2023; data NOAA NCEI accession 0260455,
  doi:10.25921/9jfw-ph50, v3.3): monthly series (1950–2022) and monthly climatologies of USGS river chemistry for 140 rivers
  (52 on the East Coast), netCDF, in the model's units (nutrients, DIC, O₂ in mmol/m³, Alk in meq/m³): alk, dic (calculated
  from Alk, pH and T; measured `dicm` is sparse), no3, nox, nh4, don, po4, do, doc, sio2, pH, temperature, discharge. Files on
  triton: `/t0/workdir/enrique/Data/RC4USCoast/` (mclim_19902022_chem.nc, mclim_19902022_disc.nc, series_*.nc,
  dataset_info.xlsx; 21 MB, fetched from `ncei.noaa.gov/data/oceans/archive/arc0207/0260455/3.3/data/0-data/`).
  Limits: stations are near river mouths but mostly at the tidal limit (the Delaware is at Trenton, so estuarine processing
  before the shelf is not included); gaps are not filled, and a month with too few records falls back to the annual mean
  (the Raritan is constant all year); **no Hudson station**.
  Feb–Mar means of the 1990–2022 climatology (mmol/m³; Alk meq/m³; Q m³/s, `runs/analysis/rc4_table.jl`):

  | river (mouth) | NO₃ | NH₄ | DON | PO₄ | Alk | DIC | O₂ | Q |
  |---|---|---|---|---|---|---|---|---|
  | Delaware (−75.56, 39.47) | 64.4 | 2.2 | 12.1 | 0.58 | 703 | 746 | 433 | 459 |
  | Connecticut (−72.35, 41.29) | 28.8 | 5.7 | 9.8 | 0.52 | 463 | 511 | 432 | 579 |
  | Housatonic (−73.11, 41.17) | 35.3 | 2.9 | 12.2 | 0.40 | 1683 | 1736 | 428 | 106 |
  | Raritan (−74.28, 40.49) | 84.9 | 3.3 | 17.7 | – | 1058 | 1113 | 340 | 33 |
  | Pawcatuck (−71.85, 41.32) | 26.0 | 2.4 | 15.3 | 0.36 | 186 | 251 | 407 | 26 |
  | Susquehanna (−76.09, 39.57; Chesapeake, outside the box) | 101.2 | 4.4 | 13.0 | 0.31 | 860 | 905 | 414 | 1682 |

  Against the values used in `bgc_feb19_riv`: the Delaware (50 assumed) is 64; the Connecticut and Housatonic (20 assumed
  for "other") are 29 and 35; the measured O₂ of the rivers (428–433) is above the 399 of fresh water at 5 °C (about 422 at
  3 °C). Alk and DIC are 460–1700 and are not loaded at all yet; NH₄ (2–6) and DON (10–18) are small but not zero.
  **Implemented (2026-10-20, `bgc/river_chemistry.jl`, `MAB_RIVER_CHEMISTRY=rc4uscoast`, the default):** the river cells
  of each class take the discharge-weighted monthly climatology (1990-2022) of NO₃, NH₄, DOM (= don), Alk, DIC and O₂ of the
  RC4USCoast rivers they stand for, linear between mid-months and multiplied by the runoff flux every step:
  Delaware = Delaware + Schuylkill (the GloFAS cell lies below their confluence), Hudson = PROXY (Connecticut, Housatonic,
  Delaware, Schuylkill: the dataset has no Hudson station), other = the other rivers with a mouth in the box. February
  values (NO₃ / NH₄ / DOM / Alk / DIC / O₂): Delaware 102 / 6.2 / 14.6 / 852 / 899 / 448, Hudson proxy 66 / 6.2 / 12 / 764 /
  810 / 442, other 38 / 6.0 / 11 / 680 / 729 / 425. `MAB_RIVER_CHEMISTRY=assumed` keeps the constant NO₃ and O₂ of
  `bgc_feb19_riv`. The run `bgc_feb19_rc4` (job 11946) replaced `bgc_feb19_riv` (11941, cancelled after 52 min) with
  BGC_BOUNDS="NO₃:0:150,O₂:50:500,DIC:0:2300" (the river values would sit above the old limits).
  Still open: the Hudson has no data (50 is the published ">50 µM"; its Alk and DIC need a proxy such as the median of the
  neighbouring rivers); **the Connecticut (clim. 579 m³/s in Feb–Mar) does not seem to be routed in the model** (the
  "other" rivers total 87 m³/s on 2019-02-01); the Chesapeake outflow (Susquehanna and Potomac, ~2000 m³/s) is outside the box.

## Analysis of `bgc_feb19_rc4` (RC4USCoast river loads; 59 days, finished cleanly; 2026-10-21)

Compared with `bgc_feb19_cmb` (same setup, no river loads, NO₃ bound 30): `bgc/compare_rivers.jl` (figures and summary in
`runs/analysis/rivers_bgc_feb19_rc4_vs_bgc_feb19_cmb/`), `bgc/validate_bgc.jl` for the observations, `runs/analysis/negatives.jl`.
- **The Pioneer sensors cannot test the loads.** Every CNSM/ISSM statistic of the rc4 run equals the cmb run's to the digit:
  the moorings are 200–270 km east of the Hudson and Delaware mouths, upstream of the southwestward shelf flow. (The
  satellite shelf statistics do change slightly: bias −0.25 → −0.24, r +0.20 → +0.25.)
- **The loads act in the near field only.** Surface NO₃ rises by more than 1 mmol/m³ in 52 cells of the mean map (19 above 5,
  none above 20), at most 23 km from the Hudson mouth and 45 km from the Delaware (median depth 10 m). Mean over the
  100-km regions (days ≥ 8, without / with loads): Hudson NO₃ 0.04 / 0.22, chl 1.00 / 1.09; Delaware NO₃ 0.70 / 2.23, chl 1.24 / 1.39 mg/m³.
  Without loads the plume regions' surface NO₃ is gone by day 20–35; with loads it lasts 2–4 weeks longer (Delaware still
  0.6 on day 56), and chlorophyll is higher by 0.1–0.4 mg/m³ from about day 30.
- **Satellite chlorophyll in the plume regions** (satellite 4.4–4.6 mg/m³, model 0.9–1.3; log10 bias / r without → with loads):
  Hudson −0.58 / 0.31 → −0.56 / 0.43; Delaware −0.51 / 0.43 → −0.47 / 0.62. The pattern improves, the model is still 3–4× too low.
  Salinity in the plumes is identical in both runs (same physics).
- **Ammonium goes negative in both runs (not caused by the loads).** On day 56, 21.9 % of all wet cells have NH₄ < 0
  (50 000 cells < −0.1, 14 000 < −0.5, minimum −2.44 mmol/m³), and 4 915 of 11 951 surface columns are below −0.1 (day 20:
  3.3 % of cells, 1 464 surface columns). P, Z, DOM, the particles and O₂ stay non-negative (bounded WENO); NO₃ dips to −0.1.
  LOBSTER's NH₄ uptake and nitrification act as sources and sinks, which the advection bounds cannot constrain, and
  OceanBioME's `scale_negatives` (conserving the nitrogen and carbon of each cell) is off by default and was not used.
  The chlorophyll and nutrient results of all prototype runs are affected until this is fixed.
- **The wide NO₃ bound brought the bottom-cell overshoot back**: maximum NO₃ 38 on day 56 against 23.9 with the bound at 30.
  The river values (NO₃ up to 103) need a bound above 30 in the mouth cells only.
- Still open: the Connecticut River is not routed; Hudson chemistry is a proxy.

## `bgc_feb19_sc`: RC4USCoast rivers + scaled negative tracers, 59 days from 2019-02-01 (finished 2026-10-23; 2026-10-23)

Run with BGC_BOUNDS="NO₃:0:150,O₂:50:500,DIC:0:2300", Copernicus data, `MAB_BGC_SCALE_NEGATIVES=true`; 4 triton24 nodes, ~33 h (460–550 s
per 6 model hours). Analyses: `bgc/validate_bgc.jl`, `bgc/compare_scaling.jl` (RUN_A=rc4, RUN_B=sc), `bgc/chl_comparison.jl`
(figures and animation of three runs; figures and the mp4 in `~/Data/mab_analysis/bgc/`, not committed).
- **Negative ammonium is gone**: 0.0 % of the wet cells and surface columns on every day (the unscaled run goes from 0.2 % of
  the wet cells on day 3 to 21.9 % on day 56, and from 0.9 % to 41 % of the surface columns); minimum NH₄ ≈ 0.
  Nitrogen total differs from the unscaled run by +4e-7 (day 3) growing to +2.6e-4 (day 56).
- The bottom-cell overshoot is still there (maximum NO₃ 38.0 → 40.1, minimum O₂ 88.6 → 86.5 on day 56).
- **Chlorophyll against the satellite** (matched cell-days, log10 bias / rms / r; model → satellite mean in mg/m³), no loads /
  rivers / rivers + scaling: shelf −0.25 / 0.42 / 0.20, −0.24 / 0.41 / 0.25, −0.31 / 0.45 / 0.30 (model 1.0, 1.0, 0.86 vs 2.39);
  slope +0.14 / 0.29 / 0.30, +0.14 / 0.29 / 0.30, +0.07 / 0.23 / 0.21 (1.0, 1.0, 0.80 vs 0.70); deep +0.32 / 0.40 / 0.48 (twice),
  +0.31 / 0.39 / 0.45 (0.57, 0.57, 0.53 vs 0.30). Region means fall by 15 % (shelf), 17 % (slope), 5 % (deep) with
  scaling over the whole period (the restart test had said 28 %, 30 %, 8 %: part of it was the one-time removal of the accumulated negatives).
  The offshore excess and the weak shelf-to-ocean contrast remain (inherited from the PISCES phytoplankton carbon).
- **Pioneer sites** (CNSM / ISSM, unscaled → scaled): chlorophyll mean 1.12 → 0.91 / 0.95 → 0.81 (fluorometer 0.96 / 0.62), rms 0.58 → 0.30
  / 0.41 → 0.27, r −0.60 → −0.23 / −0.17 → −0.15; the late bloom spike (2.3 mg/m³ at CNSM on days 37–42) is gone but the model is
  too flat (std model/obs 0.29 / 0.35). Nitrate at 7 m: 1.36 → 1.96 / 0.54 → 0.73 (SUNA 6.5 / 1.5); the unscaled run is driven to zero after
  day 40, the scaled one keeps 1.5–4.7 at CNSM. O₂ unchanged (model 96–98 % saturated, sensors 125–131 %).
- **Are the river nitrate inputs too low? Not the main problem (budget, 2026-10-23).** Over the 59 days of `bgc_feb19_sc` the rivers delivered
  6.5e8 mol NO₃ (mean 126 mol/s; plus 4.6e7 mol NH₄ and 1.2e8 mol DON), which spread over the shelf (< 200 m: 9.4e12 m³; < 100 m: 6.8e12 m³) is
  at most 0.07-0.10 mmol/m³ of nitrate, against a shortfall of 0.8 (ISSM) to 4.6 (CNSM) mmol/m³ at the Pioneer sites; closing it would need
  rivers 10-60× larger. The Pioneer sites are upstream of the rivers anyway, and the shortfall is there from the start (day 3: CNSM 1.8 vs 6.4
  SUNA), so it is inherited from the initial and boundary nitrate. Missing river sources would add little: the Connecticut (unrouted) ~17 mol/s,
  the Chesapeake (Susquehanna 1682 m³/s × 101 mmol/m³ ≈ 170 mol/s, mouth just outside the box) ≈ 8.7e8 mol over the run, +0.09 mmol/m³.
  Locally the loads may still matter (plume chlorophyll 3-4× below the satellite) — untested.
- **Where the shelf nitrate shortfall comes from: the Copernicus input (`bgc/nitrate_sources.jl`, 2026-10-23).** Nitrate at 7 m at the Pioneer sites
  (mmol/m³), SUNA / raw Copernicus (nearest 1/4° column, 108 m deep) / the model's pipeline (inpainted, interpolated) / WOA 2023 monthly (× 1.025):
  CNSM: Jan 25-31 6.30 / Jan 2.71 / 2.81 / 6.46; Feb 5.53 / 1.09 / 1.09 / 4.64; Mar 7.48 / 0.90 / 1.18 / 3.17; Apr – / 0.06 / 0.10 / 2.74.
  ISSM: Jan 1.61 / 2.72 / 2.00 / 6.96; Feb 1.51 / 1.03 / 0.92 / 4.91; Mar 1.76 / 1.04 / 0.72 / 3.38; Apr – / 0.05 / 0.06 / 2.82.
  Model initial value (Feb 1, 7 m): CNSM 1.87, ISSM 1.41. Mid-Feb profile at 7, 25, 50, 100 m: CNSM Copernicus 1.1, 1.2, 1.7, 2.5 vs WOA 4.6, 4.7, 5.3, 8.2;
  ISSM Copernicus 1.0, 1.2, 1.6, 2.8 vs WOA 4.9, 4.9, 5.6, 8.6.
  - The pipeline (inpainting, interpolation) does not change the data (1.09 → 1.09); the shortfall at CNSM is in the Copernicus hindcast itself.
  - The hindcast has almost no cross-shelf nitrate gradient (CNSM 1.09, ISSM 1.03 in February) where the sensors show 5.5 against 1.5; WOA matches CNSM
    in Jan-Feb (6.5, 4.6) but is far too high at the inshore ISSM (4.9-7.0 vs 1.5): at 1° it mixes in slope water. Neither product resolves the shelf-break nitrate gradient.
  - The hindcast's nitrate is drawn to ~0 by April (2.7 → 1.1 → 0.9 → 0.06) while SUNA stays at 7.5 (CNSM) in March; at 100 m it has 2.5-2.8 where WOA has 8-8.6.
  - So the model's nitrate follows its Copernicus initial and boundary values; the river loads and the scaling cannot repair this. Candidates: nitrate
    from WOA below the surface layer, a blend of the two, or nitrate from a nitrate-density relation fitted to WOA/SUNA.
- **Nitrate from density does not work (`bgc/nitrate_density_fit.jl`, 2026-10-23).** Nitrate as a non-decreasing function of σ = ρ(Θ, Sᴬ, 0) − 1000 from WOA 2023
  monthly medians (σ bins of 0.1; box -76..-64°E 34..42°N, shallower than 800 m): in-sample r² 0.81-0.86 (rms 2.3-3.3 mmol/m³), but that skill is in the deep water
  (σ > 27: 4-12 mmol/m³); at σ = 26 the table gives 1.0-2.8 (Feb 1.0). Test at 7 m, 2019-02-01 → 03-31, against SUNA (CNSM / ISSM):
  mean SUNA 6.59 / 1.53; from GLORYS T,S 1.57 / 1.56 (bias −5.0 / +0.04, r −0.19 / +0.45); from the optode's own T,S 1.85 / 1.73 (bias −4.7 / +0.2, r +0.59 / +0.46).
  CNSM has 5.5-7.5 mmol/m³ at σ = 26.07, where WOA's box-wide relation gives 1-2.3 (Gulf Stream and Sargasso water of that density has no nitrate). Restricting the fit to
  the shelf columns (WOA data ending shallower than 500 m, 145 samples) gives an almost flat 4-5 in February: CNSM bias −2.5 (r −0.8), ISSM +2.4 (r −0.1/−0.3).
  Nitrate is not a function of density here: at CNSM it is highest when the water is colder and fresher (daily r with T −0.74, with S −0.57; mean NO₃ 6.6 ± 1.2, T 6.9, S 33.3), at ISSM
  the sign flips (r with S +0.43; NO₃ 1.7 ± 0.7). GLORYS σ at the sites also differs from the mooring's (CNSM 25.90 ± 0.39 vs 26.07 ± 0.09, r between the series 0.32). Not
  implemented in the model; the Copernicus nitrate stays. Remaining options: nitrate from WOA (right at CNSM, 4.6-6.5; too high at ISSM, 4.9-7.0), a shelf-water mixing
  line (nitrate vs T and S for the cold fresh shelf water), or nudging the shelf nitrate to the SUNA seasonal means.
- **WOA nitrate in the model (`bgc_feb19_woa`, finished 2026-10-24; as `bgc_feb19_sc` but `MAB_BGC_NO3_SOURCE=woa`; figures in `~/Data/mab_analysis/bgc/validation_woa` and `chl_comparison_woa_vs_copernicus_nitrate`).**
  Copernicus nitrate → WOA nitrate, 59 days: shelf chlorophyll 0.86 → 1.05 mg/m³ (satellite 2.39), log10 bias −0.31 → −0.22, rms 0.45 → 0.40, r 0.30 → 0.29; slope 0.80 → 0.83 and deep 0.53 → 0.52 unchanged;
  whole domain r 0.63 → 0.72, rms 0.40 → 0.37. CNSM nitrate at 7 m 1.96 → 4.59 (SUNA 6.52): bias −4.56 → −1.92, rms 4.76 → 2.53, but r −0.04 → −0.49; CNSM chlorophyll 0.91 → 1.00 (fluorometer 0.96).
  ISSM nitrate 0.73 → 5.44 (SUNA 1.54): bias −0.81 → +3.90, rms 1.19 → 4.01; ISSM chlorophyll 0.81 → 1.01 (0.62). Maximum NO₃ 36.2 on the last day (overshoot unchanged).
  So WOA fixes CNSM and breaks ISSM, as the data comparison predicted. Raising the CNSM nitrate by 2.6 mmol/m³ raised its chlorophyll by 10 % and the shelf mean by 22 %: the shelf chlorophyll deficit
  (satellite 2.3× the model) is not mainly a nitrate deficit in Feb-Mar (light and mixing are the untested candidates; near the coast the satellite is also high).
  Time structure (figure `pioneer_runs.png`): at CNSM the WOA-nitrate run follows SUNA well from day 3 to day ~28 (model 5.5 → 4.9, SUNA 6.4 → 4.5 → 5), then SUNA rises to 7.5-8 (days 28-40) while the model stays near 5 and falls to 2.4-4 after day 45: the late
  nitrate-rich intrusion is missed. At ISSM the model decreases slowly from 6.1 to 4.5 while SUNA oscillates between 0.3 and 2.6 (peaks near days 11, 21, 40, 52): the cross-shelf gradient (SUNA CNSM 5.5-7.5 vs ISSM 1.5 only 35 km apart) is in neither data product and is not created by the model in 8 weeks.
  Off the shelf the WOA and Copernicus nitrate runs have identical chlorophyll (slope, deep): the offshore excess (model 0.5 vs satellite 0.3 mg/m³) is not controlled by nitrate; on the shelf the WOA run is higher by a nearly constant 0.1-0.15 mg/m³ from day 3 on.
  CNSM chlorophyll: neither run reproduces the fluorometer's rise to 1.5 mg/m³ on days 10-25 or its fall to 0.65 afterwards (the WOA run rises to 1.2-1.3 on days 37-49: r −0.42); ISSM chlorophyll stays near 1.0 where the fluorometer falls from 0.8 to 0.3.
  Next diagnostic: the limitation terms of LOBSTER's growth (light, NO₃, NH₄) and the PAR and mixed-layer depth at CNSM, to see what limits the shelf phytoplankton in Feb-Mar.
- Still open: shelf nitrogen supply (nitrate 3× too low at CNSM), Hudson chemistry is a proxy, the Connecticut River is not routed,
  the offshore chlorophyll excess, the bottom-cell overshoot, the first-month oxygen spin-up.
- The forks' `everything` moved to Oceananigans 0.114 on 2026-10-21; `mpi05_bgc` still has the earlier manifest (these runs used it).

- **Scaling negative tracers (2026-10-21; `MAB_BGC_SCALE_NEGATIVES`, `bgc_feb19_rc4sn` = restart of `bgc_feb19_rc4` from its day-35
  checkpoint to day 50 with ScaleNegativeTracers, invalid cells set to zero; `bgc/compare_scaling.jl`, 9 shared frames, days 39–47):**
  - NH₄ < 0 in 0.0 % of the wet cells and 0 % of the surface columns (unscaled: 6.9 % rising to 15.4 % of the wet cells and
    20 % to 34 % of the surface columns; minimum NH₄ −1.9 → −2.3); NO₃ has no negatives either (minimum 5e-4 on day 48).
  - Nitrogen is conserved: the total of NO₃+NH₄+P+Z+DOM+sPOM+bPOM differs from the unscaled run by +6e-6 to +2e-5 (relative).
  - The bottom-cell overshoot is untouched: maximum NO₃ 42.2 and minimum O₂ 78.6 on the last shared day, in both runs.
  - Surface chlorophyll is lower with scaling (shelf 1.07 → 0.77, slope 1.26 → 0.88, deep 0.63 → 0.58 mg/m³): the negative NH₄
    had left more nitrogen in the positive pools (the total, negatives included, is conserved). Against the satellite (log10
    bias / r, unscaled → scaled): shelf −0.20 / +0.03 → −0.35 / −0.11, slope +0.24 / +0.35 → +0.11 / +0.22, deep +0.30 / +0.45 → +0.28 / +0.39.
    The gap is already there at the first shared frame (the restart removed the accumulated negatives at once), so the
    size of the steady effect needs a scaled run from day 0.
  - Cost: about 2× per step (460–550 s per 6 model hours against 207–215 s), so a 59-day run takes about 30 h instead of 15.
  - The earlier chlorophyll conclusions (shelf 1.3–2.4× too low) were flattered by the negative NH₄: with scaling the shelf is
    up to 3× too low, which points again at the shelf nitrogen supply.

## BGC environment update (2026-10-25)

`/t0/workdir/enrique/mpi05_bgc` (triton) updated from Oceananigans 0.113.5 / NumericalEarth 0.8.1 / OceanBioME 0.19.4 / ClimaSeaIce 0.5.10 (fork branch `substep-halo-exchange`) to
**Oceananigans 0.114.0** (fork `everything`, tree 0d6cafe0, includes the always-fma change), **NumericalEarth 0.8.2** (fork `everything` 5cee74ae, tree cff71a8e: contains the sea-ice fix d9dfcbf6 and the Copernicus reader
fe421393, upstream `main` merged), **OceanBioME 0.20.2** (0.19.4 only allows Oceananigans 0.113.3), **ClimaSeaIce 0.5.12** (fork branch `substep-halo-exchange-0.5.12`, as in `mpi05_bgc_114`; NumericalEarth now needs ≥ 0.5.12).
WorldOceanAtlasTools 0.6.5 and CopernicusMarine 0.2.0 unchanged. Old files: `Project.toml.bak-2026-10-25`, `Manifest.toml.bak-2026-10-25` (and the 10-15, 10-20 ones): every `bgc_feb19_*` run used the old stack.
- Only for triton24/triton48 nodes: triton16 has no FMA hardware (the year-long run keeps its frozen `env_old` for that reason).
- Checks on the new stack: precompile and load of all packages, the LOBSTER + scale-negatives column test (identical numbers), the Copernicus/WOA readers (`bgc/nitrate_sources.jl`, identical numbers); the full model script: see the test below.
- `mpi05_bgc_114` (made by the main session on 2026-10-21) and `mpi05_bgc_cmems` (the reader branch before it was on the fork) are untouched; `mpi05_bgc_cmems` is obsolete.

**Full-model test on the new stack (job 12210, `bgc_env114_test2`, 4 triton24 nodes, 1 day from 2019-02-01, same flags as `bgc_feb19_woa`; 2026-10-25).** The first attempt (12208) failed in `dic_from_pH` because OceanBioME 0.20 renamed the `CarbonChemistry` keyword `P` to `water_pressure`; commit 2198a19 makes the keyword version-dependent. The second runs the whole day. At 2019-02-01 00:00 the ranges equal the old stack's (NO₃ 0.0877-25.13, P 0.00185-2.744, O₂ 131-366.4); DIC and Alk differ slightly (DIC 233-2167 against 231-2169, Alk 246-2424 against 244-2419) because the new pressure keyword and package versions move the carbonate solution at the sub-percent level; the 6-hourly ranges of NO₃, P, Z and O₂ afterwards agree to within about 1 %, and the wall time per 6 h of model time is similar (about 410-440 s on both). The new stack is fine for new runs; runs to be compared with the `bgc_feb19_*` set should keep using `mpi05_bgc_old` (= the `.bak-2026-10-25` Project and Manifest, made 2026-10-25).

## Growth limitation diagnostic (2026-10-25)

`MAB_BGC_LIMITATION=true` (`LIMITATION=true` in `bgc_run.sbatch`; `bgc/bgc_limitation.jl`) writes, as daily means of hourly samples in the top `MAB_BGC_LIMITATION_LEVELS` (12, down to 49 m) levels, PAR, the light term PAR/(k_PAR+PAR), the nitrate term NO₃ exp(-ψNH₄)/(NO₃+k_NO₃), the ammonium term, the temperature factor (1 in LOBSTER) and the growth rate μ₀ L_light L_N (1/day) of the phytoplankton, to `<TAG>_bgc_limitation_rank<r>.jld2`. Restart `bgc_feb19_woa_lim` (job 12211): from the day-35 checkpoint of `bgc_feb19_woa` on the old stack (`mpi05_bgc_old`), to day 50. `bgc/limitation_analysis.jl` reads it with the daily T, S (mixed-layer depth from Δσ = 0.03 kg/m³ against 10 m, Roquet simplified EOS) and prints means at CNSM, ISSM and the shelf/slope/deep regions, the share of light- against nutrient-limited cell-days, and maps/time series/profiles.

## Practical notes

- Sandbox: `~/Dropbox/Models/NumericalEarthSandbox` (Mac) and `/t0/workdir/enrique/Models/NumericalEarthSandbox`
  (triton). Run on triton via antares (`ssh antares`, then `ssh triton`), Slurm; exclude triton24 node[38,41]
  (`sbatch --exclude=node[38,41]`). Environments: `/t0/workdir/enrique/mpi05_ib` (MPI over InfiniBand) and
  `mpi05_rivers`; `JULIA_PKG_GC_AUTO=false`. OceanBioME must be added to the environment used.
- NumericalEarth and Oceananigans come from the fork `everything` branches (enri66), which carry the open-boundary and
  distributed fixes.
- Upstream code (Oceananigans, NumericalEarth, OceanBioME) must follow NumericalEarth's restraint rules; PRs need
  Enri's OK before opening.
