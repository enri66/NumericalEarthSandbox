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
