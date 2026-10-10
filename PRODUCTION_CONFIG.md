# MAB production configuration (script 05)

The settings agreed for production runs of `scripts/05_mab_glorys_tides_reservoirs_distributed.jl`, with the
reasons and the follow-ups still open. Agreed 2026-10-09/10.

## Target simulation

The 1/60° run with 100 levels (`hr100_fx3` on triton: 25 triton16 nodes, 5 × 5 ranks, dt 50 s) is the target
simulation. Analysis: `scripts/analysis_package.sh` (surface, MLD/thermocline sections, Argo, Pioneer moorings).

## Settings

| setting | value | why |
|---|---|---|
| vertical coordinate | `MAB_ZSTAR=true` | freshwater adds volume; improved the slope/deep MLD bias in the 1/12° tests |
| bottom | `MAB_BOTTOM=partial` | partial bottom cells (cost about 1.2× per step) |
| open boundaries with z* | fixed boundary treatment (sandbox 73a0081) | GLORYS transport and the consistent normal velocity account for the z* stretching |
| tracer advection | `MAB_TRACER_ADVECTION=weno7` | harmless; it does **not** fix the bottom-cell overshoot (see below) |
| rivers | `MAB_RIVERS=true` (`mpi05_rivers` env) | GloFAS daily discharge, hand-placed Hudson and Delaware mouths |
| surface salinity restoring | `MAB_SSS_PISTON=0.5` (m/day) toward GLORYS | holds the shelf, slope and deep salinity |
| restoring near rivers | `MAB_SSS_RIVER_MASK=30,100` (km), `MAB_SSS_RIVER_QMIN=20` (m³/s) | no restoring within 30 km of mouths with ≥ 20 m³/s, full beyond 100 km: uniform restoring weakened the plumes by 20–30 %, the masked plume equals the unrestored one near the mouths |
| river-mouth mixing | default (0.1 m²/s over the top 10 m) | keeps a plume in one surface cell from driving the salinity to zero |
| triton24 nodes | `sbatch --exclude=node[38-41]` | node38's socket-0 memory runs at ~1/45 of normal bandwidth; 39–41 not yet diagnosed (`runs/NODES_README.txt`) |

## Open items

- **Masked-restoring check:** the corrected no-river control `sr_ctl_msk2` should confirm that the shelf salinity
  away from the mouths does not drift where the restoring is weaker; if it does, reduce the outer radius.
- **Mask radii by season:** 30/100 km suits September's small plumes; for spring freshets (Hudson plume 100+ km
  along the New Jersey coast) widen the mask, e.g. 50/150 km, or scale it with discharge.
- **High-discharge test:** September discharge is low (Hudson ~250 m³/s), so it cannot discriminate between river
  options; repeat the river tests in a spring-freshet period (global GloFAS archive on antares, GLORYS/ERA5 for the
  period to be brought to triton).
- **Rivers with momentum:** with z*, add rivers as a lateral mass and momentum injection at the mouth instead of a
  surface freshwater flux (the z* volume input alone carries no momentum and barely changes the plumes).
- **Land mask near the mouths:** do the production land-mask clean-up first (including the suspected channel at
  Cape Cod, 69.94°W 41.81°N, where the 1/60° run has |v| up to 3.9 m/s and z* runs η down to −2.4 m).
- **Bottom-cell overshoot:** shelf bottom cells heat and salt themselves with WENO tracer advection where vertical
  mixing at the bottom is weak (CATKE; worse with a weak constant diffusivity, absent with first-order upwind
  advection or the Ri-based closure). Not yet seen in the 1/60° run. Parked: implications for the WENO scheme to be
  considered.
