# NumericalEarthSandbox

A single environment that combines five in-progress Oceananigans.jl features —
open-boundary radiation, tracer reservoirs, tides, low-pass-filtered output,
and a split-explicit substep-clock fix — with
[NumericalEarth.jl](https://github.com/NumericalEarth/NumericalEarth.jl)'s
realistic regional-ocean machinery (ETOPO bathymetry, GLORYS, ERA5, sea ice),
so they can be tried together before each lands upstream on its own.

None of this is meant for upstream merge. It exists to let a colleague clone
one thing and get the whole combined feature set, while the individual
pieces work their way through review:

| feature | Oceananigans PR |
|---|---|
| `ObliqueRadiation` | [CliMA/Oceananigans.jl#5962](https://github.com/CliMA/Oceananigans.jl/pull/5962) |
| `TracerReservoir` | [CliMA/Oceananigans.jl#5964](https://github.com/CliMA/Oceananigans.jl/pull/5964) |
| `TidalHarmonics` (planet-agnostic), `tidal_forcing`, `tidal_boundary_conditions` | [CliMA/Oceananigans.jl#5970](https://github.com/CliMA/Oceananigans.jl/pull/5970) |
| `LowPassFilter` | [CliMA/Oceananigans.jl#5971](https://github.com/CliMA/Oceananigans.jl/pull/5971) |
| Split-explicit barotropic substep clock fix | [CliMA/Oceananigans.jl#5982](https://github.com/CliMA/Oceananigans.jl/pull/5982) |
| `TPXO10Atlas`, `earth_tidal_harmonics` (Earth's tidal astronomy) | [NumericalEarth/NumericalEarth.jl#681](https://github.com/NumericalEarth/NumericalEarth.jl/pull/681) |

The five Oceananigans branches are merged together, conflict-free, at
[`enri66/Oceananigans.jl#everything`](https://github.com/enri66/Oceananigans.jl/tree/everything).
[`enri66/NumericalEarth.jl#everything`](https://github.com/enri66/NumericalEarth.jl/tree/everything)
adds `TPXO10Atlas` and points its own `Oceananigans` dependency at that
branch (and its `ClimaSeaIce` dependency at
[`enri66/ClimaSeaIce.jl#widen-oceananigans-compat`](https://github.com/enri66/ClimaSeaIce.jl/tree/widen-oceananigans-compat),
since the registered ClimaSeaIce caps its Oceananigans compat below what
the combined branch needs).

## Running

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

This environment's own `[sources]` already point at the right branches, so
no manual `Pkg.develop` is needed.

## Scripts

| script | what it does |
|---|---|
| `03_mab_m2_tide_smoke.jl` | Mid-Atlantic Bight, M2 tide only: astronomical body force + TPXO10 boundary forcing from one `TidalHarmonics`, so they can't drift out of phase. Harmonically analyses the model's free surface and reports skill against TPXO and against NOAA tide gauges. `MAB_TANGENTIAL=oblique` swaps `ObliqueRadiation` in on the tangential velocity component. |
| `02_mab_glorys_obc.jl` | Same MAB box, GLORYS-driven open boundaries (Flather/Chapman fed real exterior data instead of zero). Copied from `NumericalEarth/scripts/02_mab_glorys_obc.jl` — keep in sync there, see that repo's `CLAUDE.md` for the fix history. Unlike `03`, this one doesn't need any of the five bundled branches (it only uses mainline NumericalEarth); it's here because this environment is the easiest way to run it without hand-assembling the package combo, not because it's testing an in-progress feature. `MAB_UEXT=native\|masked\|legacy` selects how the Flather barotropic exterior is computed (`native`, the default, is the most precise). |
| `04_mab_glorys_tides_reservoirs.jl` | `02`'s GLORYS boundaries, but now WITH the M2..Mm tide added on top (additively — see the script's header for why `tidal_boundary_conditions` isn't used directly), `TracerReservoir` on T/S instead of `NormalRadiation`, and `LowPassFilter` de-tided daily/pentad output alongside raw hourly. This one DOES need three of the five bundled branches (`TracerReservoir`, tides, `LowPassFilter`) plus `TPXO10Atlas`, so unlike `02` it can only run here. First validated 2026-09-20 on a 14-day run — see `NumericalEarth/CLAUDE.md`'s session log for the results and the known east-edge spike still worth another look. Also supports checkpoint/restart (`MAB_CHECKPOINT_EVERY`, `MAB_PICKUP`) — see the script's own comment above the `Checkpointer` for the offset needed to keep the de-tided output continuous across a restart (validated on a 30-day run, `mab_1month_continuous`, 2026-09-21). |

Its includes (`tidal_harmonics.jl`, `tpxo.jl`, `tpxo_boundaries.jl`,
`noaa_harcon_mab.jl`, `variable_bottom_drag.jl`,
`kepsilon_tuple_closure_patch.jl`, `glorys_profiles.jl`, `glorys_bathymetry.jl`)
are local helpers, not part of either package — the skill-analysis machinery
in particular (harmonic analysis, atlas readers) is deliberately not
something the package provides.

Four animation scripts, all reading a run's saved `.jld2` output via
`MAB_TAG`: `animate_ssh.jl` (raw hourly SSH — will visibly carry the tide if
`04` produced it), `animate_ssh_detided.jl`, `animate_sst_detided.jl`, and
`animate_velocity_detided.jl` (filled contours + arrows2d for surface speed
and direction) — the latter three are filled-contour/vector plots from
`04`'s `LowPassFilter` daily output (`{TAG}_eta_daily.jld2`/
`{TAG}_surface_daily.jld2`), so they only cover the days the daily filter
actually produced valid output for. That span depends on whether/how the run
was restarted: a single unbroken run only trims ~3 days off each end (a
14-day run gets days 3–11); a checkpoint/restart run gets a real gap unless
the restart checkpoint was taken far enough before the split (see `04`'s
comment on this — a 30-day run restarted at day 9 for a day-15 split covers
the full days 3–27 with no gap, `mab_1month_continuous`).

Plus `snapshot_fields.jl`, a static SST/SSH/surface-speed three-panel PNG at a
few chosen days (`MAB_SNAPSHOT_DAYS=34,60,87`, default first/middle/last of
the de-tided record) — faster than scrubbing an animation when spot-checking
a run. Extended `mab_1month_continuous` to 90 days as `mab_3months`
(resuming from its day-30 checkpoint, since that's the one with the
validated-continuous offset — `mab_1month_restart`'s day-12.5 checkpoint
predates that fix and has a real gap); snapshots at days 34/60/87 show a
coherent Gulf Stream meandering and shedding rings the whole way, and the
restart seam checks out clean against the raw hourly record.

Run with `julia -t 8 --project=. scripts/03_mab_m2_tide_smoke.jl`,
`scripts/02_mab_glorys_obc.jl`, or `scripts/04_mab_glorys_tides_reservoirs.jl`;
see each script's own comment block for its environment-variable knobs
(`03`: `MAB_DAYS`, `MAB_NZ`, `MAB_CLOSURE`, `MAB_TANGENTIAL`, ...; `02`/`04`:
`MAB_DAYS`, `MAB_UEXT`, `MAB_TANGENTIAL`, `MAB_TAU_IN`, `MAB_MATCH_BATHY`, ...;
`04` additionally: `MAB_TIDE_CONSTITUENTS`, `MAB_TIDE_RAMP`,
`MAB_RESERVOIR_L_IN`, `MAB_RESERVOIR_L_OUT`).

### Data

TPXO10-atlas-v2 is licensed and cannot be redistributed here: request it at
[tpxo.net](https://www.tpxo.net), then point `ENV["TPXO_DIR"]` at the
unpacked files (default `~/Data/TPXO10_atlas_v2_nc`). Nothing else script
`03`'s default smoke test needs requires a download — the NOAA gauge
constants are inlined, and the optional GLORYS stratification profiles
(`MAB_STRAT=deep`/`shelf`) are off by default. `04` needs TPXO too, for its
tidal boundary constants.

`02_mab_glorys_obc.jl` and `04_mab_glorys_tides_reservoirs.jl` default
`DATA_DIR` to `~/Data/NumericalEarth` (override with `MAB_DATA_DIR`) — the
same out-of-Dropbox cache the `NumericalEarth` repo's own scripts use (moved
there 2026-09-21; it used to live under that repo's `data/`, synced through
Dropbox), so this doesn't re-download the several-GB GLORYS/ETOPO/ERA5 data
that's already sitting there.
