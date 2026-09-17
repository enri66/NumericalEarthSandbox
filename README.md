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

Its includes (`tidal_harmonics.jl`, `tpxo.jl`, `tpxo_boundaries.jl`,
`noaa_harcon_mab.jl`, `variable_bottom_drag.jl`,
`kepsilon_tuple_closure_patch.jl`, `glorys_profiles.jl`) are local helpers,
not part of either package — the skill-analysis machinery in particular
(harmonic analysis, atlas readers) is deliberately not something the
package provides.

Run with `julia -t 8 --project=. scripts/03_mab_m2_tide_smoke.jl`; see the
comment block at the top of the script for its environment-variable knobs
(`MAB_DAYS`, `MAB_NZ`, `MAB_CLOSURE`, `MAB_TANGENTIAL`, ...).

### Data

TPXO10-atlas-v2 is licensed and cannot be redistributed here: request it at
[tpxo.net](https://www.tpxo.net), then point `ENV["TPXO_DIR"]` at the
unpacked files (default `~/Data/TPXO10_atlas_v2_nc`). Nothing else the
default smoke test needs requires a download — the NOAA gauge constants
are inlined, and the optional GLORYS stratification profiles
(`MAB_STRAT=deep`/`shelf`) are off by default.
