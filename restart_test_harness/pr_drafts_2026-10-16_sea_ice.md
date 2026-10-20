
## Opened 2026-10-20 (from the main MAB session, Enri's OK)

Each branch cherry-picked onto current upstream main; for each, the new test fails on main without the fix and passes
with it (bounds checking on); quality-assurance and the touched test files pass.
- CliMA/Oceananigans.jl#6291 zstar-pressure-halos-pr (was `zstar-pressure-vertical-halos`): z* pressure vertical halos
- CliMA/Oceananigans.jl#6292 bounds-limiter-halos-pr (was `bounds-limiter-tracer-halos`): tracer halos before bounded WENO limiters
- CliMA/Oceananigans.jl#6293 clocked-fts-window-pr (was `clocked-halo-fill-fts-window`): clocked halo fills update BC series
- NumericalEarth/NumericalEarth.jl#790 glorys-sea-ice-zero-pr (was `glorys-sea-ice-open-water`): missing GLORYS ice = 0
- NumericalEarth/NumericalEarth.jl#791 interpolated-no-clock-pr (was `cyclic-fts-boundary-window`): clockless Interpolated read
Still held: N1/N2 (after NumericalEarth #746), N3 frazil (Enri to decide). ClimaSeaIce#188 open.
