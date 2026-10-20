
## Opened 2026-10-20 (from the main MAB session, Enri's OK)

Each branch cherry-picked onto current upstream main; for each, the new test fails on main without the fix and passes
with it (bounds checking on); quality-assurance and the touched test files pass.
- CliMA/Oceananigans.jl#6291 zstar-pressure-halos-pr (was `zstar-pressure-vertical-halos`): z* pressure vertical halos
- CliMA/Oceananigans.jl#6292 bounds-limiter-halos-pr (was `bounds-limiter-tracer-halos`): tracer halos before bounded WENO limiters
- CliMA/Oceananigans.jl#6293 clocked-fts-window-pr (was `clocked-halo-fill-fts-window`): clocked halo fills update BC series
- NumericalEarth/NumericalEarth.jl#790 glorys-sea-ice-zero-pr (was `glorys-sea-ice-open-water`): missing GLORYS ice = 0
- NumericalEarth/NumericalEarth.jl#791 interpolated-no-clock-pr (was `cyclic-fts-boundary-window`): clockless Interpolated read
Still held: N1/N2 (after NumericalEarth #746), N3 frazil (Enri to decide). ClimaSeaIce#188 open.

## N3 (frazil update_state!) checked 2026-10-20: not submitted as drafted
`default_sea_ice()` is `FreezingLimitedOceanTemperature`, so the frazil adjustment runs in every default coupled model, the
MAB production runs included. N3 adds a second `update_state!(ocean.model)` per coupled step: measured 19.6 ms on a 133.6 ms
coupled step (+15 %, CPU, 120 × 80 × 50, CATKE, `frazil_cost.jl` in the session scratchpad). It also runs the boundary-condition
updates twice per step (radiation schemes anchor at stage ≤ 1; the per-iteration anchoring of Oceananigans #6092 is not merged).
The gain is the ~1e-8 serial-vs-distributed difference at rank edges. Kept on the fork (`sea-ice-halos-frazil`) for exactness tests.

## N1, N2 opened and N3 raised as an issue (2026-10-20, Enri's OK)
- NumericalEarth/NumericalEarth.jl#792 (N1, branch `sea-ice-startup-halos-pr`) and #793 (N2, `sea-ice-ocean-drag-pr`): stacked on #746,
  built on current main with #746's commits underneath. N1's test fails without N1's two source changes (4 checks) and
  passes; N2's fails without its change (4 checks) and passes (4 ranks).
- Nudge posted on #746.
- NumericalEarth/NumericalEarth.jl#794: N3 as an issue with a 55-line reproducer (`frazil_repro.jl`, 2 ranks). On the Oceananigans
  fork stack (everything-gsl-exact): u 2.2e-6, v 1.3e-7, T 1.8e-7, S 6.6e-8 without the extra `update_state!`, ~1e-17 with it.
  On plain upstream Oceananigans the NE-Ice harness shows u 6e-6, v 5e-5 at the rank edge with or without N3 (other distributed
  differences), so N3's effect only shows with the other fixes applied.
