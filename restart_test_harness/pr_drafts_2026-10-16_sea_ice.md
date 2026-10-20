# Sea-ice PR drafts (2026-10-16)

> Superseded for status: these drafts were written on 2026-10-16. What was opened from them, and what was held back or
> raised as an issue instead, is recorded in the dated sections after the horizontal rule below.


Prepared from the enri66 fork fixes. Each PR body ends with the Claude Code footer only.

## C1 — ClimaSeaIce wall halos at rank-edge corners: no PR needed

The fork fix (enri66/ClimaSeaIce.jl `extended-halo-corners`, 81b196b) filled the wall halos across the connected halos
during the split-explicit substeps. Upstream already does this: Oceananigans #6157 (released in 0.114.0) makes lateral
halo fills cover the connected halos (`connected_fill_extent`), and ClimaSeaIce `main` requires Oceananigans 0.114.
A serial-vs-distributed probe with walls meeting the rank edges (topology Periodic × Bounded, 2 × 1) differs by
6.6e-2 m/s in u on ClimaSeaIce 0.5.10 + Oceananigans 0.113.6 and is exact on ClimaSeaIce `main` + Oceananigans
0.114.0. The fork commit also defines a `connected_range` that would now clash with `Oceananigans.Grids.connected_range`.
Optional: a test-only PR adding a Bounded topology to ClimaSeaIce's distributed tests (they are all periodic now).

## C2 — ClimaSeaIce: `SplitExplicitSolver(grid; extend_halos = false)`

- Repo / base: CliMA/ClimaSeaIce.jl `main` (4e7958d)
- Head: enri66/ClimaSeaIce.jl `split-explicit-halo-exchange` (879cd00; local worktree ~/dev/ClimaSeaIce-pr-halo-exchange)
- Diff: +41 / −17 (src +23/−12, tests +18/−5)
- Tests run on the Mac (4 ranks, Oceananigans 0.114.0): TEST_GROUP distributed, timestepping, checkpointing,
  quality_assurance pass. Probe: extend_halos = false exact vs serial for 2×1, 1×2, 2×2, Periodic×Bounded and
  Bounded×Bounded, 10 and 120 substeps. Control: 120 substeps on `main` fails with "halo=243 must be ≤ size=40".

**Title:** `SplitExplicitSolver: extend_halos = false fills the velocity halos at every update`

**Body:**

On a distributed grid the split-explicit solver extends the halos in `Connected` directions to `2substeps + 3`
cells and substeps into them without communication. Those halos have to fit inside each rank, so the substep
counts used in practice (120–240) cannot run on regional partitions: 120 substeps need 243-cell halos.

This adds `extend_halos` to `SplitExplicitSolver`, with the same meaning as in Oceananigans' `SplitExplicitFreeSurface`:

- `extend_halos = true` (default): unchanged.
- `extend_halos = false`: the halos are not extended and the velocity halos are filled, including communication,
  after every velocity update.

The new distributed test runs the existing 4-rank equality test with 120 substeps and `extend_halos = false`
(25 cells per rank). The probe below also checks walls meeting rank edges (Periodic × Bounded and Bounded × Bounded,
2×2) against serial; all exact.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

## N1 — NumericalEarth: fill the sea-ice and temperature halos before the first coupled fluxes

- Repo / base: NumericalEarth/NumericalEarth.jl `main`, but only after #746 merges (see below)
- Head: enri66/NumericalEarth.jl `distributed-sea-ice-startup-halos` (196e4f49; local worktree
  ~/dev/NumericalEarth-pr-startup-halos), currently stacked on #746's head `exchanger-temperature-halo` (6aeae1e5)
- Diff over #746: src +4/−1 (SeaIces.jl, earth_system_model.jl), test +43/−1
- Depends on #746: without it the coupled model's initial fluxes already differ from serial at rank edges with no
  sea ice at all (τx 3e-6 in an ocean-only probe), which hides this fix. Once #746 merges, rebase onto `main`.
- Tests run on the Mac (4 ranks, Oceananigans 0.113.6, ClimaSeaIce 0.5.12): new testset fails on #746 alone
  (initial τx off by 1.4e-4 of 2.9e-4) and passes with N1 (6/6 on every rank).

**Title:** `Fill the sea-ice and ocean temperature halos before the first coupled fluxes`

**Body:**

The coupled model computes its first interface fluxes when it is built. On a distributed grid two fields reach that
computation with stale halos next to the rank edges:

- The ice concentration and thickness: `set!(sea_ice.model, ...)` fills the interior only, and the momentum fluxes
  average the ice concentration to the velocity faces. The sea-ice exchanger now calls `update_state!(sea_ice.model)`
  when it is initialized, which fills them.
- The ocean temperature after `above_freezing_ocean_temperature!` raises supercooled water to the freezing point over
  the interior: its halos are refilled.

The test builds a supercooled ocean under an ice patch that crosses the rank edges and checks that the initial net
ocean stresses equal the serial ones on 4×1, 1×4 and 2×2. Before this change they differ by up to half of their
magnitude at the rank edges.

Depends on #746.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

## N2 — NumericalEarth: refresh the ice-ocean drag's ocean velocities

- Repo / base: NumericalEarth/NumericalEarth.jl, base = N1's branch until N1 merges (one level deep), then `main`
- Head: enri66/NumericalEarth.jl `distributed-sea-ice-ocean-drag` (d5f049f0; local worktree
  ~/dev/NumericalEarth-pr-ocean-drag)
- Diff over N1: src +2/−1 (sea_ice_ocean_fluxes.jl), test +22/−3
- Tests run on the Mac (4 ranks): new testset fails with N1 alone and passes with N2.

**Title:** `Refresh the ice-ocean drag's ocean velocities before computing the interface stress`

**Body:**

With split-explicit sea-ice dynamics on a distributed grid, the ice-ocean stress keeps halo-extended copies of the
ocean surface velocities, which `update_external_stress!` refreshes. ClimaSeaIce calls it inside the ice momentum
step, after the coupled model has computed the ice-ocean drag coefficients, so the first coefficients read unfilled
halos at the rank edges. `compute_sea_ice_ocean_fluxes!` now calls `update_external_stress!` before computing them.

The test checks that the ice-ocean drag coefficients equal the serial ones on 4×1, 1×4 and 2×2, on the velocity faces
off the walls. (On the wall faces at the domain corners the halo-extended copies still differ with Oceananigans 0.113;
the ice velocity there is zero, and Oceananigans 0.114 fills those corners.)

Depends on the previous PR (N1).

🤖 Generated with [Claude Code](https://claude.com/claude-code)

## Not in these PRs

- N3 (frazil: `update_state!(ocean.model)` after the frazil adjustment, fork `sea-ice-halos-frazil`): held back;
  cost and open-boundary side effects unchecked.
- I1, I2 (separate session): fork branches `glorys-sea-ice-open-water` (7870435d), `cyclic-fts-boundary-window`
  (4482e8f5) on NumericalEarth; Oceananigans `clocked-halo-fill-fts-window` (b7788078).
- I5 (GSL ocean 1- vs 2-rank): under investigation in a separate session.

---

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
