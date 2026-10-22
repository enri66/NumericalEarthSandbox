# PR draft (2026-10-22), NOT SUBMITTED. Submit only if the Oceananigans developers agree on CliMA/Oceananigans.jl#6296.

Repo / base: CliMA/Oceananigans.jl `main` (77f96ece9 when prepared)
Head: enri66/Oceananigans.jl `fused-multiply-add-option` (ad3f135b2, rebased on main 77f96ece9; local worktree ~/dev/Oceananigans-pr-fma)
Diff: Project.toml +2, src/Advection/Advection.jl +30/-1, test/unit/fused_multiply_add.jl +34 (new), docs/src/simulation_tips.md +16

**Title:** `Add a fused_multiply_add preference: muladd or fma in the advection kernels`

**Body:**

## What this adds

A `fused_multiply_add` preference for Oceananigans with the values `"muladd"` (the default, so nothing changes unless it
is set) and `"fma"`. The Advection module's `@muladd` expands as `MuladdMacro.@muladd` does; with `"fma"` every `muladd`
it generates is replaced by a fused multiply-add that always fuses. `BFloat16`, which has no `fma`, fuses in `Float32`,
and mixed argument types are promoted first. An unknown value is an error. It is set with
`Preferences.set_preferences!(Oceananigans, "fused_multiply_add" => "fma")` and read when Oceananigans is compiled; the
choice is documented in the simulation tips.

Adds `Preferences` as a dependency.

## Why

Refs #6296. `muladd` leaves it to the compiler whether to fuse. On the CPU, the kernels compiled for a single rank and
for a rank of a distributed run can make different choices, and a few cells differ in the last bit; the bounds-preserving
limiter and the following steps turn that into differences of up to 1e-13. With `"fma"` the layouts agree exactly.

Regional hydrostatic case (z*, partial cells, `BoundsPreservingWENO` tracers, split-explicit free surface, open boundaries,
tides, sea ice), 3 steps, Float64, Xeon Gold 6248R, comparing every prognostic field with the single-rank run:

| `fused_multiply_add` | 1×2 vs 1 rank | 2×1 vs 1 rank |
|---|---|---|
| `"muladd"` | 26 of 52 field-step comparisons differ (up to 7e-14) | not run |
| `"fma"` | identical (52 of 52) | identical (52 of 52) |

Timing, 30 steps of the ocean-only case on one rank: 3.17 and 3.22 s per step with `muladd`, 3.23 and 3.44 s with `fma`.

## Tests

`test/unit/fused_multiply_add.jl` checks that the fused operation rounds once in Float32 and Float64, accepts mixed and
BFloat16 arguments, and that the `muladd` calls of the macro's expansion are replaced. The weno_smoothness,
bounds_preserving_advection, materialize_advection and immersed_advection tests pass with each value of the preference.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
