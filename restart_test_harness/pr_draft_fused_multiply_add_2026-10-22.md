# PR draft (2026-10-22). SUBMITTED 2026-10-22 as CliMA/Oceananigans.jl#6313 (Enri's OK), before the maintainers had replied on #6296.

Repo / base: CliMA/Oceananigans.jl `main` (77f96ece9 when prepared)
Head: enri66/Oceananigans.jl `fused-multiply-add-option` (ad3f135b2, rebased on main 77f96ece9; local worktree ~/dev/Oceananigans-pr-fma)
Diff: Project.toml +2, src/Advection/Advection.jl +30/-1, test/unit/fused_multiply_add.jl +34 (new), docs/src/simulation_tips.md +16

**Title:** `Add a fused_multiply_add preference: muladd or fma in the advection kernels`

**Body:**

## What this adds

A `fused_multiply_add` preference for Oceananigans with the values `"muladd"` (the default, so nothing changes unless it
is set) and `"fma"`. The Advection module's `@fused_muladd` expands as `MuladdMacro.@muladd` does; with `"fma"` every `muladd`
it generates is replaced by `fma`, which always fuses. `BFloat16` takes its `fma` from BFloat16s.jl
(JuliaMath/BFloat16s.jl#106, merged and not yet released), so with `"fma"` it needs a BFloat16s newer than 0.6.2.
An unknown value is an error. It is set with
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

`test/unit/fused_multiply_add.jl` checks that the `muladd` calls of the macro's expansion are replaced by `fma` or kept,
and that an unknown preference is an error. The weno_smoothness,
bounds_preserving_advection, materialize_advection and immersed_advection tests pass with each value of the preference.

🤖 Generated with [Claude Code](https://claude.com/claude-code)


