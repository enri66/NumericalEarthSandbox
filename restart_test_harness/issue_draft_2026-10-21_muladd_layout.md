# Issue draft (2026-10-21), for Enri to read before anything is posted

Repo: CliMA/Oceananigans.jl. Type: issue, not a PR. Nothing here has been posted.

## Title

`@muladd` in the Advection kernels gives last-bit differences between rank layouts on the CPU; explicit `fma` removes them

## Body

**Summary.** A regional hydrostatic model (z*, partial cells, `BoundsPreservingWENO` tracers, split-explicit free surface,
open boundaries, tides) run on 1 rank and on a 1×2 rank layout differs in a few cells by about one ulp after the first step,
growing to 1e-13 within two steps. A 2×1 layout is identical to the single rank. Making every `@muladd` in the Advection module
plain arithmetic, or replacing `muladd` by `fma`, makes all layouts identical (see below). Is bit-identity across rank layouts
a goal, and would `fma` be an acceptable change?

**Where it starts.** I merged upstream main (0.114.0) step by step into a fork and compared the 1-rank and 1×2 runs at each
point. They are identical up to and including #6258. They differ at the next point I tested, which includes #6253
(regridding) and #6209 ("fused multiply-adds in the WENO reconstruction"); I did not separate those two. #6209 added `@muladd`
to `biased_p`, the global smoothness indicator and `zweno_ratio_term`; other `@muladd` sites in the Advection module predate it.
Making only those three sites plain removes most of the difference but not all of it, which fits a compiler-dependent
effect rather than one bad site.

**What I think is happening (not proven).** `muladd` lets the compiler choose between a fused and an unfused operation. The 1×2
run compiles different kernel specializations than the 1-rank run (different grid topology type parameters), and the choice
can differ, so the same inputs round differently in a few cells.

**Evidence.** Same configuration, 1 rank vs 2 ranks (CPU, Float64, Xeon Gold 6248R, 2-3 steps), differences in the
prognostic fields (η, u, v, T, S, e, and for the 3-step cases sea-ice fields):

| code | 1×2 vs 1 rank | 2×1 vs 1 rank |
|---|---|---|
| upstream `muladd` (state at #6209/#6253) | step 1: 1 η, 69 u, 6 S cells differ (S up to 4e-14); step 2: 59 η, 1224 u (up to 1.4e-13), 171 S cells | identical |
| only the three #6209 sites made plain | step 1: one T cell, 8.9e-16 (1 ulp); by step 3 about 150 u cells, 12 T, 11 S, up to 1.4e-14 | identical |
| same, bounds-preserving limiter switched off | identical (1 step) | not run |
| same, z* off, or grid-fitted bottom | still one T cell at a different place | not run |
| every `@muladd` in Advection plain | identical (3 steps, 140-158 comparisons) | identical |
| every `@muladd` in Advection replaced by `fma` (macro below) | identical (3 steps, 140 comparisons) | identical |

So the remaining single-cell difference needs the limiter, which evaluates WENO reconstructions inside its own kernel.

```julia
fma_expr(x) = x
fma_expr(x::Symbol) = x === :muladd ? :fma : x
fma_expr(x::Function) = x === Base.muladd ? Base.fma : x
fma_expr(x::GlobalRef) = x.name === :muladd ? GlobalRef(Base, :fma) : x
fma_expr(ex::Expr) = Expr(ex.head, map(fma_expr, ex.args)...)
macro muladd(ex)
    expanded = macroexpand(__module__, Expr(:macrocall, GlobalRef(MuladdMacro, Symbol("@muladd")), __source__, ex))
    return esc(fma_expr(expanded))
end
```

**Caveats, please read.**
- The case is a regional production configuration, not a small script. I could not reproduce the difference in a small
  distributed box test (an immersed strip, bounded WENO, 240 ice substeps was identical either way).
- The fork I tested contains halo fixes that are not in upstream yet (CliMA/Oceananigans.jl#6291, #6292, #6293), without
  which the layouts differ by far more than a last bit. So I cannot say that bare upstream main shows exactly this.
- `fma` has no hardware support on older CPUs (Sandy Bridge), where it is a slow software routine; the results are still
  identical across layouts.
- I have not measured the speed of `fma` versus `muladd`.

**Possible relevance to GPUs (untested).** Multi-rank GPU runs use the same domain decomposition. If different kernel
specializations are fused differently there, they could show the same layout dependence. I have no GPU data; the same
1-rank vs 2-rank comparison on a GPU would settle it. `fma` is native on NVIDIA GPUs.

**Not covered.** Global reductions (CG solver, `all_reduce`) change summation order with the layout and cannot be
bit-identical; this issue is about elementwise kernels only.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
