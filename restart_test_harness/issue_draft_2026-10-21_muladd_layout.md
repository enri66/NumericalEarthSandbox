# Issue draft (2026-10-21), for Enri to read before anything is posted

Repo: CliMA/Oceananigans.jl. Type: issue, not a PR. Nothing here has been posted.

## Title

`@muladd` in the Advection kernels gives last-bit differences between rank layouts on the CPU; explicit `fma` removes them

## Body

**Summary.** A regional hydrostatic model (z*, partial cells, `BoundsPreservingWENO` tracers, split-explicit free surface,
open boundaries, tides) run on 1 rank and on a 1×2 rank layout differs in a handful of cells by one ulp after the first step
(T: 8.9e-16 at one cell of 1.78M); by the third step 150-ish cells differ at 1e-16 to 1e-14. A 2×1 layout is identical to the
single rank. Making `@muladd` plain arithmetic in the Advection module, or replacing `muladd` by `fma`, makes all layouts
identical (see below). Is bit-identity across rank layouts a goal, and would `fma` be an acceptable change?

**Where it starts.** I merged upstream main (0.114.0) step by step into a fork and compared the 1-rank and 1×2 runs at every
point. They are identical up to and including #6258. The first difference appears with #6209 ("fused multiply-adds in the
WENO reconstruction"), which added `@muladd` to `biased_p`, the global smoothness indicator and `zweno_ratio_term`. Other
`@muladd` sites in `centered_reconstruction.jl`, `upwind_biased_reconstruction.jl` and `weno_interpolants.jl` predate it.

**What I think is happening (not proven).** `muladd` lets the compiler choose between a fused and an unfused operation. The 1×2
run compiles different kernel specializations than the 1-rank run (different grid topology type parameters), and the choice
can differ, so the same inputs round differently in a few cells.

**Evidence.** Same configuration, 1 rank vs 2 ranks (CPU, Float64, Xeon Gold 6248R, 3 steps with sea ice, 140-158
field comparisons per case):

| change | 1×2 vs 1 rank | 2×1 vs 1 rank |
|---|---|---|
| upstream code (`muladd`) | 1 cell at step 1, ~1e-14 after 3 steps | identical |
| bounds limiter switched off | identical | n/a |
| `@muladd` replaced by plain arithmetic | identical | identical |
| `@muladd` replaced by `fma` (macro below) | identical | identical |

The single-cell differences also appear with z* off and with grid-fitted bottom, at other cells, and only with the
bounds-preserving limiter active.

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
