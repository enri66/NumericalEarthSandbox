# Comment on CliMA/Oceananigans.jl#6296 (2026-10-21). POSTED with Enri's OK: https://github.com/CliMA/Oceananigans.jl/issues/6296#issuecomment-6327377416

Two follow-ups after trying `fma` in the fork.

**1. `fma` needs a BFloat16 fallback.** `Base.fma` has no method for `BFloat16` (it errors with "fma not defined for
BFloat16"), and Oceananigans supports that type. With a plain `muladd` -> `fma` replacement, `test/unit/weno_smoothness.jl`
errors for BFloat16 (Float32, Float64 and BigFloat pass). Mixed argument types (for example an `Int` coefficient with
BFloat16 operands) also reach `fma` before promotion. What I use in the fork:

```julia
@inline fused_multiply_add(a, b, c) = fused_multiply_add(promote(a, b, c)...)
@inline fused_multiply_add(a::T, b::T, c::T) where T = fma(a, b, c)
@inline fused_multiply_add(a::BFloat16, b::BFloat16, c::BFloat16) = BFloat16(fma(Float32(a), Float32(b), Float32(c)))
```

and the macro replaces `muladd` in MuladdMacro's expansion by `fused_multiply_add`. With this,
`unit/weno_smoothness`, `unit/bounds_preserving_advection`, `unit/materialize_advection` and `time_stepping/immersed_advection`
pass (CPU, including BFloat16).

**2. Layouts and speed with `fma`.** Same regional configuration as above, with sea ice, 3 steps: the 1-rank, 2×1 and 1×2
runs are identical in all 104 field comparisons. For speed I timed 30 steps after 5 warm-up steps of the ocean-only
configuration on one rank (Xeon Gold 6248R, 24 threads, Float64), two runs each:

| `@muladd` | s per step |
|---|---|
| upstream (`muladd`) | 3.17, 3.22 |
| `fma` | 3.23, 3.44 |

That is about 4% slower on average, but the spread between the two runs of the same code is about as large (0.05-0.2 s),
so I would not read more into it than "no large slowdown here". I did not measure on a GPU.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
