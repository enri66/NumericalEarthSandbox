# Open-boundary regression cases (2026-10-20)

Checks of Oceananigans' open-boundary behaviour that our PRs fixed (#6092, #6113, #6185) and of inflow accuracy, written to
run on both the current API (`OPEN` unset: `Bounded` sides) and Simone's refactor stack (#6275 -> #6276 -> #6277, branch
`ss/open-state-steps`; `OPEN=true`: `Open` sides). `a_*` needs `mpiexec -n 4` (`ONLY=22` runs the 2×2 layout alone).

| case | main | #6185 | refactor stack |
|---|---|---|---|
| a distributed oblique vs serial | 1.15e-5 (both layouts) | exact | 1.15e-5 (both layouts) |
| b repeated `update_state!` unchanged | fails (QAB2) | – | passes (QAB2, SRK3, both scheme pairs) |
| c time-varying inflow, WENO-5 `Value` | 0.010 | – | 0.010 |
| d restart exact | QAB2 1.2e-7, SRK3 exact | – | exact |
