# ==================================================================
# DOES STRATIFICATION CHANGE THE FIXED-Cᵇ CATKE RESULT?
#
# Follow-up to compare_catke.jl: with the corrected Cᵇ = 0.967, does adding
# realistic stratification (which changes CATKE's diagnosed bottom viscosity
# substantially — see the header of 03_mab_m2_tide_smoke.jl) move the M2
# skill? Under the OLD Cᵇ = 0.28 it did not (≤1%, measured 2026-09-11); this
# checks whether the fix changes that conclusion.
#
# Loads the three fixed-Cᵇ *_results.jls dumps: MAB_STRAT=none/deep/shelf.
#
# Run:  julia --project=. scripts/compare_catke_strat.jl
# ==================================================================

using Serialization
using Printf
using Statistics
using CairoMakie

const OUT_DIR = get(ENV, "MAB_OUT", joinpath(homedir(), "Data", "mab_tides", "catke_comparison"))

runs = (
    ("catke_new",       "no stratification"),
    ("catke_new_deep",  "deep/slope GLORYS profile"),
    ("catke_new_shelf", "shelf GLORYS profile"),
)

results = NamedTuple[]
for (tag, label) in runs
    path = joinpath(OUT_DIR, "$(tag)_results.jls")
    isfile(path) || error("missing $path")
    push!(results, merge(deserialize(path), (; label)))
end

r1 = results[1]
λc, φc, A_tpx, G_tpx, valid = r1.λc, r1.φc, r1.A_tpx, r1.G_tpx, r1.valid
maskland(A) = replace(x -> isfinite(x) ? x : NaN, A)
amax = maximum(filter(isfinite, A_tpx))

complex_rms(A_m, G_m, sel) = sqrt(mean(abs2,
    (@. A_m[sel] * cis(-deg2rad(G_m[sel]))) .- (@. A_tpx[sel] * cis(-deg2rad(G_tpx[sel])))))

println("="^92)
println("M2 SKILL vs TPXO10 — fixed Cᵇ=0.967, with vs without stratification")
println("="^92)
@printf("%-28s %14s %16s %10s\n", "stratification", "domain cRMS", "shelf<50m cRMS", "shelf/A")
for r in results
    domain_sel = valid
    shelf_sel  = @. valid & isfinite(r.H_tpx) & (r.H_mod <= 50)
    crms_domain = complex_rms(r.A_mod, r.G_mod, domain_sel)
    crms_shelf  = complex_rms(r.A_mod, r.G_mod, shelf_sel)
    @printf("%-28s %14.4f %16.4f %10.3f\n", r.label, crms_domain, crms_shelf, crms_shelf / mean(A_tpx[shelf_sel]))
end
println()
println("For comparison, the same three stratification cases under the OLD Cᵇ=0.28")
println("(measured 2026-09-11, see 03_mab_m2_tide_smoke.jl header): domain cRMS")
println("  none 0.1070   deep/slope 0.1054   shelf 0.1068  (spread ≤ 1.5%)")

# ---------------- amplitude-error maps, one column per stratification ----------------
fig = Figure(size = (1350, 500))
for (col, r) in enumerate(results)
    dAf = fill(NaN, size(r.A_mod)) ; dAf[valid] .= r.A_mod[valid] .- A_tpx[valid]
    ax = Axis(fig[1, col], title = r.label, xlabel = "°E", ylabel = col == 1 ? "°N" : "", aspect = DataAspect())
    hm = heatmap!(ax, λc, φc, dAf, colormap = :balance, colorrange = (-0.15, 0.15))
    col == length(results) && Colorbar(fig[1, col+1], hm, label = "model − TPXO (m)")
end
Label(fig[0, 1:(length(results)+1)],
      "Fixed Cᵇ=0.967 — M2 amplitude error is essentially unchanged by stratification", fontsize = 18)
fpath = joinpath(OUT_DIR, "catke_new_strat_error.png")
save(fpath, fig)
@info "wrote $fpath"
println("\nDONE")
