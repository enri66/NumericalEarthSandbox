# ==================================================================
# COMPARE CATKE'S Cᵇ FIX (Oceananigans PR #6024) AGAINST THE OLD DEFAULT
# AND THE FIXED-ν_z PROBE, ON THE M2 SMOKE TEST.
#
# Loads the three `*_results.jls` dumps that 03_mab_m2_tide_smoke.jl writes
# (one run each, MAB_CB=0.28 / MAB_CB=0.967 / MAB_CLOSURE=nuz) and builds:
#   - a cotidal comparison: amplitude, phase and amplitude-error maps, one
#     column per closure, plus TPXO
#   - a gauge time series overlay: all three model curves against NOAA and
#     TPXO at the open-water MAB gauges
#   - a skill summary table (domain and by depth class)
#
# Run:  julia --project=. scripts/compare_catke.jl
# ==================================================================

using Serialization
using Printf
using Statistics
using Oceananigans.Units: days
using CairoMakie
using Dates

const OUT_DIR = get(ENV, "MAB_OUT", joinpath(homedir(), "Data", "mab_tides", "catke_comparison"))

runs = (
    ("catke_old",  "CATKE, old Cᵇ=0.28"),
    ("catke_new",  "CATKE, fixed Cᵇ=0.967 (PR #6024)"),
    ("nuz_fixed",  "fixed high ν_z = 0.1 m²/s"),
)

results = NamedTuple[]
for (tag, label) in runs
    path = joinpath(OUT_DIR, "$(tag)_results.jls")
    isfile(path) || error("missing $path — run 03_mab_m2_tide_smoke.jl with MAB_TAG=$tag first")
    push!(results, merge(deserialize(path), (; label)))
end

# TPXO/NOAA truth and the analysis grid are identical across runs (same box, resolution,
# start date and constituent) — take them from the first result.
r1 = results[1]
λc, φc, A_tpx, G_tpx, valid, wet = r1.λc, r1.φc, r1.A_tpx, r1.G_tpx, r1.valid, r1.wet
maskland(A) = replace(x -> isfinite(x) ? x : NaN, A)
amax = maximum(filter(isfinite, A_tpx))

# ---------------- skill summary ----------------
complex_rms(A_m, G_m, sel) = sqrt(mean(abs2,
    (@. A_m[sel] * cis(-deg2rad(G_m[sel]))) .- (@. A_tpx[sel] * cis(-deg2rad(G_tpx[sel])))))

println("="^96)
println("M2 SKILL vs TPXO10 — CATKE Cᵇ comparison")
println("="^96)
@printf("%-34s %8s %14s %14s %10s\n", "run", "Cᵇ", "domain cRMS", "shelf<50m cRMS", "shelf/A")
for r in results
    domain_sel = valid
    shelf_sel  = @. valid & isfinite(r.H_tpx) & (r.H_mod <= 50)
    crms_domain = complex_rms(r.A_mod, r.G_mod, domain_sel)
    crms_shelf  = complex_rms(r.A_mod, r.G_mod, shelf_sel)
    shelf_ratio = crms_shelf / mean(A_tpx[shelf_sel])
    @printf("%-34s %8s %14.4f %14.4f %10.3f\n",
            r.label, r.Cᵇ === missing ? "n/a" : @sprintf("%.3f", r.Cᵇ), crms_domain, crms_shelf, shelf_ratio)
end

# ---------------- cotidal comparison: amplitude, phase, error ----------------
fig = Figure(size = (1650, 1150))

for (col, r) in enumerate(results)
    ax = Axis(fig[1, col], title = r.label, xlabel = "°E", ylabel = col == 1 ? "°N" : "", aspect = DataAspect())
    hm = heatmap!(ax, λc, φc, maskland(r.A_mod), colormap = :viridis, colorrange = (0, amax))
    col == length(results) && Colorbar(fig[1, col+1], hm, label = "M2 amplitude (m)")

    ax2 = Axis(fig[2, col], xlabel = "°E", ylabel = col == 1 ? "°N" : "", aspect = DataAspect())
    hm2 = heatmap!(ax2, λc, φc, maskland(r.G_mod), colormap = :twilight, colorrange = (0, 360))
    col == length(results) && Colorbar(fig[2, col+1], hm2, label = "M2 phase (°)")

    dAf = fill(NaN, size(r.A_mod)) ; dAf[valid] .= r.A_mod[valid] .- A_tpx[valid]
    ax3 = Axis(fig[3, col], xlabel = "°E", ylabel = col == 1 ? "°N" : "", aspect = DataAspect())
    hm3 = heatmap!(ax3, λc, φc, dAf, colormap = :balance, colorrange = (-0.15, 0.15))
    col == length(results) && Colorbar(fig[3, col+1], hm3, label = "model − TPXO (m)")
end

axT = Axis(fig[1, length(results)+2], title = "TPXO10 (truth)", xlabel = "°E", aspect = DataAspect())
heatmap!(axT, λc, φc, maskland(A_tpx), colormap = :viridis, colorrange = (0, amax))
axT2 = Axis(fig[2, length(results)+2], xlabel = "°E", aspect = DataAspect())
heatmap!(axT2, λc, φc, maskland(G_tpx), colormap = :twilight, colorrange = (0, 360))

Label(fig[0, 1:(length(results)+2)], "M2 cotidal comparison — CATKE Cᵇ fix vs old default vs fixed ν_z", fontsize = 20)
cpath = joinpath(OUT_DIR, "catke_comparison_cotidal.png")
save(cpath, fig)
@info "wrote $cpath"

# ---------------- co-phase overlay: all three + TPXO on one map ----------------
signed_phase_diff(G, θ) = @. mod(G - θ + 180, 360) - 180
function cophase_field(G, θ)
    d = signed_phase_diff(G, θ)
    return map(x -> (isfinite(x) && abs(x) < 90) ? x : NaN, d)
end

colors = (:steelblue, :seagreen, :darkorange)
figo = Figure(size = (900, 750))
axo = Axis(figo[1, 1], title = "M2 co-phase, every 30° — model runs vs TPXO (grey, dashed)",
          xlabel = "°E", ylabel = "°N", aspect = DataAspect())
heatmap!(axo, λc, φc, maskland(A_tpx), colormap = (:grays, 0.3), colorrange = (0, amax))
θs = 0:30:330
for θ in θs
    contour!(axo, λc, φc, cophase_field(G_tpx, θ); levels = [0.0], color = (:black, 0.6),
            linewidth = 1.2, linestyle = :dash)
end
for (r, c) in zip(results, colors)
    for θ in θs
        contour!(axo, λc, φc, cophase_field(r.G_mod, θ); levels = [0.0], color = (c, 0.9), linewidth = 1.3)
    end
end
elems = [LineElement(color = c, linewidth = 2) for c in colors]
push!(elems, LineElement(color = :black, linewidth = 2, linestyle = :dash))
Legend(figo[1, 2], elems, [[r.label for r in results]; "TPXO10"]; framevisible = false)
opath = joinpath(OUT_DIR, "catke_comparison_cophase.png")
save(opath, figo)
@info "wrote $opath"

# ---------------- gauge time series overlay ----------------
# Open-water gauges only — the fair comparison (see noaa_harcon_mab.jl).
show_names = ("Duck, NC", "Atlantic City, NJ", "Montauk, NY", "Sandy Hook, NJ", "Sewells Point, VA")
include(joinpath(@__DIR__, "tidal_harmonics.jl"))  # reconstruct, reconstruction_parameters

harmonics_for_recon = TidalHarmonics(r1.start_date; constituents = (:M2,))
p_rec = reconstruction_parameters(harmonics_for_recon)

hits = [(name, i, j, A_n, G_n) for (name, i, j, A_n, G_n) in r1.gauge_hits if name in show_names]
if !isempty(hits)
    ng = length(hits)
    figg = Figure(size = (1500, 280 * ng))
    for (row, (name, i, j, A_n, G_n)) in enumerate(hits)
        times1 = r1.times
        tshow = findall(t -> t > last(times1) - 4days, times1)
        tt = times1[tshow] ./ days

        noaa = [reconstruct([A_n], [deg2rad(G_n)], times1[n], p_rec) for n in tshow]
        tpxo = [reconstruct([A_tpx[i,j]], [deg2rad(G_tpx[i,j])], times1[n], p_rec) for n in tshow]

        ax = Axis(figg[row, 1],
                  title = name,
                  xlabel = row == ng ? "days since $(Dates.format(r1.start_date, "yyyy-mm-dd"))" : "",
                  ylabel = "η (m)")
        for (r, c) in zip(results, colors)
            model = [r.η_store[n][i, j] for n in tshow]
            lines!(ax, tt, model, color = c, linewidth = 2.2, label = r.label)
        end
        lines!(ax, tt, noaa, color = :black, linestyle = :dash, linewidth = 1.8, label = "NOAA M2")
        lines!(ax, tt, tpxo, color = :gray40, linestyle = :dot, linewidth = 1.8, label = "TPXO M2")
        row == 1 && axislegend(ax, position = :rt, framevisible = false, labelsize = 11)
    end
    Label(figg[0, 1], "M2 sea level at MAB tide gauges — CATKE Cᵇ comparison", fontsize = 18)
    gpath = joinpath(OUT_DIR, "catke_comparison_gauges.png")
    save(gpath, figg)
    @info "wrote $gpath"
end

println("\nDONE")
