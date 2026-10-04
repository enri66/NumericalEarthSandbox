# Where in the water column is the near-inertial energy, in the runs and at the OOI Pioneer current profilers? If the
# model's missing upper-ocean near-inertial energy shows up deeper, it is radiated down; if it is missing at every
# depth, it is lost. For each mooring and source (ADCP, and every run in MAB_TAGS with MAB_MOORINGS=pioneer columns),
# the near-inertial velocity (16-22 h band-pass at each depth, ooi_common.jl) gives KE(z, t), smoothed over one
# inertial period, on a common depth grid from the top of the usable ADCP range to the shallowest of the ADCP's
# range reached in three quarters of the hours and the runs' cell depths (minus three bins). Reported, per layer:
#   the depth-integrated near-inertial energy ∫ KE dz (m³/s², i.e. J/m² per unit density), record mean and in storm
#   windows, and model / ADCP;
#   the share of the energy below the top layer;
#   after the 17 October storm, the time each layer's energy peaks, relative to the top layer's peak (hours).
# The figure shows depth-time sections of log₁₀ KE.
# Usage:
#   MAB_TAGS=/t0/.../catke_fall/moor,/t0/.../res_test/r12,/t0/.../res_test/r24 OOI_DIR=/t0/workdir/enrique/Data/OOI/pioneer \
#       julia --project=. scripts/ni_budget.jl
using CairoMakie
include(joinpath(@__DIR__, "ooi_common.jl"))

const TAGS   = filter(!isempty, split(get(ENV, "MAB_TAGS", ""), ","))
const AXIS   = collect(DateTime(2019, 8, 30):Hour(1):DateTime(2019, 10, 27))
const STORMS = (("Dorian", DateTime(2019, 9, 6), DateTime(2019, 9, 13)),
                ("mid-Sep", DateTime(2019, 9, 14), DateTime(2019, 9, 21)),
                ("17 Oct", DateTime(2019, 10, 17), DateTime(2019, 10, 24)))
const LAYERS = Dict(8.0 => ((40, 100), (100, 200), (200, 1000)),     # slope moorings (8 m grid)
                    4.0 => ((12, 40), (40, 80), (80, 1000)))         # shelf moorings (4 m grid)
const OUTDIR = get(ENV, "NI_OUTDIR", OOI_DIR)

function ni_ke(src, grid)
    ti = Dict(t => n for (n, t) in enumerate(src.times))
    take(A) = [haskey(ti, t) ? A[m, ti[t]] : NaN for m in axes(A, 1), t in AXIS]
    U = take(regrid(src.depths, src.U, grid)); V = take(regrid(src.depths, src.V, grid))
    KE = fill(NaN, length(grid), length(AXIS))
    w = round(Int, T_INERTIAL ÷ 2)
    for m in eachindex(grid)
        ub = bandpass(U[m, :]); vb = bandpass(V[m, :])
        k = 0.5 .* (ub .^ 2 .+ vb .^ 2)
        KE[m, :] = [nanmean(k[max(1, n - w):min(end, n + w)]) for n in eachindex(k)]
    end
    return KE
end

layer_rows(grid, (a, b)) = findall(z -> a <= z < b, grid)
integral(KE, rows, Δz) = [sum(x -> isfinite(x) ? x : 0.0, KE[rows, n]) * Δz for n in axes(KE, 2)]
window(a, b) = findall(t -> a <= t < b, AXIS)

results = Dict{String, Any}()
for name in ("OSSM", "PMUO", "PMCO", "CNSM")
    dataset, Δz, ztop = MOORINGS[name]
    adcp = read_adcp(dataset)
    runs = [(basename(String(t)), read_model(String(t), name)) for t in TAGS]
    runs = filter(r -> !isnothing(r[2]), runs)
    # The bins' depths move with the mooring's tilt and the tide, so take the deepest valid bin of each hour and keep
    # the depth that three quarters of the hours reach
    deepest = [maximum(adcp.depths[isfinite.(adcp.U[:, n])]; init = 0.0) for n in axes(adcp.U, 2)]
    deep_adcp = quantile(filter(>(0), deepest), 0.25)
    floors = [maximum(r.depths[[any(x -> isfinite(x) && x != 0, r.U[k, :]) for k in eachindex(r.depths)]]) for (_, r) in runs]
    zmax = min(deep_adcp, minimum(floors) - 3Δz)
    grid = collect(ztop:Δz:zmax)
    sources = [("ADCP", adcp); runs]
    KE = Dict(label => ni_ke(src, grid) for (label, src) in sources)
    layers = [(a, min(b, grid[end] + Δz)) for (a, b) in LAYERS[Δz] if a < grid[end]]
    results[name] = (; grid, KE, layers, labels = first.(sources))

    @printf("\n== %s: common grid %g-%g m every %g m (ADCP well sampled to %g m; run columns to %s m)\n", name, grid[1], grid[end],
            Δz, deep_adcp, join(round.(Int, floors), "/"))
    E = Dict(label => [integral(KE[label], layer_rows(grid, l), Δz) for l in layers] for label in first.(sources))
    println("   depth-integrated near-inertial energy (10⁻³ m³/s²) by layer: record mean / Dorian / mid-Sep / 17 Oct")
    for (i, l) in enumerate(layers)
        @printf("   %4d-%-4d m", l[1], l[2])
        for label in first.(sources)
            e = E[label][i]
            @printf("  %s %5.2f/%5.2f/%5.2f/%5.2f", label, 1e3 * mean(e), [1e3 * mean(e[window(a, b)]) for (_, a, b) in STORMS]...)
        end
        println()
    end
    println("   model / ADCP, record mean by layer, and share of the energy below $(layers[1][2]) m:")
    total(label) = sum(mean.(E[label]))
    @printf("   %-6s (reference)   share below %.0f%%\n", "ADCP", 100 * (1 - mean(E["ADCP"][1]) / total("ADCP")))
    for (label, _) in runs
        ratios = [mean(E[label][i]) / mean(E["ADCP"][i]) for i in eachindex(layers)]
        @printf("   %-6s %s   share below %.0f%%   (whole column %.2f of ADCP)\n", label, join([@sprintf("%5.2f", r) for r in ratios], " "),
                100 * (1 - mean(E[label][1]) / total(label)), total(label) / total("ADCP"))
    end
    # Peak timing after the 17 October storm
    w = window(DateTime(2019, 10, 16), DateTime(2019, 10, 25))
    print("   17 Oct: hours from the top layer's peak to each layer's peak:")
    for label in first.(sources)
        t = [AXIS[w[argmax(E[label][i][w])]] for i in eachindex(layers)]
        @printf("  %s %s", label, join([string(round(Int, Dates.value(t[i] - t[1]) / 3_600_000)) for i in eachindex(layers)], "/"))
    end
    println()
end

# ---------------- figure ----------------
fignames = ("OSSM", "PMUO", "CNSM")
ncols = maximum(length(results[n].labels) for n in fignames)
fig = Figure(size = (380 * ncols + 120, 330 * length(fignames) + 80), fontsize = 13)
Label(fig[0, 1:ncols], "near-inertial kinetic energy, log₁₀ KE (J/kg): Pioneer ADCPs and model runs", fontsize = 17)
td = [Dates.value(t - AXIS[1]) / 86_400_000 for t in AXIS]
hm = nothing
for (row, name) in enumerate(fignames)
    r = results[name]
    for (col, label) in enumerate(r.labels)
        ax = Axis(fig[row, col], title = "$name: $label", yreversed = true, ylabel = col == 1 ? "depth (m)" : "",
                  xlabel = row == length(fignames) ? "days since 2019-08-30" : "")
        global hm = heatmap!(ax, td, r.grid, log10.(max.(r.KE[label], 1e-7))'; colormap = :thermal, colorrange = (-5.5, -1.5))
        for (a, b) in r.layers[2:end]; hlines!(ax, a; color = :white, linestyle = :dot, linewidth = 1); end
    end
end
Colorbar(fig[1:length(fignames), ncols + 1], hm; label = "log₁₀ KE")
out = joinpath(OUTDIR, "ni_budget.png")
save(out, fig; px_per_unit = 1.2)
println("\nsaved ", out)
