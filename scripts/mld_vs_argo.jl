# Mixed-layer and seasonal-thermocline depths of a run and of GLORYS against Argo profiles in the MAB box. Each Argo
# profile is compared with the model's de-tided daily column and GLORYS's column at the nearest model cell and the
# nearest daily frame (argo_common.jl), and all three use the same definitions (column_depths in mab_analysis_common.jl).
# Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/sponge_year_mpi2/spy ARGO_DIR=/t0/workdir/enrique/Data/Argo/mab \
#       julia --project=. scripts/mld_vs_argo.jl
# Argo's in-situ temperature and practical salinity go into the same σ₀ as the model's and GLORYS's fields, which
# changes the 0.03 kg/m³ mixed-layer threshold by far less than its own uncertainty in the upper few hundred metres.
using CairoMakie
include(joinpath(@__DIR__, "argo_common.jl"))

const TAG      = get(ENV, "MAB_TAG", "mab_H90")
const PREFIX   = run_prefix(TAG)
const ARGO_DIR = get(ENV, "ARGO_DIR", joinpath(homedir(), "Data", "Argo", "mab"))

rows = []
for m in matched_argo_profiles(PREFIX, ARGO_DIR)
    a_mld, a_tcl = column_depths(m.argo...)
    isfinite(a_mld) || continue
    m_mld, m_tcl = column_depths(m.model...)
    g_mld, g_tcl = column_depths(m.glorys...)
    push!(rows, (; m.time, m.lat, m.lon, m.region, m.month,
                 argo = a_mld, model = m_mld, glorys = g_mld, argo_tcl = a_tcl, model_tcl = m_tcl, glorys_tcl = g_tcl))
end
@printf("%d profiles with an Argo mixed-layer depth\n\n", length(rows))

stat(x, y) = (f = isfinite.(x) .& isfinite.(y); n = count(f);
              n == 0 ? (NaN, NaN, NaN, 0) : (median(x[f] .- y[f]), mean(x[f] .- y[f]), sqrt(mean((x[f] .- y[f]) .^ 2)), n))
col(rs, name) = [getproperty(r, name) for r in rs]

println("Mixed-layer depth bias against Argo (m): median / mean / rms (profiles)")
@printf("%-20s %-28s %-28s %8s\n", "group", "model − Argo", "GLORYS − Argo", "Argo median")
function report(label, rs)
    isempty(rs) && return
    m = stat(col(rs, :model), col(rs, :argo)); g = stat(col(rs, :glorys), col(rs, :argo))
    @printf("%-20s %+6.1f / %+6.1f / %5.1f (%3d)   %+6.1f / %+6.1f / %5.1f (%3d)   %6.1f\n", label, m[1:3]..., m[4], g[1:3]..., g[4],
            median(col(rs, :argo)))
end
report("all", rows)
for r in 1:3; report(region_names[r], filter(x -> x.region == r, rows)); end
println()
for month in sort(unique(col(rows, :month))); report(month, filter(x -> x.month == month, rows)); end

println("\nSeasonal-thermocline depth bias against Argo (m): median / mean / rms (profiles)")
let rs = rows
    m = stat(col(rs, :model_tcl), col(rs, :argo_tcl)); g = stat(col(rs, :glorys_tcl), col(rs, :argo_tcl))
    @printf("%-20s %+6.1f / %+6.1f / %5.1f (%3d)   %+6.1f / %+6.1f / %5.1f (%3d)\n", "all", m[1:3]..., m[4], g[1:3]..., g[4])
end

# ---------------- figure ----------------
fig = Figure(size = (1700, 600), fontsize = 15)
Label(fig[0, 1:3], "$(basename(PREFIX)): mixed-layer depth against Argo ($(length(rows)) profiles)", fontsize = 19)
for (c, (name, label)) in enumerate(((:model, "model"), (:glorys, "GLORYS")))
    ax = Axis(fig[1, c], xlabel = "Argo MLD (m)", ylabel = "$(label) MLD (m)", title = "$(label) vs Argo",
              xscale = log10, yscale = log10)
    for r in 1:3
        rs = filter(x -> x.region == r && isfinite(getproperty(x, name)), rows)
        isempty(rs) || scatter!(ax, col(rs, :argo), col(rs, name); markersize = 6, label = region_names[r])
    end
    lines!(ax, [5, 500], [5, 500]; color = :black)
    c == 1 && axislegend(ax; position = :lt)
end
ax = Axis(fig[1, 3], xlabel = "month", ylabel = "median bias (m)", title = "monthly median MLD bias against Argo")
months = sort(unique(col(rows, :month)))
for (name, label) in ((:model, "model"), (:glorys, "GLORYS"))
    b = [stat(col(filter(x -> x.month == m, rows), name), col(filter(x -> x.month == m, rows), :argo))[1] for m in months]
    scatterlines!(ax, 1:length(months), b; label)
end
hlines!(ax, 0; color = :gray)
ax.xticks = (1:length(months), [m[3:end] for m in months]); ax.xticklabelrotation = π / 4
axislegend(ax; position = :lt)
out = PREFIX * "_mld_vs_argo.png"
save(out, fig; px_per_unit = 1.2)
println("saved ", out)
