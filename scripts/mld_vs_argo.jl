# Mixed-layer and seasonal-thermocline depths of a run and of GLORYS against Argo profiles in the MAB box. Each Argo
# profile is compared with the model's de-tided daily column and GLORYS's column at the nearest model cell and the
# nearest daily frame, and all three use the same definitions (column_depths in mab_analysis_common.jl). Argo levels
# are used only with good quality flags (1 or 2); delayed-mode and adjusted profiles use the adjusted values.
# The profiles come from the Ifremer ERDDAP ArgoFloats dataset as monthly CSV files (one header row and one units row).
# Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/sponge_year_mpi2/spy ARGO_DIR=/t0/workdir/enrique/Data/Argo/mab \
#       julia --project=. scripts/mld_vs_argo.jl
# Argo's in-situ temperature and practical salinity go into the same σ₀ as the model's and GLORYS's fields, which
# changes the 0.03 kg/m³ mixed-layer threshold by far less than its own uncertainty in the upper few hundred metres.
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const TAG      = get(ENV, "MAB_TAG", "mab_H90")
const PREFIX   = run_prefix(TAG)
const ARGO_DIR = get(ENV, "ARGO_DIR", joinpath(homedir(), "Data", "Argo", "mab"))

# ---------------- Argo ----------------
good(flag) = flag in ("1", "2")
parsefloat(s) = (x = tryparse(Float64, s); isnothing(x) ? NaN : x)

function read_argo(dir)
    profiles = Dict{Tuple{String, String, String}, Any}()
    for f in sort(filter(endswith(".csv"), readdir(dir; join = true)))
        lines = readlines(f)
        length(lines) < 3 && continue
        header = split(lines[1], ",")
        col = Dict(name => i for (i, name) in enumerate(header))
        for line in lines[3:end]
            v = split(line, ",")
            good(v[col["position_qc"]]) || continue
            adjusted = v[col["data_mode"]] in ("D", "A")
            p, pq = adjusted ? (v[col["pres_adjusted"]], v[col["pres_adjusted_qc"]]) : (v[col["pres"]], v[col["pres_qc"]])
            T, Tq = adjusted ? (v[col["temp_adjusted"]], v[col["temp_adjusted_qc"]]) : (v[col["temp"]], v[col["temp_qc"]])
            S, Sq = adjusted ? (v[col["psal_adjusted"]], v[col["psal_adjusted_qc"]]) : (v[col["psal"]], v[col["psal_qc"]])
            (good(pq) && good(Tq) && good(Sq)) || continue
            key = (v[col["platform_number"]], v[col["cycle_number"]], v[col["direction"]])
            prof = get!(profiles, key) do
                (time = DateTime(v[col["time"]][1:19]), lat = parsefloat(v[col["latitude"]]), lon = parsefloat(v[col["longitude"]]),
                 p = Float64[], T = Float64[], S = Float64[])
            end
            push!(prof.p, parsefloat(p)); push!(prof.T, parsefloat(T)); push!(prof.S, parsefloat(S))
        end
    end
    out = []
    for (key, prof) in profiles
        o = sortperm(prof.p)
        push!(out, (; key, prof.time, prof.lat, prof.lon, p = prof.p[o], T = prof.T[o], S = prof.S[o]))
    end
    return sort(out; by = x -> x.time)
end

argo = read_argo(ARGO_DIR)
@printf("%d Argo profiles with good data, %s to %s\n", length(argo), first(argo).time, last(argo).time)

# ---------------- model ----------------
vol = PREFIX * "_volume_daily.jld2"
TV = open_series(vol, "T"; backend = OnDisk()); SV = open_series(vol, "S"; backend = OnDisk())
grid = TV.grid; ug = grid.underlying_grid
Nx, Ny, Nz = size(ug)
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center())); zc = collect(znodes(ug, Center()))
B = grid.immersed_boundary.bottom_height
bh = [B[i, j, 1] for i in 1:Nx, j in 1:Ny]
wet = bh .< 0
depth_m = reverse(-zc)
tdays = TV.times ./ 86400
region(i, j) = -bh[i, j] < 200 ? 1 : -bh[i, j] < 1000 ? 2 : 3
region_names = ("shelf (<200 m)", "slope (200-1000 m)", "deep (>1000 m)")

# ---------------- compare ----------------
rows = []
# The profiles are in time order, so only the current day's model and GLORYS fields are kept
current = Dict{Symbol, Any}(:n => 0, :day => -1)
for prof in argo
    (λ[1] <= prof.lon <= λ[end] && φ[1] <= prof.lat <= φ[end]) || continue
    i = argmin(abs.(λ .- prof.lon)); j = argmin(abs.(φ .- prof.lat))
    wet[i, j] || continue
    d = Dates.value(prof.time - start_date) / 86_400_000          # days since the start, as a real number
    n = argmin(abs.(tdays .- d))
    abs(tdays[n] - d) <= 1 || continue                            # outside the run's daily frames
    a_mld, a_tcl = column_depths(prof.p, prof.T, prof.S)
    isfinite(a_mld) || continue
    if current[:n] != n
        current[:T] = Array(interior(TV[n])); current[:S] = Array(interior(SV[n])); current[:n] = n
    end
    T = current[:T][i, j, :]; S = current[:S][i, j, :]
    k = [zc[k] > bh[i, j] for k in 1:Nz]
    T[.!k] .= NaN; S[.!k] .= NaN
    m_mld, m_tcl = column_depths(depth_m, reverse(T), reverse(S))
    day = round(Int, tdays[n])
    if current[:day] != day
        lon, lat, dep, GT = glorys_at("thetao", day); _, _, _, GS = glorys_at("so", day)
        current[:glorys] = (; lon, lat, dep, GT, GS); current[:day] = day
    end
    g = current[:glorys]
    gT = glorys_profile(g.GT, g.lon, g.lat, g.dep, prof.lon, prof.lat, -bh[i, j])
    gS = glorys_profile(g.GS, g.lon, g.lat, g.dep, prof.lon, prof.lat, -bh[i, j])
    g_mld, g_tcl = column_depths(g.dep, gT, gS)
    push!(rows, (; prof.time, prof.lat, prof.lon, region = region(i, j), month = Dates.format(prof.time, "yyyy-mm"),
                 argo = a_mld, model = m_mld, glorys = g_mld, argo_tcl = a_tcl, model_tcl = m_tcl, glorys_tcl = g_tcl))
end
@printf("%d profiles inside the run's grid and period\n\n", length(rows))

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
