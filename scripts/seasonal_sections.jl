# Temperature sections of a run, one row per season, with the mixed-layer and seasonal-thermocline depths of the model
# and of GLORYS drawn on. Reads the de-tided daily full-depth output (`<tag>_volume_daily.jld2`), including script 05's
# per-rank output. Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/sponge_year_mpi2/spy julia --project=. scripts/seasonal_sections.jl
# MAB_DATES:    one date per row (default mid-April, July and October 2019 and mid-January 2020)
# MAB_SECTIONS: "lat:<φ>" (zonal) or "lon:<λ>" (meridional) sections, comma separated
# MAB_ZMAX:     depth shown in the sections (m)
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const TAG    = get(ENV, "MAB_TAG", "mab_H90")
const PREFIX = run_prefix(TAG)
const ZMAX   = parse(Float64, get(ENV, "MAB_ZMAX", "300"))
const DATES  = Date.(split(get(ENV, "MAB_DATES", "2019-04-15,2019-07-15,2019-10-15,2020-01-15"), ","))
const SECTIONS = [(Symbol(split(s, ":")[1]), parse(Float64, split(s, ":")[2]))
                  for s in split(get(ENV, "MAB_SECTIONS", "lat:35.0,lat:38.0,lat:40.5,lon:-70.5"), ",")]

vol = PREFIX * "_volume_daily.jld2"
TV = open_series(vol, "T"; backend = OnDisk()); SV = open_series(vol, "S"; backend = OnDisk())
grid = TV.grid; ug = grid.underlying_grid
Nx, Ny, Nz = size(ug)
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center())); zc = collect(znodes(ug, Center()))
B = grid.immersed_boundary.bottom_height
bh = [B[i, j, 1] for i in 1:Nx, j in 1:Ny]
wet = bh .< 0
wet3 = [zc[k] > bh[i, j] for i in 1:Nx, j in 1:Ny, k in 1:Nz]
tdays = TV.times ./ 86400
depth_m = reverse(-zc)                                  # model depths, top first

season(d) = month(d) in (3, 4, 5) ? "spring" : month(d) in (6, 7, 8) ? "summer" : month(d) in (9, 10, 11) ? "fall" : "winter"

function section_columns(kind, v)
    if kind === :lat
        j = argmin(abs.(φ .- v))
        return [(i, j) for i in 1:Nx], λ, "longitude", @sprintf("%.2f°N", φ[j])
    else
        i = argmin(abs.(λ .- v))
        return [(i, j) for j in 1:Ny], φ, "latitude", @sprintf("%.2f°E", λ[i])
    end
end

fig = Figure(size = (560 * length(SECTIONS) + 120, 330 * length(DATES) + 120), fontsize = 15)
Label(fig[0, 1:length(SECTIONS)], "$(basename(PREFIX)): model temperature by season, with mixed-layer base (solid) and " *
      "seasonal thermocline (dashed), model (white) and GLORYS (cyan)", fontsize = 19)
hm = nothing
println("median depth along each section (m):          mixed layer            seasonal thermocline")
println("date        section                           model   GLORYS          model   GLORYS")
for (row, date) in enumerate(DATES)
    dwant = Dates.value(Date(date) - Date(start_date))
    n = argmin(abs.(tdays .- dwant)); d = tdays[n]
    T = Array(interior(TV[n])); S = Array(interior(SV[n]))
    T[.!wet3] .= NaN; S[.!wet3] .= NaN
    lon, lat, dep, GT = glorys_at("thetao", d)
    _, _, _, GS = glorys_at("so", d)
    for (col, (kind, v)) in enumerate(SECTIONS)
        cols, x, xl, lbl = section_columns(kind, v)
        Tm = [T[i, j, k] for (i, j) in cols, k in 1:Nz]
        depths = map(cols) do (i, j)
            wet[i, j] || return (NaN, NaN, NaN, NaN)
            m = column_depths(depth_m, reverse(T[i, j, :]), reverse(S[i, j, :]))
            g = column_depths(dep, glorys_profile(GT, lon, lat, dep, λ[i], φ[j], -bh[i, j]),
                                   glorys_profile(GS, lon, lat, dep, λ[i], φ[j], -bh[i, j]))
            return (m..., g...)
        end
        mml, mtl, gml, gtl = ([p[q] for p in depths] for q in 1:4)
        med(a) = (f = filter(isfinite, a); isempty(f) ? NaN : median(f))
        @printf("%s  %-6s %-10s %-10s        %6.1f   %6.1f         %6.1f   %6.1f\n", date, String(kind), lbl, season(date),
                med(mml), med(gml), med(mtl), med(gtl))

        ax = Axis(fig[row, col], title = @sprintf("%s, %s (%s)", lbl, Dates.format(Date(date), "d u yyyy"), season(date)),
                  xlabel = row == length(DATES) ? xl : "", ylabel = col == 1 ? "z (m)" : "")
        p = sortperm(zc)
        global hm = heatmap!(ax, x, zc[p], Tm[:, p]; colormap = :thermal, colorrange = (4, 30), nan_color = :gray85)
        lines!(ax, x, -mml; color = :white, linewidth = 2)
        lines!(ax, x, -mtl; color = :white, linewidth = 2, linestyle = :dash)
        lines!(ax, x, -gml; color = :cyan, linewidth = 1.5)
        lines!(ax, x, -gtl; color = :cyan, linewidth = 1.5, linestyle = :dash)
        lines!(ax, x, [wet[i, j] ? bh[i, j] : 0.0 for (i, j) in cols]; color = :black, linewidth = 1)
        ylims!(ax, -ZMAX, 0)
    end
end
Colorbar(fig[1:length(DATES), length(SECTIONS) + 1], hm; label = "T (°C)")
out = PREFIX * "_seasonal_sections.png"
save(out, fig; px_per_unit = 1.2)
println("saved ", out)
