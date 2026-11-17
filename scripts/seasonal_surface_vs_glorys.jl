# Seasonal evolution of the surface against GLORYS for a long run: every de-tided daily frame of SST, SSH anomaly and
# surface speed against the mean of the two GLORYS daily means around it (as surface_vs_glorys.jl), as daily time series
# of bias, rms difference and pattern correlation over the shelf, the slope and the deep water, and monthly mean maps of
# model − GLORYS for SST and SSH anomaly. SSH is compared as each one's anomaly from its own mean over the wet cells.
# Usage (a run in segments):
#   MAB_TAGS=/t0/.../mabyear_s01,/t0/.../mabyear_s02,... MAB_OUT=/t0/.../analysis/year/mabyear MAB_START_DATE=2019-04-01 \
#   MAB_DATA_DIR=/t0/workdir/enrique/Data/NumericalEarth julia --project=<run env> scripts/seasonal_surface_vs_glorys.jl
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const PREFIXES = run_prefixes()
const OUT = get(ENV, "MAB_OUT", PREFIXES[1])

# ---------------- model ----------------
grid = global_grid(PREFIXES[1] * "_surface_daily.jld2"); ug = grid.underlying_grid
Nx, Ny, _ = size(ug)
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center()))
B = grid.immersed_boundary.bottom_height
bh = [B[i, j, 1] for i in 1:Nx, j in 1:Ny]
wet = bh .< 0

Tm, times = joined_surface_series(PREFIXES, "_surface_daily", "T")
um, _ = joined_surface_series(PREFIXES, "_surface_daily", "u")
vm, _ = joined_surface_series(PREFIXES, "_surface_daily", "v")
ηm, ηtimes = joined_surface_series(PREFIXES, "_eta_daily", "η")
tdays = times ./ 86400
@printf("%d segments, %d daily frames, days %.1f to %.1f\n", length(PREFIXES), length(tdays), tdays[1], tdays[end])

centre_u(u) = (u[1:Nx, :] .+ u[2:Nx+1, :]) ./ 2
centre_v(v) = (v[:, 1:Ny] .+ v[:, 2:Ny+1]) ./ 2
mask(X) = ifelse.(wet, X, NaN)
anomaly(X) = X .- mean(filter(isfinite, X))

# ---------------- GLORYS on the model's wet cells ----------------
const GFILES = Dict{Tuple{String, Date}, String}()
for f in readdir(GDIR)
    m = match(r"^(thetao|zos|uo|vo)_GLORYSDaily_(\d{4}-\d{2}-\d{2})T.*_-76\.0_-64\.0_34\.0_42\.0\.nc$", f)
    isnothing(m) || (GFILES[(String(m[1]), Date(m[2]))] = joinpath(GDIR, f))
end

function read_glorys_surface(var, date)
    ds = NCD.Dataset(GFILES[(var, Date(date))])
    lon = Float64.(ds["longitude"][:]); lat = Float64.(ds["latitude"][:])
    raw = ds[var]
    A = ndims(raw) == 4 ? raw[:, :, 1, 1] : raw[:, :, 1]
    close(ds)
    return lon, lat, [ismissing(x) ? NaN : Float64(x) for x in A]
end

# bilinear weights from the GLORYS grid to each wet model cell, computed once
glon, glat, _ = read_glorys_surface("thetao", start_date)
cells = [(i, j) for i in 1:Nx, j in 1:Ny if wet[i, j]]
stencil = map(cells) do (i, j)
    x, y = λ[i], φ[j]
    a = clamp(searchsortedlast(glon, x), 1, length(glon) - 1); b = clamp(searchsortedlast(glat, y), 1, length(glat) - 1)
    p = (x - glon[a]) / (glon[a+1] - glon[a]); q = (y - glat[b]) / (glat[b+1] - glat[b])
    (a, b, ((1 - p) * (1 - q), p * (1 - q), (1 - p) * q, p * q))
end
function on_model(A)
    out = fill(NaN, Nx, Ny)
    for (c, (i, j)) in enumerate(cells)
        a, b, w = stencil[c]
        v = (A[a, b], A[a+1, b], A[a, b+1], A[a+1, b+1])
        num = 0.0; den = 0.0
        for n in 1:4
            isfinite(v[n]) && (num += w[n] * v[n]; den += w[n])
        end
        den > 0 && (out[i, j] = num / den)
    end
    return out
end

# the mean of the two daily means around the frame's midnight
function glorys_surface_at(var, d)
    day0 = start_date + Day(floor(Int, d)) - Day(1)
    _, _, A0 = read_glorys_surface(var, day0)
    _, _, A1 = read_glorys_surface(var, day0 + Day(1))
    return on_model((A0 .+ A1) ./ 2)
end

# ---------------- daily statistics and monthly sums ----------------
region(i, j) = -bh[i, j] < 200 ? 1 : -bh[i, j] < 1000 ? 2 : 3
const RNAMES = ("shelf (<200 m)", "slope (200-1000 m)", "deep (>1000 m)")
regmask = [[wet[i, j] && region(i, j) == r for i in 1:Nx, j in 1:Ny] for r in 1:3]
function statistics(X, Y, r)
    sel = regmask[r] .& isfinite.(X) .& isfinite.(Y)
    count(sel) < 10 && return (NaN, NaN, NaN)
    x = X[sel]; y = Y[sel]
    return (mean(x .- y), sqrt(mean((x .- y) .^ 2)), cor(x, y))
end

const VARS = ("SST (°C)", "SSH anomaly (m)", "surface speed (m/s)")
stats = Dict((v, r) => fill(NaN, 3, length(tdays)) for v in 1:3, r in 1:3)   # rows: bias, rms, correlation
months = Date[]; msum = Dict{Tuple{Int, Date}, Matrix{Float64}}(); mcount = Dict{Tuple{Int, Date}, Matrix{Int}}()

for n in eachindex(tdays)
    d = tdays[n]
    date = start_date + Second(round(Int, d * 86400))
    nη = argmin(abs.(ηtimes ./ 86400 .- d))
    SST = mask(Tm[:, :, n]); SSH = anomaly(mask(ηm[:, :, nη]))
    speed = mask(sqrt.(centre_u(um[:, :, n]) .^ 2 .+ centre_v(vm[:, :, n]) .^ 2))
    gSST = glorys_surface_at("thetao", d); gSSH = anomaly(glorys_surface_at("zos", d))
    gspeed = sqrt.(glorys_surface_at("uo", d) .^ 2 .+ glorys_surface_at("vo", d) .^ 2)
    for (v, (X, Y)) in enumerate(((SST, gSST), (SSH, gSSH), (speed, gspeed))), r in 1:3
        stats[(v, r)][:, n] .= statistics(X, Y, r)
    end
    mon = Date(year(date), Dates.month(date), 1)
    mon in months || push!(months, mon)
    for (v, X, Y) in ((1, SST, gSST), (2, SSH, gSSH))
        D = X .- Y; ok = isfinite.(D)
        s = get!(msum, (v, mon), zeros(Nx, Ny)); c = get!(mcount, (v, mon), zeros(Int, Nx, Ny))
        s[ok] .+= D[ok]; c[ok] .+= 1
    end
    n % 20 == 0 && @printf("  frame %d of %d (%s)\n", n, length(tdays), Dates.format(date, "yyyy-mm-dd"))
end

# ---------------- tables ----------------
dates = [start_date + Second(round(Int, d * 86400)) for d in tdays]
nanmean(x) = (f = filter(isfinite, x); isempty(f) ? NaN : mean(f))
for (v, name) in enumerate(VARS)
    println("\n", name, ": monthly mean of the daily bias / rms / correlation against GLORYS")
    @printf("  %-9s", "month"); for r in 1:3; @printf(" %-26s", RNAMES[r]); end; println()
    for m in months
        sel = [Date(year(t), Dates.month(t), 1) == m for t in dates]
        @printf("  %-9s", Dates.format(m, "yyyy-mm"))
        for r in 1:3
            s = stats[(v, r)][:, sel]
            @printf(" %+7.3f %7.3f %6.2f     ", nanmean(s[1, :]), nanmean(s[2, :]), nanmean(s[3, :]))
        end
        println()
    end
end

# ---------------- figures ----------------
colors = (:seagreen, :darkorange, :navy)
fig = Figure(size = (1700, 1100), fontsize = 15)
Label(fig[0, 1:3], "Year run vs GLORYS: daily bias, rms difference and pattern correlation by region", fontsize = 19)
for (v, name) in enumerate(VARS), (row, what) in enumerate(("bias", "rms difference", "correlation"))
    ax = Axis(fig[row, v], title = "$name: $what", xticklabelrotation = π / 5)
    for r in 1:3
        lines!(ax, 1:length(dates), stats[(v, r)][row, :]; color = colors[r], label = RNAMES[r])
    end
    row == 1 && hlines!(ax, 0; color = :gray)
    ticks = [findfirst(t -> Date(year(t), Dates.month(t), 1) == m, dates) for m in months]
    ax.xticks = (ticks, Dates.format.(months, "u"))
    (row, v) == (1, 1) && axislegend(ax; position = :lb)
end
save(OUT * "_seasonal_surface_stats.png", fig; px_per_unit = 1.2); println("saved ", OUT * "_seasonal_surface_stats.png")

for (v, name, file, lim) in ((1, "SST, model − GLORYS (°C)", "_seasonal_sst_bias.png", 3.0),
                              (2, "SSH anomaly, model − GLORYS (m)", "_seasonal_ssh_bias.png", 0.4))
    ncol = min(4, length(months)); nrow = cld(length(months), ncol)
    fig = Figure(size = (450 * ncol + 100, 330 * nrow + 80), fontsize = 14)
    Label(fig[0, 1:ncol], "Year run: monthly mean $name", fontsize = 18)
    local hm
    for (k, m) in enumerate(months)
        D = msum[(v, m)] ./ mcount[(v, m)]; D[mcount[(v, m)] .== 0] .= NaN
        ax = Axis(fig[cld(k, ncol), mod1(k, ncol)], title = Dates.format(m, "yyyy-mm") * @sprintf(" (%d days)", maximum(mcount[(v, m)])),
                  aspect = DataAspect())
        hm = heatmap!(ax, λ, φ, D; colormap = :balance, colorrange = (-lim, lim), nan_color = :gray85)
    end
    Colorbar(fig[1:nrow, ncol+1], hm)
    save(OUT * file, fig; px_per_unit = 1.2); println("saved ", OUT * file)
end
