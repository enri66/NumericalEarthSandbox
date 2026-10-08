# Model vs GLORYS at the surface: sea-surface temperature, sea-surface height and surface currents, as maps (model,
# GLORYS, model − GLORYS) with arrows, and a table of bias, rms difference and pattern correlation over the shelf, the
# slope and the deep water. Uses the de-tided daily output (`<tag>_surface_daily.jld2`, `<tag>_eta_daily.jld2`), so the
# model is compared with the mean of the two GLORYS daily means around each midnight (mab_analysis_common.jl). SSH is
# compared as the anomaly from each one's own mean over the wet cells, since the model's η and GLORYS's zos have
# different references. Works with script 05's per-rank output. Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/res_test/hr100_fx3 MAB_START_DATE=2019-08-29 MAB_DAYS=5,10,20 \
#       julia --project=/t0/workdir/enrique/mpi05_ib scripts/surface_vs_glorys.jl
# MAB_TAG:      a run tag in ~/Data/mab_glorys_obc, or an absolute path prefix
# MAB_DAYS:     days to analyse (default: the last daily frame)
# MAB_ARROW_DEG: spacing of the current arrows in degrees (default 0.4)
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const TAG = get(ENV, "MAB_TAG", "mab_H90")
const PREFIX = run_prefix(TAG)
const FIGDIR = dirname(PREFIX)
const ARROW_DEG = parse(Float64, get(ENV, "MAB_ARROW_DEG", "0.4"))

# ---------------- model ----------------
surface = PREFIX * "_surface_daily.jld2"
grid = global_grid(surface); ug = grid.underlying_grid
Nx, Ny, _ = size(ug)
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center()))
B = grid.immersed_boundary.bottom_height
bh = [B[i, j, 1] for i in 1:Nx, j in 1:Ny]
wet = bh .< 0

Tm, times = surface_series(surface, "T")
um, _ = surface_series(surface, "u")
vm, _ = surface_series(surface, "v")
ηm, ηtimes = surface_series(PREFIX * "_eta_daily.jld2", "η")
tdays = times ./ 86400
days = haskey(ENV, "MAB_DAYS") ? parse.(Float64, split(ENV["MAB_DAYS"], ",")) : [tdays[end]]

# u and v at the cell centres
centre_u(u) = (u[1:end-1, :] .+ u[2:end, :]) ./ 2
centre_v(v) = (v[:, 1:end-1] .+ v[:, 2:end]) ./ 2

# ---------------- GLORYS ----------------
function read_glorys_surface(var, date)
    ds = NCD.Dataset(gfile(var, date))
    lon = Float64.(ds["longitude"][:]); lat = Float64.(ds["latitude"][:])
    raw = ds[var]
    A = ndims(raw) == 4 ? raw[:, :, 1, 1] : raw[:, :, 1]
    close(ds)
    return lon, lat, [ismissing(x) ? NaN : Float64(x) for x in A]
end

# The mean of the two daily means around the de-tided frame's midnight, on the model's cell centres
function glorys_surface_at(var, d)
    day0 = start_date + Day(floor(Int, d)) - Day(1)
    lon, lat, A0 = read_glorys_surface(var, day0)
    _, _, A1 = read_glorys_surface(var, day0 + Day(1))
    A = (A0 .+ A1) ./ 2
    return [wet[i, j] ? bilinear(A, lon, lat, λ[i], φ[j]) : NaN for i in 1:Nx, j in 1:Ny]
end

# ---------------- comparison ----------------
region(i, j) = -bh[i, j] < 200 ? 1 : -bh[i, j] < 1000 ? 2 : 3
region_names = ("shelf (<200 m)", "slope (200-1000 m)", "deep (>1000 m)")
anomaly(X) = X .- mean(filter(isfinite, X))

function statistics(X, Y, r)
    pairs = [(X[i, j], Y[i, j]) for i in 1:Nx, j in 1:Ny if wet[i, j] && region(i, j) == r && isfinite(X[i, j]) && isfinite(Y[i, j])]
    isempty(pairs) && return (NaN, NaN, NaN, 0)
    x = first.(pairs); y = last.(pairs)
    return (mean(x .- y), sqrt(mean((x .- y) .^ 2)), cor(x, y), length(pairs))
end

for dwant in days
    n = argmin(abs.(tdays .- dwant)); d = tdays[n]
    nη = argmin(abs.(ηtimes ./ 86400 .- d))
    date = start_date + Second(round(Int, d * 86400))
    @printf("\n==== %s, day %.1f (%s) ====\n", basename(PREFIX), d, Dates.format(date, "yyyy-mm-dd HH:MM"))

    mask(X) = ifelse.(wet, X, NaN)
    SST = mask(Tm[:, :, n]); SSH = anomaly(mask(ηm[:, :, nη]))
    U = mask(centre_u(um[:, :, n])); V = mask(centre_v(vm[:, :, n]))
    gSST = glorys_surface_at("thetao", d); gSSH = anomaly(glorys_surface_at("zos", d))
    gU = glorys_surface_at("uo", d); gV = glorys_surface_at("vo", d)
    speed = sqrt.(U .^ 2 .+ V .^ 2); gspeed = sqrt.(gU .^ 2 .+ gV .^ 2)

    println("                        bias      rms   correlation   cells")
    for (name, X, Y, unit) in (("SST", SST, gSST, "°C"), ("SSH anomaly", SSH, gSSH, "m"), ("surface speed", speed, gspeed, "m/s"),
                               ("surface u", U, gU, "m/s"), ("surface v", V, gV, "m/s"))
        println(name, " (", unit, ")")
        for r in 1:3
            b, e, c, m = statistics(X, Y, r)
            @printf("  %-20s %+8.3f %8.3f %8.2f %9d\n", region_names[r], b, e, c, m)
        end
    end

    # arrows every ARROW_DEG degrees
    stride = max(1, round(Int, ARROW_DEG / (λ[2] - λ[1])))
    ia = 1:stride:Nx; ja = 1:stride:Ny
    pts = [(i, j) for i in ia, j in ja if wet[i, j]]
    ax_λ = [λ[i] for (i, j) in pts]; ax_φ = [φ[j] for (i, j) in pts]
    scale = ARROW_DEG                                      # 1 m/s spans ARROW_DEG degrees

    fig = Figure(size = (1500, 1300), fontsize = 16)
    Label(fig[0, 1:3], @sprintf("%s day %.1f (%s): model, GLORYS and model − GLORYS", basename(PREFIX), d,
                                Dates.format(date, "yyyy-mm-dd")), fontsize = 20)
    rows = (("SST (°C)", SST, gSST, :thermal, extrema(filter(isfinite, gSST)), 3.0),
            ("SSH anomaly (m)", SSH, gSSH, :balance, (-0.8, 0.8), 0.3),
            ("surface speed (m/s)", speed, gspeed, :speed, (0.0, 1.5), 0.5))
    for (row, (title, X, Y, cmap, crange, drange)) in enumerate(rows)
        for (col, (Z, who)) in enumerate(((X, "model"), (Y, "GLORYS"), (X .- Y, "model − GLORYS")))
            ax = Axis(fig[2row - 1, col], title = "$title: $who", aspect = DataAspect())
            hm = heatmap!(ax, λ, φ, Z; colormap = col == 3 ? :balance : cmap, colorrange = col == 3 ? (-drange, drange) : crange,
                          nan_color = :gray85)
            if row == 3 && col < 3
                Uc, Vc = col == 1 ? (U, V) : (gU, gV)
                arrows2d!(ax, ax_λ, ax_φ, [scale * Uc[i, j] for (i, j) in pts], [scale * Vc[i, j] for (i, j) in pts];
                          color = :black, shaftwidth = 1, tipwidth = 5, tiplength = 5)
            end
            (col == 1 || col == 3) && Colorbar(fig[2row, col == 1 ? (1:2) : 3], hm; vertical = false, flipaxis = false)
        end
    end
    out = joinpath(FIGDIR, @sprintf("%s_surface_vs_glorys_day%02d.png", basename(PREFIX), round(Int, d)))
    save(out, fig; px_per_unit = 1.2); println("saved ", out)
end
