# Close-ups of a run's raw surface output around chosen points, cell by cell: ζ/f, T, u and v, to tell a real front from
# grid-scale noise. Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/res_test/hr100_fx3 MAB_START_DATE=2019-08-29 MAB_POINTS="106:55,535:322" \
#       julia --project=/t0/workdir/enrique/mpi05_ib scripts/streak_zoom.jl
# MAB_POINTS: "i:j,..." centres of the close-ups; MAB_HALF: half-width in cells (default 20); MAB_FRAME: frame (default middle)
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const PREFIX = run_prefix(get(ENV, "MAB_TAG", "mab_H90"))
const HALF = parse(Int, get(ENV, "MAB_HALF", "20"))
points = [Tuple(parse.(Int, split(p, ":"))) for p in split(ENV["MAB_POINTS"], ",")]

path = PREFIX * ".jld2"
grid = global_grid(path); ug = grid.underlying_grid
Nx, Ny = ug.Nx, ug.Ny
B = grid.immersed_boundary.bottom_height
wet = [B[i, j, 1] < 0 for i in 1:Nx, j in 1:Ny]
φc = collect(φnodes(ug, Center())); φf = collect(φnodes(ug, Face()))
Δλ = deg2rad(λnodes(ug, Center())[2] - λnodes(ug, Center())[1]); Δφ = deg2rad(φc[2] - φc[1])

T, times = surface_series(path, "T"); u, _ = surface_series(path, "u"); v, _ = surface_series(path, "v")
n = haskey(ENV, "MAB_FRAME") ? parse(Int, ENV["MAB_FRAME"]) : length(times) ÷ 2
U = u[:, :, n]; V = v[:, :, n]
Ro = fill(NaN, Nx, Ny)
for i in 2:Nx, j in 2:Ny
    (wet[i-1, j-1] & wet[i, j-1] & wet[i-1, j] & wet[i, j]) || continue
    ζ = ((V[i, j] - V[i-1, j]) / Δλ - (U[i, j] * cosd(φc[j]) - U[i, j-1] * cosd(φc[j-1])) / Δφ) / (6.371e6 * cosd(φf[j]))
    Ro[i, j] = ζ / (2 * 7.292115e-5 * sind(φf[j]))
end

fig = Figure(size = (1700, 430 * length(points)), fontsize = 14)
for (r, (ic, jc)) in enumerate(points)
    is = max(2, ic - HALF):min(Nx - 1, ic + HALF); js = max(2, jc - HALF):min(Ny - 1, jc + HALF)
    for (c, (name, A, cmap)) in enumerate((("ζ/f", Ro, :balance), ("T (°C)", T[:, :, n], :thermal),
                                           ("u (m/s)", U[1:Nx, :], :balance), ("v (m/s)", V[:, 1:Ny], :balance)))
        X = A[is, js]
        lim = name == "T (°C)" ? extrema(filter(isfinite, X)) : (m = maximum(abs, filter(isfinite, X)); (-m, m))
        ax = Axis(fig[r, c], title = "$name around ($ic, $jc)", xlabel = "i", ylabel = "j", aspect = DataAspect())
        hm = heatmap!(ax, collect(is), collect(js), X; colormap = cmap, colorrange = lim, nan_color = :gray80)
        Colorbar(fig[r, c][1, 2], hm)
    end
end
Label(fig[0, :], @sprintf("%s, frame %d (day %.2f), cell-by-cell close-ups", basename(PREFIX), n, times[n] / 86400), fontsize = 18)
out = PREFIX * "_streak_zoom.png"
save(out, fig; px_per_unit = 1); println("saved ", out)
