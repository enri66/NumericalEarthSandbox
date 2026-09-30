# Animation of a run's de-tided daily surface fields: SST, SSH (domain mean removed) and surface speed with velocity
# vectors, one frame per day. Reads `<tag>_surface_daily.jld2` and `<tag>_eta_daily.jld2`, including script 05's
# per-rank output. Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/sponge_year_mpi2/spy julia --project=. scripts/animate_surface_year.jl
# MAB_ARROW_STRIDE: plot every n-th velocity vector in each direction (default 5)
# MAB_FRAMERATE:    frames (days) per second of video (default 12)
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const TAG    = get(ENV, "MAB_TAG", "mab_H90")
const PREFIX = run_prefix(TAG)
const STRIDE = parse(Int, get(ENV, "MAB_ARROW_STRIDE", "5"))
const FPS    = parse(Int, get(ENV, "MAB_FRAMERATE", "12"))

surf = PREFIX * "_surface_daily.jld2"
T, times = surface_series(surf, "T")
u, _ = surface_series(surf, "u"); v, _ = surface_series(surf, "v")
η, ηtimes = surface_series(PREFIX * "_eta_daily.jld2", "η")
@assert ηtimes ≈ times "surface and SSH daily frames differ"
grid = open_series(PREFIX * "_volume_daily.jld2", "T"; backend = OnDisk()).grid
ug = grid.underlying_grid
Nx, Ny = ug.Nx, ug.Ny
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center()))
B = grid.immersed_boundary.bottom_height
wet = [B[i, j, 1] < 0 for i in 1:Nx, j in 1:Ny]
mask(A) = [wet[i, j] ? A[i, j] : NaN for i in 1:Nx, j in 1:Ny]
Nt = length(times)
@printf("%s: %d daily frames, days %.1f to %.1f\n", basename(PREFIX), Nt, times[1] / 86400, times[end] / 86400)

# u and v on their faces, averaged to cell centres
function velocities(n)
    uc = [(u[i, j, n] + u[i+1, j, n]) / 2 for i in 1:Nx, j in 1:Ny]
    vc = [(v[i, j, n] + v[i, j+1, n]) / 2 for i in 1:Nx, j in 1:Ny]
    return mask(uc), mask(vc)
end
function ssh(n)
    h = mask(η[1:Nx, 1:Ny, n])
    return h .- mean(filter(isfinite, h))
end

# arrows on a coarser lattice; the arrow for 1 m/s spans `scale` degrees
ia = 1:STRIDE:Nx; ja = 1:STRIDE:Ny
λa = [λ[i] for i in ia, j in ja]; φa = [φ[j] for i in ia, j in ja]
const scale = 0.5

n = Observable(1)
SST = @lift mask(T[1:Nx, 1:Ny, $n])
SSH = @lift ssh($n)
UV  = @lift velocities($n)
speed = @lift sqrt.($UV[1] .^ 2 .+ $UV[2] .^ 2)
ua = @lift [isfinite($UV[1][i, j]) ? scale * $UV[1][i, j] : 0.0 for i in ia, j in ja]
va = @lift [isfinite($UV[2][i, j]) ? scale * $UV[2][i, j] : 0.0 for i in ia, j in ja]
date = @lift Dates.format(start_date + Second(round(Int, times[$n])), "yyyy-mm-dd")

fig = Figure(size = (2100, 760), fontsize = 18)
Label(fig[0, 1:3], @lift("$(basename(PREFIX)): de-tided daily surface fields, $($date)"), fontsize = 24)
panels = ((SST, "SST (°C)", :thermal, (2, 30)), (SSH, "SSH, domain mean removed (m)", :balance, (-1, 1)),
          (speed, "surface speed (m/s) and velocity", :deep, (0, 2)))
for (c, (obs, title, cmap, crange)) in enumerate(panels)
    ax = Axis(fig[1, c], title = title, aspect = DataAspect(), xlabel = "longitude", ylabel = c == 1 ? "latitude" : "")
    hm = heatmap!(ax, λ, φ, obs; colormap = cmap, colorrange = crange, nan_color = :gray80)
    c == 3 && arrows2d!(ax, vec(λa), vec(φa), @lift(vec($ua)), @lift(vec($va)); color = :black,
                        shaftwidth = 1, tipwidth = 5, tiplength = 5)
    Colorbar(fig[2, c], hm; vertical = false, flipaxis = false)
end

out = PREFIX * "_surface_year.mp4"
record(fig, out, 1:Nt; framerate = FPS) do k
    n[] = k
end
println("saved ", out)
