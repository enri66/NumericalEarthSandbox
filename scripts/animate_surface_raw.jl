# Animation of a run's raw (tidal) surface output `<tag>.jld2` (3-hourly in script 05): SST and surface speed with
# velocity vectors, side by side, one frame per output time. Works with script 05's per-rank output. Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/res_test/hr100_fx3 MAB_START_DATE=2019-08-29 \
#       julia --project=/t0/workdir/enrique/mpi05_ib scripts/animate_surface_raw.jl
# MAB_ARROW_DEG:    spacing of the velocity vectors in degrees (default 0.4)
# MAB_FRAME_STRIDE: use every n-th output time (default 1)
# MAB_FRAMERATE:    frames per second of video (default 12)
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const PREFIX = run_prefix(get(ENV, "MAB_TAG", "mab_H90"))
const ARROW_DEG = parse(Float64, get(ENV, "MAB_ARROW_DEG", "0.4"))
const STRIDE = parse(Int, get(ENV, "MAB_FRAME_STRIDE", "1"))
const FPS = parse(Int, get(ENV, "MAB_FRAMERATE", "12"))

path = PREFIX * ".jld2"
grid = global_grid(path); ug = grid.underlying_grid
Nx, Ny = ug.Nx, ug.Ny
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center()))
B = grid.immersed_boundary.bottom_height
wet = [B[i, j, 1] < 0 for i in 1:Nx, j in 1:Ny]
mask(A) = ifelse.(wet, A, NaN)

T, times = surface_series(path, "T")
u, _ = surface_series(path, "u"); v, _ = surface_series(path, "v")
# a restart from a checkpoint rewrites the frames after it: keep each output time once, in order
order = sortperm(times); keep = [order[k] for k in eachindex(order) if k == length(order) || times[order[k+1]] > times[order[k]]]
frames = keep[1:STRIDE:end]
@printf("%s: %d frames, days %.2f to %.2f\n", basename(PREFIX), length(frames), times[frames[1]] / 86400, times[frames[end]] / 86400)

uc(n) = mask((u[1:Nx, :, n] .+ u[2:Nx+1, :, n]) ./ 2)
vc(n) = mask((v[:, 1:Ny, n] .+ v[:, 2:Ny+1, n]) ./ 2)

stride = max(1, round(Int, ARROW_DEG / (λ[2] - λ[1])))
pts = [(i, j) for i in 1:stride:Nx, j in 1:stride:Ny if wet[i, j]]
λa = [λ[i] for (i, j) in pts]; φa = [φ[j] for (i, j) in pts]

k = Observable(1)
n = @lift frames[$k]
sst = @lift mask(T[:, :, $n])
U = @lift uc($n); V = @lift vc($n)
speed = @lift sqrt.($U .^ 2 .+ $V .^ 2)
ua = @lift [ARROW_DEG * $U[i, j] for (i, j) in pts]; va = @lift [ARROW_DEG * $V[i, j] for (i, j) in pts]
title = @lift @sprintf("%s  %s (day %.2f)", basename(PREFIX),
                       Dates.format(start_date + Second(round(Int, times[$n])), "yyyy-mm-dd HH:MM"), times[$n] / 86400)

fig = Figure(size = (1800, 760), fontsize = 18)
Label(fig[0, 1:2], title, fontsize = 22)
ax1 = Axis(fig[1, 1], title = "SST (°C)", aspect = DataAspect())
hm1 = heatmap!(ax1, λ, φ, sst; colormap = :thermal, colorrange = (14, 30), nan_color = :gray80)
Colorbar(fig[2, 1], hm1; vertical = false, flipaxis = false)
ax2 = Axis(fig[1, 2], title = "surface speed (m/s), arrows: 1 m/s = $(ARROW_DEG)°", aspect = DataAspect())
hm2 = heatmap!(ax2, λ, φ, speed; colormap = :speed, colorrange = (0, 2), nan_color = :gray80)
arrows2d!(ax2, λa, φa, ua, va; color = :black, shaftwidth = 1, tipwidth = 4, tiplength = 4)
Colorbar(fig[2, 2], hm2; vertical = false, flipaxis = false)

out = PREFIX * "_surface_raw.mp4"
record(fig, out, eachindex(frames); framerate = FPS) do f
    k[] = f
end
println("saved ", out)
