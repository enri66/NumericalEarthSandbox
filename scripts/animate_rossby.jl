# Animation of the surface Rossby number ζ/f from a run's raw (tidal) surface output `<tag>.jld2` (3-hourly in script
# 05). ζ = ∂v/∂x − ∂u/∂y at the cell corners of the latitude-longitude grid,
#   ζ = [∂v/∂λ − ∂(u cos φ)/∂φ] / (R cos φ),
# and f = 2Ω sin φ there; corners next to land are left out. Works with script 05's per-rank output. Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/res_test/hr100_fx3 MAB_START_DATE=2019-08-29 \
#       julia --project=/t0/workdir/enrique/mpi05_ib scripts/animate_rossby.jl
# MAB_RO_MAX:       colour range ±MAB_RO_MAX (default 1)
# MAB_FRAME_STRIDE: use every n-th output time (default 1)
# MAB_FRAMERATE:    frames per second of video (default 12)
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const PREFIX = run_prefix(get(ENV, "MAB_TAG", "mab_H90"))
const RO_MAX = parse(Float64, get(ENV, "MAB_RO_MAX", "1"))
const STRIDE = parse(Int, get(ENV, "MAB_FRAME_STRIDE", "1"))
const FPS = parse(Int, get(ENV, "MAB_FRAMERATE", "12"))
const R = 6.371e6
const Ω = 7.292115e-5

path = PREFIX * ".jld2"
grid = global_grid(path); ug = grid.underlying_grid
Nx, Ny = ug.Nx, ug.Ny
λc = collect(λnodes(ug, Center())); φc = collect(φnodes(ug, Center()))
λf = collect(λnodes(ug, Face()));   φf = collect(φnodes(ug, Face()))
Δλ = deg2rad(λc[2] - λc[1]); Δφ = deg2rad(φc[2] - φc[1])
B = grid.immersed_boundary.bottom_height
wet = [B[i, j, 1] < 0 for i in 1:Nx, j in 1:Ny]
# interior corners (i, j) between cells i-1, i and j-1, j; kept only where all four cells are wet
corner_wet = [wet[i-1, j-1] & wet[i, j-1] & wet[i-1, j] & wet[i, j] for i in 2:Nx, j in 2:Ny]
f = [2Ω * sind(φf[j]) for i in 2:Nx, j in 2:Ny]
cosφf = [cosd(φf[j]) for i in 2:Nx, j in 2:Ny]

u, times = surface_series(path, "u")
v, _ = surface_series(path, "v")
order = sortperm(times); keep = [order[k] for k in eachindex(order) if k == length(order) || times[order[k+1]] > times[order[k]]]
frames = keep[1:STRIDE:end]
@printf("%s: %d frames, days %.2f to %.2f\n", basename(PREFIX), length(frames), times[frames[1]] / 86400, times[frames[end]] / 86400)

# u is (Nx+1) × Ny on x-faces, v is Nx × (Ny+1) on y-faces
function rossby(n)
    U = view(u, :, :, n); V = view(v, :, :, n)
    ζ = [((V[i, j] - V[i-1, j]) / Δλ - (U[i, j] * cosd(φc[j]) - U[i, j-1] * cosd(φc[j-1])) / Δφ) for i in 2:Nx, j in 2:Ny]
    Ro = ζ ./ (R .* cosφf) ./ f
    return ifelse.(corner_wet, Ro, NaN)
end

k = Observable(1)
Ro = @lift rossby(frames[$k])
title = @lift @sprintf("%s  surface Rossby number ζ/f  %s (day %.2f)", basename(PREFIX),
                       Dates.format(start_date + Second(round(Int, times[frames[$k]])), "yyyy-mm-dd HH:MM"), times[frames[$k]] / 86400)

# the figure takes the map's shape (degrees of longitude shortened by cos φ), with the colour bar at the side
aspect = (λf[end] - λf[1]) * cosd((φf[1] + φf[end]) / 2) / (φf[end] - φf[1])
fig = Figure(size = (round(Int, 1000 * aspect) + 220, 1100), fontsize = 18)
Label(fig[0, 1:2], title, fontsize = 22)
ax = Axis(fig[1, 1], aspect = DataAspect())
hm = heatmap!(ax, λf[2:Nx], φf[2:Ny], Ro; colormap = :balance, colorrange = (-RO_MAX, RO_MAX), nan_color = :gray80)
Colorbar(fig[1, 2], hm; label = "ζ / f")
resize_to_layout!(fig)

out = PREFIX * "_rossby_raw.mp4"
record(fig, out, eachindex(frames); framerate = FPS) do i
    k[] = i
end
println("saved ", out)
