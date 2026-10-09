# A quick look at a run in progress: SST, SSS and surface speed from the latest frame of its raw (tidal) surface output
# `<tag>.jld2` (3-hourly in script 05), which exists long before the de-tided daily frames. Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/res_test/hr100_fx3 MAB_START_DATE=2019-08-29 \
#       julia --project=/t0/workdir/enrique/mpi05_ib scripts/surface_snapshot.jl
# MAB_FRAME: the frame to plot (default: the last)
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const PREFIX = run_prefix(get(ENV, "MAB_TAG", "mab_H90"))
path = PREFIX * ".jld2"
grid = global_grid(path); ug = grid.underlying_grid
Nx, Ny, _ = size(ug)
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center()))
B = grid.immersed_boundary.bottom_height
wet = [B[i, j, 1] < 0 for i in 1:Nx, j in 1:Ny]

T, times = surface_series(path, "T")
S, _ = surface_series(path, "S")
u, _ = surface_series(path, "u")
v, _ = surface_series(path, "v")
n = haskey(ENV, "MAB_FRAME") ? parse(Int, ENV["MAB_FRAME"]) : length(times)
mask(X) = ifelse.(wet, X, NaN)
speed = mask(sqrt.(((u[1:end-1, :, n] .+ u[2:end, :, n]) ./ 2) .^ 2 .+ ((v[:, 1:end-1, n] .+ v[:, 2:end, n]) ./ 2) .^ 2))
date = start_date + Second(round(Int, times[n]))

fig = Figure(size = (1800, 620), fontsize = 16)
Label(fig[0, 1:3], @sprintf("%s, %s (day %.2f), raw surface output", basename(PREFIX), Dates.format(date, "yyyy-mm-dd HH:MM"),
                            times[n] / 86400), fontsize = 20)
for (c, (X, title, cmap, crange)) in enumerate(((mask(T[:, :, n]), "SST (°C)", :thermal, (14, 30)),
                                                (mask(S[:, :, n]), "SSS (psu)", :haline, (30, 36.5)),
                                                (speed, "surface speed (m/s)", :speed, (0, 2))))
    ax = Axis(fig[1, c], title = title, aspect = DataAspect())
    hm = heatmap!(ax, λ, φ, X; colormap = cmap, colorrange = crange, nan_color = :gray85)
    Colorbar(fig[2, c], hm; vertical = false, flipaxis = false)
end
out = PREFIX * @sprintf("_snapshot_day%05.2f.png", times[n] / 86400)
save(out, fig; px_per_unit = 1); println("saved ", out)
