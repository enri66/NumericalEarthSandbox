# De-tided (FilteredTimeInterval daily) surface velocity animation, from a run of
# 04_mab_glorys_tides_reservoirs.jl: speed as filled contours, direction as
# subsampled arrows. Usage:
#   MAB_TAG=mab_glorys_tides_14d julia --project=. animate_velocity_detided.jl
using Oceananigans, CairoMakie, Printf, Statistics, Dates

const DATA   = joinpath(ENV["HOME"], "Data", "mab_glorys_obc")
const TAG    = get(ENV, "MAB_TAG", "mab_glorys_tides_14d")
const STRIDE = parse(Int, get(ENV, "MAB_QUIVER_STRIDE", "4"))   # every STRIDE-th cell gets an arrow
const start_date = DateTime(2019, 4, 1)

uts = FieldTimeSeries(joinpath(DATA, TAG * "_surface_daily.jld2"), "u")
vts = FieldTimeSeries(joinpath(DATA, TAG * "_surface_daily.jld2"), "v")
t    = uts.times
grid = uts.grid
ug   = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
λ = collect(λnodes(ug, Center()))
φ = collect(φnodes(ug, Center()))
Nλ, Nφ = length(λ), length(φ)

# u is Face,Center (Nλ+1, Nφ); v is Center,Face (Nλ, Nφ+1). Average neighbouring faces onto
# the Center,Center grid so a single arrow represents both components at the same point.
to_center_u(A) = 0.5 .* (A[1:Nλ, :] .+ A[2:Nλ+1, :])
to_center_v(A) = 0.5 .* (A[:, 1:Nφ] .+ A[:, 2:Nφ+1])

uall = [to_center_u(Array(interior(uts[n]))[:, :, 1]) for n in eachindex(t)]
vall = [to_center_v(Array(interior(vts[n]))[:, :, 1]) for n in eachindex(t)]

# Land mask: a cell dry at every frame in BOTH components (see animate_ssh_detided.jl — this
# output convention zeros immersed cells rather than NaN-ing them).
wet = [any(n -> uall[n][i, j] != 0 || vall[n][i, j] != 0, eachindex(t)) for i in 1:Nλ, j in 1:Nφ]

speed_all = [begin
    s = sqrt.(uall[n] .^ 2 .+ vall[n] .^ 2)
    map((v, w) -> w ? v : NaN, s, wet)
end for n in eachindex(t)]
lim = maximum(x -> isfinite(x) ? x : 0.0, vcat(vec.(speed_all)...))
@printf("%d de-tided daily velocity frames, |speed| ≤ %.2f m/s\n", length(t), lim)

iλ = 1:STRIDE:Nλ
iφ = 1:STRIDE:Nφ
λs, φs = λ[iλ], φ[iφ]

n = Observable(1)
ttl = @lift @sprintf("MAB, GLORYS open boundaries + tides — DE-TIDED surface velocity   %s",
                     Dates.format(start_date + Second(round(Int, t[$n])), "yyyy-mm-dd"))
fig = Figure(size = (1400, 1000), fontsize = 20)
ax = Axis(fig[1, 1], title = ttl, xlabel = "longitude (°E)", ylabel = "latitude (°N)",
          aspect = DataAspect())
speed_n = @lift speed_all[$n]
us = @lift [wet[i, j] ? uall[$n][i, j] : 0.0 for i in iλ, j in iφ]
vs = @lift [wet[i, j] ? vall[$n][i, j] : 0.0 for i in iλ, j in iφ]

land = [w ? NaN : 1.0 for w in wet]
heatmap!(ax, λ, φ, land, colormap = [RGBAf(0.72, 0.70, 0.66, 1)], colorrange = (1, 1))
hm = contourf!(ax, λ, φ, speed_n, levels = range(0, lim, length = 31), colormap = :speed, extendhigh = :auto)
arrows2d!(ax, λs, φs, us, vs; lengthscale = 0.6, color = (:black, 0.75), minshaftlength = 0)
Colorbar(fig[1, 2], hm, label = "de-tided |u| (m/s)")
vlines!(ax, [λ[1], λ[end]], color = (:limegreen, 0.8), linewidth = 2)
hlines!(ax, [φ[1], φ[end]], color = (:limegreen, 0.8), linewidth = 2)

out = joinpath(DATA, TAG * "_velocity_detided.mp4")
CairoMakie.record(fig, out, eachindex(t); framerate = 2, px_per_unit = 2) do i
    n[] = i
end
println("saved ", out)
