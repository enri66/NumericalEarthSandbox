# De-tided (LowPassFilter daily) SSH animation, contoured, from a run of
# 04_mab_glorys_tides_reservoirs.jl.  Usage:
#   MAB_TAG=mab_glorys_tides_14d julia --project=. animate_ssh_detided.jl
using Oceananigans, CairoMakie, Printf, Statistics, Dates

const DATA = joinpath(ENV["HOME"], "Data", "mab_glorys_obc")
const TAG  = get(ENV, "MAB_TAG", "mab_glorys_tides_14d")
const start_date = DateTime(2019, 4, 1)

ηts = FieldTimeSeries(joinpath(DATA, TAG * "_eta_daily.jld2"), "η")
uts = FieldTimeSeries(joinpath(DATA, TAG * "_surface_daily.jld2"), "u")
t    = ηts.times
grid = ηts.grid
ug   = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
λ = collect(λnodes(ug, Center()))
φ = collect(φnodes(ug, Center()))

# Land mask from the (also de-tided) surface u: Oceananigans zeros velocity in immersed
# cells rather than NaN-ing it, so a cell dry at EVERY frame is land — checked across all
# frames, not just the first, so a momentarily-calm wet cell isn't mistaken for land.
u_allframes = [Array(interior(uts[n]))[:, :, 1] for n in eachindex(t)]
Nu1, Nu2 = size(u_allframes[1])
wet = [any(f -> f[min(i, Nu1), min(j, Nu2)] != 0, u_allframes) for i in eachindex(λ), j in eachindex(φ)]

ηall = [begin
    A = Array(interior(ηts[n]))[:, :, 1]
    map((v, w) -> w ? v : NaN, A, wet)
end for n in eachindex(t)]
lim  = maximum(x -> isfinite(x) ? abs(x) : 0.0, vcat(vec.(ηall)...))
@printf("%d de-tided daily frames, |η| ≤ %.2f m\n", length(t), lim)

levels = range(-lim, lim, length = 41)
lines_levels = range(-lim, lim, length = 11)

n = Observable(1)
ttl = @lift @sprintf("MAB, GLORYS open boundaries + tides — DE-TIDED (LowPassFilter, 1-day) SSH   %s",
                     Dates.format(start_date + Second(round(Int, t[$n])), "yyyy-mm-dd"))
fig = Figure(size = (1400, 1000), fontsize = 20)
ax = Axis(fig[1, 1], title = ttl, xlabel = "longitude (°E)", ylabel = "latitude (°N)",
          aspect = DataAspect())
ηn = @lift ηall[$n]

# land as a flat fill, drawn first, then filled contours + contour lines on top
land = [w ? NaN : 1.0 for w in wet]
heatmap!(ax, λ, φ, land, colormap = [RGBAf(0.72, 0.70, 0.66, 1)], colorrange = (1, 1))
hm = contourf!(ax, λ, φ, ηn, levels = levels, colormap = :balance, extendlow = :auto, extendhigh = :auto)
contour!(ax, λ, φ, ηn, levels = lines_levels, color = (:black, 0.35), linewidth = 0.75)
Colorbar(fig[1, 2], hm, label = "de-tided η (m)")
vlines!(ax, [λ[1], λ[end]], color = (:limegreen, 0.8), linewidth = 2)
hlines!(ax, [φ[1], φ[end]], color = (:limegreen, 0.8), linewidth = 2)

out = joinpath(DATA, TAG * "_ssh_detided.mp4")
CairoMakie.record(fig, out, eachindex(t); framerate = 2, px_per_unit = 2) do i
    n[] = i
end
println("saved ", out)
