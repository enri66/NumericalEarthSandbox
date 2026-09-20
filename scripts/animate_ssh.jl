# Animate SSH from a run of 02_mab_glorys_obc.jl.  Usage:
#   MAB_TAG=mab_obc_smoke julia --project=. scripts/animate_ssh.jl
using Oceananigans, CairoMakie, Printf, Statistics, Dates
const TAG = get(ENV, "MAB_TAG", "mab_obc_smoke")
const start_date = DateTime(2019, 4, 1)

ηts = FieldTimeSeries(joinpath(@__DIR__, "..", TAG * "_eta.jld2"), "η")
uts = FieldTimeSeries(joinpath(@__DIR__, "..", TAG * ".jld2"), "u")
t   = ηts.times
grid = ηts.grid
ug = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
λ = collect(λnodes(ug, Center()))
φ = collect(φnodes(ug, Center()))
@assert λ isa AbstractVector && φ isa AbstractVector "grid nodes are not 1-D"


# land mask from the first velocity frame (immersed cells never get written)
u1 = Array(interior(uts[1]))[:, :, 1]
mask = [abs(u1[min(i,size(u1,1)), min(j,size(u1,2))]) > 0 || u1[min(i,size(u1,1)), min(j,size(u1,2))] == 0 ? 1.0 : NaN
        for i in eachindex(λ), j in eachindex(φ)]
wet = [isfinite(u1[min(i,size(u1,1)), min(j,size(u1,2))]) for i in eachindex(λ), j in eachindex(φ)]

ηall = [Array(interior(ηts[n]))[:, :, 1] for n in eachindex(t)]
lim  = maximum(x -> isfinite(x) ? abs(x) : 0.0, vcat(vec.(ηall)...))
@printf("%d frames, |η| ≤ %.2f m\n", length(t), lim)

n = Observable(1)
ttl = @lift @sprintf("MAB, GLORYS open boundaries — SSH   %s",
                     Dates.format(start_date + Second(round(Int, t[$n])), "yyyy-mm-dd HH:MM"))
fig = Figure(size = (980, 720))
ax = Axis(fig[1, 1], title = ttl, xlabel = "longitude (°E)", ylabel = "latitude (°N)",
          aspect = DataAspect())
ηn = @lift map(x -> isfinite(x) ? x : NaN, ηall[$n])
hm = heatmap!(ax, λ, φ, ηn, colormap = :balance, colorrange = (-lim, lim),
              nan_color = RGBAf(0.72, 0.70, 0.66, 1))
Colorbar(fig[1, 2], hm, label = "η (m)")
# mark the open boundaries
vlines!(ax, [λ[1], λ[end]], color = (:limegreen, 0.8), linewidth = 2)
hlines!(ax, [φ[1], φ[end]], color = (:limegreen, 0.8), linewidth = 2)

out = joinpath(@__DIR__, "..", TAG * "_ssh.mp4")
CairoMakie.record(fig, out, eachindex(t); framerate = 6) do i
    n[] = i
end
println("saved ", out)
