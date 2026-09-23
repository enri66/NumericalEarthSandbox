# De-tided (FilteredTimeInterval daily) SST animation, contoured, from a run of
# 04_mab_glorys_tides_reservoirs.jl.  Usage:
#   MAB_TAG=mab_glorys_tides_14d julia --project=. animate_sst_detided.jl
using Oceananigans, CairoMakie, Printf, Statistics, Dates

const DATA = joinpath(ENV["HOME"], "Data", "mab_glorys_obc")
const TAG  = get(ENV, "MAB_TAG", "mab_glorys_tides_14d")
const start_date = DateTime(2019, 4, 1)

Tts = FieldTimeSeries(joinpath(DATA, TAG * "_surface_daily.jld2"), "T")
uts = FieldTimeSeries(joinpath(DATA, TAG * "_surface_daily.jld2"), "u")
t    = Tts.times
grid = Tts.grid
ug   = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
λ = collect(λnodes(ug, Center()))
φ = collect(φnodes(ug, Center()))

# Land mask from u, checked across all frames (see animate_ssh_detided.jl — this output
# convention zeros immersed cells rather than NaN-ing them).
u_allframes = [Array(interior(uts[n]))[:, :, 1] for n in eachindex(t)]
Nu1, Nu2 = size(u_allframes[1])
wet = [any(f -> f[min(i, Nu1), min(j, Nu2)] != 0, u_allframes) for i in eachindex(λ), j in eachindex(φ)]

Tall = [begin
    A = Array(interior(Tts[n]))[:, :, 1]
    map((v, w) -> w ? v : NaN, A, wet)
end for n in eachindex(t)]
lo = minimum(x -> isfinite(x) ? x : Inf,  vcat(vec.(Tall)...))
hi = maximum(x -> isfinite(x) ? x : -Inf, vcat(vec.(Tall)...))
@printf("%d de-tided daily SST frames, T ∈ [%.1f, %.1f] °C\n", length(t), lo, hi)

levels = range(lo, hi, length = 41)
lines_levels = range(lo, hi, length = 13)

n = Observable(1)
ttl = @lift @sprintf("MAB, GLORYS open boundaries + tides — DE-TIDED (low-pass filtered, 1-day) SST   %s",
                     Dates.format(start_date + Second(round(Int, t[$n])), "yyyy-mm-dd"))
fig = Figure(size = (1400, 1000), fontsize = 20)
ax = Axis(fig[1, 1], title = ttl, xlabel = "longitude (°E)", ylabel = "latitude (°N)",
          aspect = DataAspect())
Tn = @lift Tall[$n]

land = [w ? NaN : 1.0 for w in wet]
heatmap!(ax, λ, φ, land, colormap = [RGBAf(0.72, 0.70, 0.66, 1)], colorrange = (1, 1))
hm = contourf!(ax, λ, φ, Tn, levels = levels, colormap = :thermal, extendlow = :auto, extendhigh = :auto)
contour!(ax, λ, φ, Tn, levels = lines_levels, color = (:black, 0.35), linewidth = 0.75)
Colorbar(fig[1, 2], hm, label = "de-tided SST (°C)")
vlines!(ax, [λ[1], λ[end]], color = (:limegreen, 0.8), linewidth = 2)
hlines!(ax, [φ[1], φ[end]], color = (:limegreen, 0.8), linewidth = 2)

out = joinpath(DATA, TAG * "_sst_detided.mp4")
CairoMakie.record(fig, out, eachindex(t); framerate = 2, px_per_unit = 2) do i
    n[] = i
end
println("saved ", out)
