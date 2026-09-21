# SST/SSH/surface-speed snapshot panels, from a run of 04_mab_glorys_tides_reservoirs.jl,
# for spot-checking the de-tided daily record rather than only watching the full animation.
# Usage:
#   MAB_TAG=mab_glorys_tides_14d julia --project=. snapshot_fields.jl
#   MAB_SNAPSHOT_DAYS=34,60,87 MAB_TAG=mab_3months julia --project=. snapshot_fields.jl
using Oceananigans, CairoMakie, Printf, Statistics

const DATA = joinpath(ENV["HOME"], "Data", "mab_glorys_obc")
const TAG  = get(ENV, "MAB_TAG", "mab_glorys_tides_14d")

ηts = FieldTimeSeries(joinpath(DATA, TAG * "_eta_daily.jld2"), "η")
Tts = FieldTimeSeries(joinpath(DATA, TAG * "_surface_daily.jld2"), "T")
uts = FieldTimeSeries(joinpath(DATA, TAG * "_surface_daily.jld2"), "u")
vts = FieldTimeSeries(joinpath(DATA, TAG * "_surface_daily.jld2"), "v")
t = ηts.times
grid = ηts.grid
ug = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
λ = collect(λnodes(ug, Center()))
φ = collect(φnodes(ug, Center()))
Nλ, Nφ = length(λ), length(φ)

# u is Face,Center (Nλ+1, Nφ); v is Center,Face (Nλ, Nφ+1) — average onto Center,Center
# (see animate_velocity_detided.jl).
to_center_u(A) = 0.5 .* (A[1:Nλ, :] .+ A[2:Nλ+1, :])
to_center_v(A) = 0.5 .* (A[:, 1:Nφ] .+ A[:, 2:Nφ+1])

# Oceananigans zeros immersed cells rather than NaN-ing them; land is dry (zero T) at
# every frame, checked across all frames so a momentarily-calm wet cell isn't mistaken
# for land.
T_allframes = [Array(interior(Tts[n]))[:, :, 1] for n in eachindex(t)]
wet = [any(f -> f[i, j] != 0, T_allframes) for i in 1:Nλ, j in 1:Nφ]
mask(A) = map((v, w) -> w ? v : NaN, A, wet)

default_days = round.(Int, [t[1], t[end ÷ 2 + 1], t[end]] ./ 86400)
days = [parse(Int, d) for d in split(get(ENV, "MAB_SNAPSHOT_DAYS", join(default_days, ",")), ",")]

for day in days
    n = argmin(abs.(t ./ 86400 .- day))
    actual_day = round(t[n] / 86400, digits = 1)

    T = mask(T_allframes[n])
    η = mask(Array(interior(ηts[n]))[:, :, 1])
    u = to_center_u(Array(interior(uts[n]))[:, :, 1])
    v = to_center_v(Array(interior(vts[n]))[:, :, 1])
    speed = mask(sqrt.(u .^ 2 .+ v .^ 2))

    Tf, ηf, sf = filter(isfinite, T), filter(isfinite, η), filter(isfinite, speed)
    @printf("day %.1f:  T ∈ [%.1f, %.1f] °C  η ∈ [%.2f, %.2f] m  |speed| ≤ %.2f m/s\n",
            actual_day, minimum(Tf), maximum(Tf), minimum(ηf), maximum(ηf), maximum(sf))

    fig = Figure(size = (1500, 400), fontsize = 16)
    ax1 = Axis(fig[1, 1], title = "SST day $actual_day", aspect = DataAspect())
    hm1 = heatmap!(ax1, λ, φ, T, colormap = :thermal)
    Colorbar(fig[1, 2], hm1, label = "°C")
    ax2 = Axis(fig[1, 3], title = "SSH day $actual_day", aspect = DataAspect())
    hm2 = heatmap!(ax2, λ, φ, η, colormap = :balance, colorrange = (-1, 1))
    Colorbar(fig[1, 4], hm2, label = "m")
    ax3 = Axis(fig[1, 5], title = "speed day $actual_day", aspect = DataAspect())
    hm3 = heatmap!(ax3, λ, φ, speed, colormap = :viridis, colorrange = (0, 1.5))
    Colorbar(fig[1, 6], hm3, label = "m/s")

    out = joinpath(DATA, TAG * "_snapshot_day$(Int(round(actual_day))).png")
    save(out, fig)
    println("  saved ", out)
end
