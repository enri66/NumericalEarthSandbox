# Sea-surface salinity of the river runs against their control: the control's SSS and, for each river run, its SSS minus
# the control's, from the de-tided daily surface output, on one day (default: the last daily frame). Usage:
#   RUNS=/t0/workdir/enrique/runs/res_test/rv1_ctl,/t0/workdir/enrique/runs/res_test/rv1_on,/t0/workdir/enrique/runs/res_test/rv1_nosss \
#   MAB_START_DATE=2019-08-29 julia --project=/t0/workdir/enrique/mpi05_ib scripts/river_sss.jl
# RUNS:     path prefixes, the control first
# MAB_DAYS: the day to plot (default: the last daily frame)
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const RUNS = split(ENV["RUNS"], ",")
names = basename.(RUNS)

grid = global_grid(RUNS[1] * "_surface_daily.jld2"); ug = grid.underlying_grid
Nx, Ny, _ = size(ug)
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center()))
B = grid.immersed_boundary.bottom_height
wet = [B[i, j, 1] < 0 for i in 1:Nx, j in 1:Ny]

series = [surface_series(r * "_surface_daily.jld2", "S") for r in RUNS]
tdays = series[1][2] ./ 86400
d = haskey(ENV, "MAB_DAYS") ? parse(Float64, first(split(ENV["MAB_DAYS"], ","))) : tdays[end]
SSS = map(series) do (S, times)
    n = argmin(abs.(times ./ 86400 .- d))
    ifelse.(wet, S[:, :, n], NaN)
end
date = start_date + Second(round(Int, d * 86400))

for (r, name) in enumerate(names[2:end])
    Δ = SSS[r + 1] .- SSS[1]
    i, j = Tuple(argmin(replace(Δ, NaN => Inf)))
    @printf("%s − %s on %s: freshest %.2f psu at %.2f°E %.2f°N; %d cells fresher by > 0.1 psu, %d by > 1 psu\n", name, names[1],
            Dates.format(date, "yyyy-mm-dd"), Δ[i, j], λ[i], φ[j], count(<(-0.1), filter(isfinite, Δ)), count(<(-1), filter(isfinite, Δ)))
end

fig = Figure(size = (650 * length(RUNS), 650), fontsize = 16)
Label(fig[0, 1:length(RUNS)], @sprintf("sea-surface salinity, de-tided daily mean, %s (day %.0f)", Dates.format(date, "yyyy-mm-dd"), d), fontsize = 20)
ax = Axis(fig[1, 1], title = "$(names[1]) SSS (psu)", aspect = DataAspect())
hm = heatmap!(ax, λ, φ, SSS[1]; colormap = :haline, colorrange = (30, 36.5), nan_color = :gray85)
Colorbar(fig[2, 1], hm; vertical = false, flipaxis = false)
for (r, name) in enumerate(names[2:end])
    ax = Axis(fig[1, r + 1], title = "$(name) − $(names[1]) (psu)", aspect = DataAspect())
    hm = heatmap!(ax, λ, φ, SSS[r + 1] .- SSS[1]; colormap = :balance, colorrange = (-2, 2), nan_color = :gray85)
    Colorbar(fig[2, r + 1], hm; vertical = false, flipaxis = false)
end
out = joinpath(dirname(RUNS[1]), @sprintf("river_sss_day%02d.png", round(Int, d)))
save(out, fig; px_per_unit = 1.2); println("saved ", out)
