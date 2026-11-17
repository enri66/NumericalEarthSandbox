# The Mid-Atlantic Bight cold pool in a long run against GLORYS: the temperature of the bottom cell on the MAB shelf
# (20-200 m deep, west of 69.5°W, south of 41.6°N) from the de-tided daily 3D output, every day, as
#   time series of the cold-pool area (bottom temperature below 10 °C) and of the shelf-mean bottom temperature, and
#   monthly mean maps of the bottom temperature (model, GLORYS and model − GLORYS).
# The model's bottom cell is the deepest one at least PARTIAL_MIN full (as script 05's partial cells); GLORYS's bottom
# temperature is its deepest valid level, interpolated to the model's cells from the two daily means around each midnight.
# Usage (a run in segments):
#   MAB_TAGS=/t0/.../mabyear_s01,... MAB_OUT=/t0/.../analysis/year/mabyear MAB_START_DATE=2019-04-01 \
#   MAB_DATA_DIR=/t0/workdir/enrique/Data/NumericalEarth julia --project=<run env> scripts/cold_pool.jl
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const PREFIXES = run_prefixes()
const OUT = get(ENV, "MAB_OUT", PREFIXES[1])
const PARTIAL_MIN = parse(Float64, get(ENV, "MAB_PARTIAL_MIN", "0.2"))
const COLD = 10.0                                         # °C
const BOX = ((-76.0, -69.5), (35.0, 41.6))                # MAB shelf, without Georges Bank and the Gulf of Maine
const HRANGE = (20.0, 200.0)

windows = segment_windows(PREFIXES, "_volume_daily")
series = [open_series(p * "_volume_daily.jld2", "T"; backend = OnDisk()) for p in PREFIXES]
grid = series[1].grid; ug = grid.underlying_grid
Nx, Ny, Nz = size(ug)
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center())); zf = collect(znodes(ug, Face()))
B = grid.immersed_boundary.bottom_height
bh = [B[i, j, 1] for i in 1:Nx, j in 1:Ny]
shelf = [BOX[1][1] <= λ[i] <= BOX[1][2] && BOX[2][1] <= φ[j] <= BOX[2][2] && HRANGE[1] <= -bh[i, j] <= HRANGE[2]
         for i in 1:Nx, j in 1:Ny]
# the bottom cell: the deepest one whose wet part is at least PARTIAL_MIN of its height
kbot = [let k = findfirst(k -> zf[k+1] - max(zf[k], bh[i, j]) >= PARTIAL_MIN * (zf[k+1] - zf[k]), 1:Nz)
            isnothing(k) ? 0 : k end for i in 1:Nx, j in 1:Ny]
const Rₑ = 6.371e6
area = [Rₑ^2 * cosd(φ[j]) * deg2rad(λ[2] - λ[1]) * deg2rad(φ[2] - φ[1]) / 1e6 for i in 1:Nx, j in 1:Ny]   # km²
@printf("MAB shelf: %d cells, %.0f km²\n", count(shelf), sum(area[shelf]))

# ---------------- GLORYS bottom temperature ----------------
const GFILES = Dict{Date, String}()
for f in readdir(GDIR)
    m = match(r"^thetao_GLORYSDaily_(\d{4}-\d{2}-\d{2})T.*_-76\.0_-64\.0_34\.0_42\.0\.nc$", f)
    isnothing(m) || (GFILES[Date(m[1])] = joinpath(GDIR, f))
end
function glorys_bottom(date)
    ds = NCD.Dataset(GFILES[Date(date)])
    lon = Float64.(ds["longitude"][:]); lat = Float64.(ds["latitude"][:])
    A = ds["thetao"][:, :, :, 1]
    close(ds)
    bottom = fill(NaN, length(lon), length(lat))
    for j in eachindex(lat), i in eachindex(lon)
        k = findlast(k -> !ismissing(A[i, j, k]) && isfinite(A[i, j, k]), axes(A, 3))
        isnothing(k) || (bottom[i, j] = Float64(A[i, j, k]))
    end
    return lon, lat, bottom
end
function glorys_bottom_at(d)
    day0 = start_date + Day(floor(Int, d)) - Day(1)
    lon, lat, A0 = glorys_bottom(day0); _, _, A1 = glorys_bottom(day0 + Day(1))
    A = (A0 .+ A1) ./ 2
    return [shelf[i, j] ? bilinear(A, lon, lat, λ[i], φ[j]) : NaN for i in 1:Nx, j in 1:Ny]
end

# ---------------- daily bottom temperature ----------------
fdays = Float64[]; model_area = Float64[]; glorys_area = Float64[]; model_mean = Float64[]; glorys_mean = Float64[]
months = Date[]; msum = Dict{Tuple{Symbol, Date}, Matrix{Float64}}(); mcount = Dict{Date, Matrix{Int}}()
for (s, TV) in enumerate(series)
    for n in eachindex(TV.times)
        t = TV.times[n]
        (windows[s][1] <= t < windows[s][2]) || continue
        (!isempty(fdays) && t / 86400 <= fdays[end]) && continue
        T = interior(TV[n])
        bottom = [shelf[i, j] && kbot[i, j] > 0 ? Float64(T[i, j, kbot[i, j]]) : NaN for i in 1:Nx, j in 1:Ny]
        d = t / 86400
        gb = glorys_bottom_at(d)
        ok = shelf .& isfinite.(bottom) .& isfinite.(gb)
        push!(fdays, d)
        push!(model_area, sum(area[ok .& (bottom .< COLD)])); push!(glorys_area, sum(area[ok .& (gb .< COLD)]))
        push!(model_mean, sum(bottom[ok] .* area[ok]) / sum(area[ok])); push!(glorys_mean, sum(gb[ok] .* area[ok]) / sum(area[ok]))
        date = start_date + Second(round(Int, t))
        mon = Date(year(date), Dates.month(date), 1)
        mon in months || push!(months, mon)
        for (key, X) in ((:model, bottom), (:glorys, gb))
            S = get!(msum, (key, mon), zeros(Nx, Ny)); S[ok] .+= X[ok]
        end
        C = get!(mcount, mon, zeros(Int, Nx, Ny)); C[ok] .+= 1
        length(fdays) % 20 == 0 && @printf("  day %.0f (%s): cold pool %.0f km² model, %.0f km² GLORYS; mean bottom T %.2f / %.2f °C\n",
                                          d, Dates.format(date, "yyyy-mm-dd"), model_area[end], glorys_area[end], model_mean[end], glorys_mean[end])
    end
end

dates = [start_date + Second(round(Int, d * 86400)) for d in fdays]
println("\nMAB shelf bottom temperature, monthly means: cold-pool area (T < $(COLD) °C, 10³ km²) and mean bottom T (°C)")
@printf("  %-9s %10s %10s   %8s %8s %8s\n", "month", "model", "GLORYS", "model", "GLORYS", "diff")
for m in months
    sel = [Date(year(t), Dates.month(t), 1) == m for t in dates]
    @printf("  %-9s %10.1f %10.1f   %8.2f %8.2f %+8.2f\n", Dates.format(m, "yyyy-mm"), mean(model_area[sel]) / 1e3,
            mean(glorys_area[sel]) / 1e3, mean(model_mean[sel]), mean(glorys_mean[sel]), mean(model_mean[sel] .- glorys_mean[sel]))
end

# ---------------- figures ----------------
fig = Figure(size = (1300, 800), fontsize = 15)
ticks = [findfirst(t -> Date(year(t), Dates.month(t), 1) == m, dates) for m in months]
ax1 = Axis(fig[1, 1], ylabel = "area (10³ km²)", title = "MAB shelf cold pool: bottom temperature below $(COLD) °C")
lines!(ax1, 1:length(fdays), model_area ./ 1e3; color = :firebrick, label = "model")
lines!(ax1, 1:length(fdays), glorys_area ./ 1e3; color = :navy, label = "GLORYS")
axislegend(ax1; position = :rt)
ax2 = Axis(fig[2, 1], ylabel = "°C", title = "MAB shelf (20-200 m, west of 69.5°W) mean bottom temperature")
lines!(ax2, 1:length(fdays), model_mean; color = :firebrick, label = "model")
lines!(ax2, 1:length(fdays), glorys_mean; color = :navy, label = "GLORYS")
for ax in (ax1, ax2); ax.xticks = (ticks, Dates.format.(months, "u")); end
save(OUT * "_cold_pool_series.png", fig; px_per_unit = 1.2); println("saved ", OUT * "_cold_pool_series.png")

ib = findall(i -> BOX[1][1] - 0.2 <= λ[i] <= BOX[1][2] + 0.2, 1:Nx); jb = findall(j -> BOX[2][1] - 0.2 <= φ[j] <= BOX[2][2] + 0.2, 1:Ny)
let fig = Figure(size = (330 * length(months) + 120, 1000), fontsize = 13), hT = nothing, hD = nothing
Label(fig[0, 1:length(months)], "MAB shelf bottom temperature, monthly means: model, GLORYS, model − GLORYS (°C)", fontsize = 17)
for (k, m) in enumerate(months)
    C = mcount[m]
    M = msum[(:model, m)] ./ C; G = msum[(:glorys, m)] ./ C
    M[C .== 0] .= NaN; G[C .== 0] .= NaN
    for (row, (X, name)) in enumerate(((M, "model"), (G, "GLORYS"), (M .- G, "model − GLORYS")))
        ax = Axis(fig[row, k], title = "$(Dates.format(m, "yyyy-mm")) $name", aspect = DataAspect())
        if row < 3
            hT = heatmap!(ax, λ[ib], φ[jb], X[ib, jb]; colormap = :thermal, colorrange = (4, 24), nan_color = :gray85)
            contour!(ax, λ[ib], φ[jb], replace(X[ib, jb], NaN => 20.0); levels = [COLD], color = :white, linewidth = 1.5)
        else
            hD = heatmap!(ax, λ[ib], φ[jb], X[ib, jb]; colormap = :balance, colorrange = (-4, 4), nan_color = :gray85)
        end
    end
end
Colorbar(fig[1:2, length(months)+1], hT); Colorbar(fig[3, length(months)+1], hD)
save(OUT * "_cold_pool_maps.png", fig; px_per_unit = 1.2); println("saved ", OUT * "_cold_pool_maps.png")
end
