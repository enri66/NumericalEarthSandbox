# Compare runs at the open-open corners (NE, SE) against GLORYS, from their de-tided daily output: sea-level bias and
# RMS, SST RMS and maximum, and maximum surface speed in a box at each corner, plus corner maps of de-tided sea level.
#   MAB_TAGS=/abs/prefix1,/abs/prefix2 MAB_LABELS=max,product MAB_MAPDAYS=42,87 julia --project=. scripts/corner_compare.jl
using Oceananigans, NumericalEarth, Printf, Statistics, Dates, CairoMakie
using Oceananigans.Grids: λnodes, φnodes
const NCD = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "NCDatasets")]

const GDIR = get(ENV, "MAB_DATA_DIR", joinpath(homedir(), "Data", "NumericalEarth"))
const TAGS = split(ENV["MAB_TAGS"], ",")
const LABELS = split(get(ENV, "MAB_LABELS", join(basename.(TAGS), ",")), ",")
const MAPDAYS = parse.(Int, split(get(ENV, "MAB_MAPDAYS", "42,87"), ","))
const B = parse(Int, get(ENV, "MAB_CORNER_BOX", "16"))     # box size in cells
const start_date = DateTime(2019, 4, 1)

# ---- GLORYS surface fields, the average of the two daily means around a de-tided frame (centred on midnight) ----
gfile(var, date) = first(filter(x -> occursin("$(var)_GLORYSDaily_$(Dates.format(date, "yyyy-mm-dd"))T", x) &&
                                     endswith(x, "-76.0_-64.0_34.0_42.0.nc"), readdir(GDIR; join = true)))
function read_surface(var, date)
    ds = NCD.Dataset(gfile(var, date))
    lon = Float64.(ds["longitude"][:]); lat = Float64.(ds["latitude"][:])
    v = ds[var]; A = ndims(v) == 4 ? v[:, :, 1, 1] : v[:, :, 1]
    close(ds)
    return lon, lat, [ismissing(x) ? NaN : Float64(x) for x in A]
end
function glorys_surface(var, d, λ, φ)
    day0 = start_date + Day(floor(Int, d)) - Day(1)
    lon, lat, A0 = read_surface(var, day0); _, _, A1 = read_surface(var, day0 + Day(1))
    A = (A0 .+ A1) ./ 2
    i(x) = clamp(searchsortedlast(lon, x), 1, length(lon) - 1); j(y) = clamp(searchsortedlast(lat, y), 1, length(lat) - 1)
    return [begin
                a = (x - lon[i(x)]) / (lon[i(x)+1] - lon[i(x)]); b = (y - lat[j(y)]) / (lat[j(y)+1] - lat[j(y)])
                w = ((1-a)*(1-b), a*(1-b), (1-a)*b, a*b)
                v = (A[i(x), j(y)], A[i(x)+1, j(y)], A[i(x), j(y)+1], A[i(x)+1, j(y)+1])
                den = sum(w[n] for n in 1:4 if isfinite(v[n]); init = 0.0)
                den > 0 ? sum(w[n] * v[n] for n in 1:4 if isfinite(v[n]); init = 0.0) / den : NaN
            end for x in λ, y in φ]
end

# ---- model runs ----
runs = map(TAGS) do tag
    ηts = FieldTimeSeries(tag * "_eta_daily.jld2", "η"; backend = OnDisk())
    Tts = FieldTimeSeries(tag * "_surface_daily.jld2", "T"; backend = OnDisk())
    uts = FieldTimeSeries(tag * "_surface_daily.jld2", "u"; backend = OnDisk())
    vts = FieldTimeSeries(tag * "_surface_daily.jld2", "v"; backend = OnDisk())
    (; ηts, Tts, uts, vts)
end
grid = runs[1].Tts.grid; ug = grid.underlying_grid
Nx, Ny = size(ug, 1), size(ug, 2)
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center()))
bh = grid.immersed_boundary.bottom_height
wet = [bh[i, j, 1] < 0 for i in 1:Nx, j in 1:Ny]
boxes = (NE = (Nx-B+1:Nx, Ny-B+1:Ny), SE = (Nx-B+1:Nx, 1:B))
days = runs[1].Tts.times ./ 86400

surface(ts, n) = Array(interior(ts[n]))[1:Nx, 1:Ny, 1]
function speed(r, n)
    u = Array(interior(r.uts[n]))[:, 1:Ny, 1]; v = Array(interior(r.vts[n]))[1:Nx, :, 1]
    return sqrt.(((u[1:Nx, :] .+ u[2:Nx+1, :]) ./ 2) .^ 2 .+ ((v[:, 1:Ny] .+ v[:, 2:Ny+1]) ./ 2) .^ 2)
end

names = ("SSH bias (m)", "SSH RMS (m)", "SST RMS (°C)", "SST max (°C)", "max speed (m/s)")
series = Dict((c, r, m) => Float64[] for c in keys(boxes), r in eachindex(runs), m in eachindex(names))
glorys_series = Dict(c => Float64[] for c in keys(boxes))            # GLORYS SST max, for reference
for (n, d) in enumerate(days)
    gη = glorys_surface("zos", d, λ, φ); gT = glorys_surface("thetao", d, λ, φ)
    for (c, (I, J)) in pairs(boxes)
        w = [wet[i, j] && isfinite(gη[i, j]) for i in I, j in J]
        push!(glorys_series[c], maximum(gT[I, J][w]))
        for (k, r) in enumerate(runs)
            η = surface(r.ηts, n)[I, J]; T = surface(r.Tts, n)[I, J]; s = speed(r, n)[I, J]
            dη = (η .- gη[I, J])[w]; dT = (T .- gT[I, J])[w]
            for (m, x) in enumerate((mean(dη), sqrt(mean(dη .^ 2)), sqrt(mean(dT .^ 2)), maximum(T[w]), maximum(s[w])))
                push!(series[(c, k, m)], x)
            end
        end
    end
end

println("\ncorner box $(B)×$(B) cells, days $(round(Int, days[1]))-$(round(Int, days[end])): time mean (and max) per run")
for c in keys(boxes), (m, name) in enumerate(names)
    @printf("%s %-16s", c, name)
    for k in eachindex(runs)
        x = series[(c, k, m)]; @printf("   %s: %7.3f (max %7.3f)", LABELS[k], mean(x), maximum(x))
    end
    m == 4 && @printf("   GLORYS max: %7.3f", maximum(glorys_series[c]))
    println()
end

fig = Figure(size = (1800, 1100), fontsize = 14)
Label(fig[0, 1:5], "Open-open corners: " * join(LABELS, " vs ") * " (de-tided daily, $(B)×$(B)-cell boxes)", fontsize = 18)
for (row, c) in enumerate(keys(boxes)), (m, name) in enumerate(names)
    ax = Axis(fig[row, m]; title = "$c $name", xlabel = "day")
    for k in eachindex(runs); lines!(ax, days, series[(c, k, m)]; label = LABELS[k]); end
    m == 4 && lines!(ax, days, glorys_series[c]; label = "GLORYS", color = :black, linestyle = :dash)
    (row == 1 && m == 1) && axislegend(ax; position = :lb)
end
# de-tided sea level near the NE and SE corners on the map days
for (row, dmap) in enumerate(MAPDAYS)
    n = argmin(abs.(days .- dmap))
    gη = glorys_surface("zos", days[n], λ, φ)
    I = Nx-2B+1:Nx; J = 1:Ny
    for (col, (A, t)) in enumerate(((surface(runs[k].ηts, n) for k in eachindex(runs))..., gη) .=> (LABELS..., "GLORYS"))
        A = copy(A); A[.!wet] .= NaN
        ax = Axis(fig[2 + row, col]; title = @sprintf("%s: SSH day %d", t, round(Int, days[n])), aspect = DataAspect(), titlesize = 12)
        heatmap!(ax, λ[I], φ[J], A[I, J]; colormap = :balance, colorrange = (-1.2, 1.2), nan_color = :gray85)
        for (I2, J2) in values(boxes)
            lines!(ax, λ[[I2[1], I2[end], I2[end], I2[1], I2[1]]], φ[[J2[1], J2[1], J2[end], J2[end], J2[1]]]; color = :black)
        end
    end
end
out = dirname(TAGS[1]) * "/corner_compare.png"
save(out, fig; px_per_unit = 1.2); println("saved ", out)
