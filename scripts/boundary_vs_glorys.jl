# Compare a run of 04_mab_glorys_tides_reservoirs.jl with GLORYS, as a function of distance from the open
# boundaries. Daily means of the model's 3-hourly surface output and hourly η are compared with GLORYS daily
# means (bilinear to the model cell centers, linear in depth to the top cell center); the de-tided daily
# volume output is compared in 3D for T and S. Usage:
#   MAB_TAG=mab_obc14 julia --project=. scripts/boundary_vs_glorys.jl
using Oceananigans, NumericalEarth, Printf, Statistics, Dates, CairoMakie
const NCD = NumericalEarth.NCDatasets

const OUT  = joinpath(homedir(), "Data", "mab_glorys_obc")
const GDIR = get(ENV, "MAB_DATA_DIR", joinpath(homedir(), "Data", "NumericalEarth"))
const TAG  = get(ENV, "MAB_TAG", "mab_obc14")
const start_date = DateTime(2019, 4, 1)
const MAXDAY = parse(Int, get(ENV, "MAB_MAXDAY", "100000"))   # limit the number of days analysed

surf = joinpath(OUT, TAG * ".jld2")
T3 = FieldTimeSeries(surf, "T"); S3 = FieldTimeSeries(surf, "S"); U3 = FieldTimeSeries(surf, "u"); V3 = FieldTimeSeries(surf, "v")
Eh = FieldTimeSeries(joinpath(OUT, TAG * "_eta.jld2"), "η")
grid = T3.grid; ug = grid.underlying_grid
Nx, Ny, Nz = size(ug)
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center())); zc = collect(znodes(ug, Center()))
bottom = Array(interior(grid.immersed_boundary.bottom_height))[:, :, 1]
wet3 = [zc[k] > bottom[i, j] for i in 1:Nx, j in 1:Ny, k in 1:Nz]
wet = wet3[:, :, Nz]

# distance (in cells) from the nearest open boundary, and which side it is
dist = [min(i - 1, Nx - i, j - 1, Ny - j) for i in 1:Nx, j in 1:Ny]
side = [argmin((i - 1, Nx - i, j - 1, Ny - j)) for i in 1:Nx, j in 1:Ny]   # 1 west, 2 east, 3 south, 4 north
bands = [("edge (0)", 0:0), ("1-2", 1:2), ("3-5", 3:5), ("6-10", 6:10), ("11-20", 11:20), ("interior (>20)", 21:10^6)]
sidenames = ("west", "east", "south", "north")

# ---------------- GLORYS on the model grid ----------------
gfile(var, d) = first(filter(f -> occursin("$(var)_GLORYSDaily_$(Dates.format(d, "yyyy-mm-dd"))T", f) && endswith(f, "-76.0_-64.0_34.0_42.0.nc"), readdir(GDIR; join = true)))

function bilinear(A, lon, lat, x, y)
    i = clamp(searchsortedlast(lon, x), 1, length(lon) - 1); j = clamp(searchsortedlast(lat, y), 1, length(lat) - 1)
    a = (x - lon[i]) / (lon[i+1] - lon[i]); b = (y - lat[j]) / (lat[j+1] - lat[j])
    w = ((1 - a) * (1 - b), a * (1 - b), (1 - a) * b, a * b); v = (A[i, j], A[i+1, j], A[i, j+1], A[i+1, j+1])
    num = sum(w[n] * v[n] for n in 1:4 if !ismissing(v[n]) && isfinite(v[n]); init = 0.0)
    den = sum(w[n] for n in 1:4 if !ismissing(v[n]) && isfinite(v[n]); init = 0.0)
    return den > 0 ? num / den : NaN
end

function glorys_on_grid(var, d; levels = [Nz])
    ds = NCD.Dataset(gfile(var, d)); lon = Float64.(ds["longitude"][:]); lat = Float64.(ds["latitude"][:])
    dep = Float64.(ds["depth"][:]); A = ds[var][:, :, :, 1]; close(ds)
    out = fill(NaN, Nx, Ny, length(levels))
    for (n, k) in enumerate(levels)
        depth = -zc[k]
        m = clamp(searchsortedlast(dep, depth), 1, length(dep) - 1); f = clamp((depth - dep[m]) / (dep[m+1] - dep[m]), 0, 1)
        L1 = A[:, :, m]; L2 = A[:, :, m+1]
        for i in 1:Nx, j in 1:Ny
            wet3[i, j, k] || continue
            a = bilinear(L1, lon, lat, λ[i], φ[j]); b = bilinear(L2, lon, lat, λ[i], φ[j])
            out[i, j, n] = isnan(b) ? a : (1 - f) * a + f * b
        end
    end
    return out
end
function glorys_zos(d)
    ds = NCD.Dataset(gfile("zos", d)); lon = Float64.(ds["longitude"][:]); lat = Float64.(ds["latitude"][:]); A = ds["zos"][:, :, 1]; close(ds)
    return [wet[i, j] ? bilinear(A, lon, lat, λ[i], φ[j]) : NaN for i in 1:Nx, j in 1:Ny]
end

# ---------------- model daily means ----------------
tday(fts) = fts.times ./ 86400
function daily_mean(fts, d; face = nothing)
    t = tday(fts); idx = [n for n in eachindex(t) if d <= t[n] < d + 1]
    isempty(idx) && return nothing
    acc = zeros(Nx, Ny)
    for n in idx
        A = Array(interior(fts[n]))[:, :, 1]
        acc .+= face === :x ? (A[1:Nx, :] .+ A[2:Nx+1, :]) ./ 2 : face === :y ? (A[:, 1:Ny] .+ A[:, 2:Ny+1]) ./ 2 : A[1:Nx, 1:Ny]
    end
    return acc ./ length(idx)
end
anom(A) = (m = mean(A[i, j] for i in 1:Nx, j in 1:Ny if wet[i, j] && isfinite(A[i, j])); A .- m)

ndays = min(floor(Int, maximum(tday(T3))), MAXDAY)
vars = ("SST", "SSS", "η", "u", "v")
stats = Dict{Tuple{String, Int, Int}, Tuple{Float64, Float64}}()   # (var, band, day) -> (bias, rms)
sidestats = Dict{Tuple{String, Int, Int}, Tuple{Float64, Float64}}()  # (var, side, day) for the edge + 1-2 bands
last_maps = nothing
for d in 0:ndays-1
    date = start_date + Day(d)
    m = Dict("SST" => daily_mean(T3, d), "SSS" => daily_mean(S3, d), "η" => anom(daily_mean(Eh, d)),
             "u" => daily_mean(U3, d; face = :x), "v" => daily_mean(V3, d; face = :y))
    g = Dict("SST" => glorys_on_grid("thetao", date)[:, :, 1], "SSS" => glorys_on_grid("so", date)[:, :, 1],
             "η" => anom(glorys_zos(date)), "u" => glorys_on_grid("uo", date)[:, :, 1], "v" => glorys_on_grid("vo", date)[:, :, 1])
    for v in vars
        D = m[v] .- g[v]
        for (b, (_, r)) in enumerate(bands)
            x = [D[i, j] for i in 1:Nx, j in 1:Ny if wet[i, j] && dist[i, j] in r && isfinite(D[i, j])]
            stats[(v, b, d)] = isempty(x) ? (NaN, NaN) : (mean(x), sqrt(mean(abs2, x)))
        end
        for s in 1:4
            x = [D[i, j] for i in 1:Nx, j in 1:Ny if wet[i, j] && dist[i, j] <= 2 && side[i, j] == s && isfinite(D[i, j])]
            sidestats[(v, s, d)] = isempty(x) ? (NaN, NaN) : (mean(x), sqrt(mean(abs2, x)))
        end
    end
    global last_maps = (d, m, g)
    @printf("day %2d done\n", d)
end

units = Dict("SST" => "°C", "SSS" => "psu", "η" => "m", "u" => "m/s", "v" => "m/s")
println("\n==== RMS model − GLORYS (daily means) by distance from the open boundary ====")
for v in vars
    @printf("\n%s (%s)   day: %s\n", v, units[v], join([@sprintf("%6d", d) for d in 0:ndays-1]))
    for (b, (name, _)) in enumerate(bands)
        @printf("  %-16s %s\n", name, join([@sprintf("%6.3f", stats[(v, b, d)][2]) for d in 0:ndays-1]))
    end
end
println("\n==== bias (model − GLORYS) within 2 cells of each open side, first and last day ====")
for v in vars
    @printf("%-4s", v)
    for s in 1:4
        @printf("  %-5s %+7.3f → %+7.3f", sidenames[s], sidestats[(v, s, 0)][1], sidestats[(v, s, ndays - 1)][1])
    end
    println()
end

# ---------------- 3D T and S from the de-tided daily volume output ----------------
vol = joinpath(OUT, TAG * "_volume_daily.jld2")
if isfile(vol)
    TV = FieldTimeSeries(vol, "T"); SV = FieldTimeSeries(vol, "S")
    Δz = diff(collect(znodes(ug, Face())))
    println("\n==== 3D RMS (depth-weighted over wet cells) of the de-tided daily T and S vs GLORYS, by band ====")
    @printf("%-8s %-6s %s\n", "day", "var", join([@sprintf("%-15s", name) for (name, _) in bands]))
    for n in eachindex(TV.times)
        d = TV.times[n] / 86400
        d <= MAXDAY || continue
        # a de-tided frame is centred on midnight, between two GLORYS daily means (centred on noon)
        day0 = start_date + Day(floor(Int, d)) - Day(1)
        for (var, fts, gv) in (("T", TV, "thetao"), ("S", SV, "so"))
            M = Array(interior(fts[n])); G = (glorys_on_grid(gv, day0; levels = 1:Nz) .+ glorys_on_grid(gv, day0 + Day(1); levels = 1:Nz)) ./ 2
            row = map(bands) do (_, r)
                num = den = 0.0
                for i in 1:Nx, j in 1:Ny, k in 1:Nz
                    (wet3[i, j, k] && dist[i, j] in r && isfinite(G[i, j, k])) || continue
                    num += Δz[k] * (M[i, j, k] - G[i, j, k])^2; den += Δz[k]
                end
                den > 0 ? sqrt(num / den) : NaN
            end
            @printf("%-8.1f %-6s %s\n", d, var, join([@sprintf("%-15.3f", x) for x in row]))
        end
    end
end

# ---------------- figure ----------------
d, m, g = last_maps
fig = Figure(size = (2000, 1400), fontsize = 18)
Label(fig[0, 1:6], "$(TAG): model vs GLORYS, day $(d) daily means (top), RMS difference by distance from the open boundary (bottom)", fontsize = 22)
mask(A) = [wet[i, j] ? A[i, j] : NaN for i in 1:Nx, j in 1:Ny]
for (c, (v, cr, cm)) in enumerate((("SST", (8, 26), :thermal), ("η", (-0.8, 0.8), :balance)))
    base = 3 * (c - 1)
    for (k, (A, t)) in enumerate(((m[v], "model"), (g[v], "GLORYS"), (m[v] .- g[v], "model − GLORYS")))
        ax = Axis(fig[1, base + k], title = "$(v): $(t)", aspect = DataAspect())
        rng = k == 3 ? (v == "SST" ? (-4, 4) : (-0.3, 0.3)) : cr
        hm = heatmap!(ax, λ, φ, mask(A); colormap = k == 3 ? :balance : cm, colorrange = rng)
        k == 3 && Colorbar(fig[2, base + k], hm; vertical = false, flipaxis = false)
        k == 1 && Colorbar(fig[2, base + 1:base + 2], hm; vertical = false, flipaxis = false)
    end
end
colors = cgrad(:viridis, length(bands); categorical = true)
for (c, v) in enumerate(vars)
    ax = Axis(fig[3, c + (c > 3 ? 1 : 0)], title = "RMS $(v) ($(units[v]))", xlabel = "day")
    for (b, (name, _)) in enumerate(bands)
        lines!(ax, 0:ndays-1, [stats[(v, b, dd)][2] for dd in 0:ndays-1]; color = colors[b], linewidth = 2.5, label = name)
    end
    c == 1 && axislegend(ax, position = :lt, labelsize = 13)
end
out = joinpath(OUT, TAG * "_boundary_vs_glorys.png")
save(out, fig; px_per_unit = 1.2)
println("\nsaved ", out)
