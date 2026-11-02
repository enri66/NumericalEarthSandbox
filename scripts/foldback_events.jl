# Follows the Gulf Stream meanders whose axis folds back (from ring_census.jl's foldback_points_<run>.csv) through a run's 3-hourly
# surface output: one panel figure (SST and surface Rossby number ζ/f, daily from 2 days before to 5 days after the first
# fold-back) and one 3-hourly video per event, in a window around it: whether the overhang is a warm crest wrapping around
# Slope Sea water (towards a warm-core ring) or a trough wrapping around it from the south (cold-core ring) shows in ζ/f.
#   MAB_TAGS, MAB_START_DATE, MAB_CELLS_PER_DEGREE, MAB_PARTITION_X, MAB_PARTITION_Y as for ring_census.jl
#   FOLDBACK_CSV=<foldback_points_<run>.csv>  FOLDBACK_OUT=<output directory>
# An event is a set of fold-back points within 150 km of each other on days at most 2 days apart.
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const PREFIXES = String.(split(ENV["MAB_TAGS"], ","))
const OUTDIR = ENV["FOLDBACK_OUT"]
mkpath(OUTDIR)
const CPD = parse(Int, get(ENV, "MAB_CELLS_PER_DEGREE", "60"))
const PX = parse(Int, get(ENV, "MAB_PARTITION_X", "7")); const PY = parse(Int, get(ENV, "MAB_PARTITION_Y", "4"))
Nx, Ny = round(Int, (12 - 4 / 12) * CPD), round(Int, (8 - 4 / 12) * CPD)
const NXR, NYR = Nx ÷ PX, Ny ÷ PY
λm = [-76 + 2 / 12 + (i - 0.5) / CPD for i in 1:Nx]
φm = [34 + 2 / 12 + (j - 0.5) / CPD for j in 1:Ny]
const Rₑ, Ω = 6.371e6, 7.292e-5
dist(a, b) = Rₑ * acos(clamp(sind(a[2]) * sind(b[2]) + cosd(a[2]) * cosd(b[2]) * cosd(a[1] - b[1]), -1, 1))

# ---------------- fold-back events ----------------
rows = [split(l, ",") for l in readlines(ENV["FOLDBACK_CSV"])[2:end]]
points = [(Date(r[1]), parse(Float64, r[2]), parse(Float64, r[3])) for r in rows]
events = Vector{Vector{eltype(points)}}()
for p in sort(points)
    k = findfirst(ev -> any(q -> 0 <= Dates.value(p[1] - q[1]) <= 2 && dist((p[2], p[3]), (q[2], q[3])) <= 150e3, ev), events)
    k === nothing ? push!(events, [p]) : push!(events[k], p)
end
@printf("%d fold-back points, %d events\n", length(points), length(events))
for (n, ev) in enumerate(events)
    @printf("  event %d: %s to %s, %d points, centre %.2fE %.2fN\n", n, ev[1][1], ev[end][1], length(ev),
            mean(p[2] for p in ev), mean(p[3] for p in ev))
end

# ---------------- 3-hourly surface output of all segments ----------------
function rank_series(prefix, names)
    files = [JLD2.jldopen("$(prefix)_rank$r.jld2", "r") for r in 0:PX*PY-1]
    iterations = sort(parse.(Int, filter(!=("serialized"), collect(keys(files[1]["timeseries/t"])))))
    times = [files[1]["timeseries/t/$it"] for it in iterations]
    fields = Dict(name => zeros(Float32, Nx, Ny, length(iterations)) for name in names)
    for (r, f) in enumerate(files), name in names, (n, it) in enumerate(iterations)
        a = f["timeseries/$name/$it"]; a = ndims(a) == 2 ? a : a[:, :, end]
        hx = (size(a, 1) - NXR) ÷ 2; hy = (size(a, 2) - NYR) ÷ 2
        bx, by = (r - 1) ÷ PY, (r - 1) % PY
        fields[name][bx*NXR+1:(bx+1)*NXR, by*NYR+1:(by+1)*NYR, n] = a[hx+1:hx+NXR, hy+1:hy+NYR]
    end
    close.(files)
    return fields, times
end
segs = [rank_series(p, ("u", "v", "T")) for p in PREFIXES]
cuts = [n < length(segs) ? segs[n][2] .< segs[n+1][2][1] : trues(length(segs[n][2])) for n in eachindex(segs)]
field(name) = cat([s[1][name][:, :, c] for (s, c) in zip(segs, cuts)]...; dims = 3)
U, V, T = field("u"), field("v"), field("T")
times = vcat([s[2][c] for (s, c) in zip(segs, cuts)]...)
stamps = [start_date + Second(round(Int, t)) for t in times]
wet = [any(T[i, j, n] != T[i, j, 1] for n in axes(T, 3)) for i in 1:Nx, j in 1:Ny]
@printf("3-hourly frames: %d, %s to %s\n", length(stamps), stamps[1], stamps[end])

# surface Rossby number ζ/f at cell corners (from u at west faces and v at south faces)
function rossby(n)
    ζ = fill(NaN32, Nx, Ny)
    dλ, dφ = deg2rad(1 / CPD), deg2rad(1 / CPD)
    for j in 2:Ny, i in 2:Nx
        (wet[i, j] && wet[i-1, j] && wet[i, j-1] && wet[i-1, j-1]) || continue
        φf = φm[j] - 0.5 / CPD
        dvdx = (V[i, j, n] - V[i-1, j, n]) / (Rₑ * cosd(φf) * dλ)
        dudy = (U[i, j, n] * cosd(φm[j]) - U[i, j-1, n] * cosd(φm[j-1])) / (Rₑ * cosd(φf) * dφ)
        ζ[i, j] = (dvdx - dudy) / (2Ω * sind(φf))
    end
    return ζ
end
sst(n) = [wet[i, j] ? T[i, j, n] : NaN32 for i in 1:Nx, j in 1:Ny]

# ---------------- one figure and one video per event ----------------
for (e, ev) in enumerate(events)
    λ0, φ0 = mean(p[2] for p in ev), mean(p[3] for p in ev)
    I = findall(λ -> abs(λ - λ0) <= 2.5, λm); J = findall(φ -> abs(φ - φ0) <= 1.75, φm)
    day0 = DateTime(ev[1][1])
    frames = findall(s -> day0 - Day(2) <= s <= day0 + Day(5), stamps)
    isempty(frames) && continue
    daily = filter(n -> hour(stamps[n]) == 0, frames)
    Tw = [sst(n)[I, J] for n in frames]
    lo, hi = quantile(filter(isfinite, reduce(vcat, vec.(Tw))), (0.02, 0.98))
    fig = Figure(size = (320 * length(daily) + 120, 720), fontsize = 14)
    Label(fig[0, 1:length(daily)], @sprintf("Fold-back event %d near %.1f°W %.1f°N (first fold-back %s): SST (top) and ζ/f (bottom)", e, -λ0, φ0, ev[1][1]), fontsize = 18)
    local hmT, hmR
    for (c, n) in enumerate(daily)
        a1 = Axis(fig[1, c], aspect = DataAspect(), title = Dates.format(stamps[n], "mm-dd HH:MM"))
        a2 = Axis(fig[2, c], aspect = DataAspect())
        hmT = heatmap!(a1, λm[I], φm[J], sst(n)[I, J]; colormap = :thermal, colorrange = (lo, hi), nan_color = :gray80)
        hmR = heatmap!(a2, λm[I], φm[J], rossby(n)[I, J]; colormap = :balance, colorrange = (-1, 1), nan_color = :gray80)
        pts = [Point2f(p[2], p[3]) for p in ev if p[1] == Date(stamps[n])]
        isempty(pts) || (scatter!(a1, pts; color = :black, markersize = 4); scatter!(a2, pts; color = :black, markersize = 4))
        c > 1 && (hideydecorations!(a1); hideydecorations!(a2))
    end
    Colorbar(fig[1, length(daily)+1], hmT; label = "SST (°C)"); Colorbar(fig[2, length(daily)+1], hmR; label = "ζ/f")
    resize_to_layout!(fig)
    save(joinpath(OUTDIR, @sprintf("foldback_event_%02d.png", e)), fig)

    k = Observable(1)
    fig = Figure(size = (1400, 620), fontsize = 16)
    title = @lift @sprintf("Fold-back event %d, %s", e, Dates.format(stamps[frames[$k]], "yyyy-mm-dd HH:MM"))
    Label(fig[0, 1:4], title, fontsize = 20)
    a1 = Axis(fig[1, 1], aspect = DataAspect(), title = "SST"); a2 = Axis(fig[1, 3], aspect = DataAspect(), title = "ζ/f")
    h1 = heatmap!(a1, λm[I], φm[J], @lift(sst(frames[$k])[I, J]); colormap = :thermal, colorrange = (lo, hi), nan_color = :gray80)
    h2 = heatmap!(a2, λm[I], φm[J], @lift(rossby(frames[$k])[I, J]); colormap = :balance, colorrange = (-1, 1), nan_color = :gray80)
    Colorbar(fig[1, 2], h1; label = "°C"); Colorbar(fig[1, 4], h2)
    record(fig, joinpath(OUTDIR, @sprintf("foldback_event_%02d.mp4", e)), eachindex(frames); framerate = 6) do i
        k[] = i
    end
end
println("saved ", OUTDIR)
