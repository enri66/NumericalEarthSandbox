# Gulf Stream rings and axis shape of a script-05 run against GLORYS and satellite altimetry (DUACS L4, 1/8°), from the
# de-tided daily SSH: warm- and cold-core ring censuses with tracking, the Gulf Stream path and its curvature, and a
# side-by-side SSH video. run_ring_census.sh sets these up for the year run on antares.
#   MAB_TAGS=<run prefixes of the restart segments, in time order, comma separated>  MAB_START_DATE=<the run's start>
#   MAB_CELLS_PER_DEGREE, MAB_PARTITION_X, MAB_PARTITION_Y = the run's (defaults 60, 7, 4: the year run)
#   SAT_SSH=<DUACS file (download_satellite_ssh_sst.sh)>  GLORYS_ROOT=<GLORYS12 daily archive>  RING_OUT=<output directory>
#   [RING_VIDEO=true|false]  [RING_MODEL_NAME=<label of the run in the output>]
# A warm-core ring is a closed SSH high (outermost closed contour found by lowering the level in 1 cm steps until the
# region leaves the valid area, meets a higher maximum or exceeds rmax) whose centre lies north of the Gulf Stream: outside
# the region above the Gulf Stream level that is connected to the southern edge. A cold-core ring is the same for a closed
# low outside the region below that level connected to the northern edge. The Gulf Stream level of each map is the SSH
# level along which the SSH gradient is largest. The model SSH is averaged onto the altimetry grid first, so all three are
# searched at the same resolution.
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const PREFIXES = String.(split(ENV["MAB_TAGS"], ","))
const SATFILE  = ENV["SAT_SSH"]
const OUTDIR   = ENV["RING_OUT"]
mkpath(OUTDIR)

const AMIN = 0.10      # minimum ring amplitude (m)
const RMIN = 30e3      # minimum equivalent radius of the outermost closed contour (m)
const RMAX = 250e3     # maximum equivalent radius (m)
const DH   = 0.01      # contour step (m)
const Rₑ   = 6.371e6

# ---------------- model: de-tided daily SSH of all segments ----------------
# The ranks' serialized grids hold an Open MPI communicator, which an MPICH build of MPI.jl cannot read back, so the
# layout comes from the run's configuration instead: MAB_CELLS_PER_DEGREE over 75.83-64.17W, 34.17-41.83N (script 05's
# "mab" region less 2 GLORYS cells), MAB_PARTITION_X x MAB_PARTITION_Y ranks, rank r at x-block r ÷ PY + 1, y-block r % PY + 1.
const CPD = parse(Int, get(ENV, "MAB_CELLS_PER_DEGREE", "60"))
const PX = parse(Int, get(ENV, "MAB_PARTITION_X", "7")); const PY = parse(Int, get(ENV, "MAB_PARTITION_Y", "4"))
Nx, Ny = round(Int, (12 - 4 / 12) * CPD), round(Int, (8 - 4 / 12) * CPD)
const NXR, NYR = Nx ÷ PX, Ny ÷ PY
λm = [-76 + 2 / 12 + (i - 0.5) / CPD for i in 1:Nx]
φm = [34 + 2 / 12 + (j - 0.5) / CPD for j in 1:Ny]
const VIDEO = get(ENV, "RING_VIDEO", "true") == "true"
const MODEL_NAME = get(ENV, "RING_MODEL_NAME", "model")

function rank_series(prefix)
    files = [JLD2.jldopen("$(prefix)_eta_daily_rank$r.jld2", "r") for r in 0:PX*PY-1]
    iterations = sort(parse.(Int, filter(!=("serialized"), collect(keys(files[1]["timeseries/η"])))))
    times = [files[1]["timeseries/t/$it"] for it in iterations]
    A = zeros(Nx, Ny, length(iterations))
    for (r, f) in enumerate(files), (n, it) in enumerate(iterations)
        a = f["timeseries/η/$it"]; a = ndims(a) == 2 ? a : a[:, :, end]
        hx = (size(a, 1) - NXR) ÷ 2; hy = (size(a, 2) - NYR) ÷ 2
        bx, by = (r - 1) ÷ PY, (r - 1) % PY
        A[bx*NXR+1:(bx+1)*NXR, by*NYR+1:(by+1)*NYR, n] = a[hx+1:hx+NXR, hy+1:hy+NYR]
    end
    close.(files)
    return A, times
end

# each restart segment rewrites the days after its checkpoint, so a segment is cut where the next one begins
function model_segment(prefix)
    A, t = rank_series(prefix)
    order = sortperm(t); keep = [order[k] for k in eachindex(order) if k == length(order) || t[order[k+1]] > t[order[k]]]
    return A[1:Nx, 1:Ny, keep], t[keep]
end
segs = [model_segment(p) for p in PREFIXES]
for n in 1:length(segs)-1
    cut = segs[n][2] .< segs[n+1][2][1]
    segs[n] = (segs[n][1][:, :, cut], segs[n][2][cut])
end
Hm = cat([s[1] for s in segs]...; dims = 3)
mdates = [Date(start_date + Second(round(Int, t))) for t in vcat([s[2] for s in segs]...)]
# land: the SSH of a land column never changes
wet = [any(Hm[i, j, n] != Hm[i, j, 1] for n in axes(Hm, 3)) for i in 1:Nx, j in 1:Ny]
# check of the rank layout: SSH jumps across the rank edges against those between ordinary neighbouring cells
function seam_ratio()
    d(i, j) = (wet[i, j] && wet[i+1, j]) ? abs(Hm[i+1, j, end] - Hm[i, j, end]) : NaN
    seams = filter(isfinite, [d(b * NXR, j) for b in 1:PX-1 for j in 1:Ny])
    others = filter(isfinite, [d(i, j) for i in 1:Nx-1 for j in 1:Ny if i % NXR != 0])
    return mean(seams) / mean(others)
end
@printf("land cells %.1f%%, rank-edge / interior SSH step ratio %.2f (about 1 when the layout is right)\n",
        100 * count(!, wet) / length(wet), seam_ratio())
@printf("model: %d daily maps, %s to %s, grid %d x %d (%.3f..%.3f E, %.3f..%.3f N)\n",
        length(mdates), mdates[1], mdates[end], Nx, Ny, λm[1], λm[end], φm[1], φm[end])

# ---------------- satellite ----------------
ds = NCD.Dataset(SATFILE)
λs = Float64.(ds["longitude"][:]); φs = Float64.(ds["latitude"][:])
sdates = Date.(ds["time"][:])
dims = NCD.dimnames(ds["adt"])
adt = Float64.(coalesce.(Array(ds["adt"][:, :, :]), NaN))
adt = permutedims(adt, [findfirst(==(d), dims) for d in ("longitude", "latitude", "time")])
close(ds)
Nλ, Nφ = length(λs), length(φs)
dλ, dφ = λs[2] - λs[1], φs[2] - φs[1]
@printf("satellite: %d daily maps, %s to %s, grid %d x %d\n", length(sdates), sdates[1], sdates[end], Nλ, Nφ)

# ---------------- the model on the altimetry grid ----------------
cellI(λ) = floor(Int, (λ - (λs[1] - dλ / 2)) / dλ) + 1
cellJ(φ) = floor(Int, (φ - (φs[1] - dφ / 2)) / dφ) + 1
ntot = zeros(Int, Nλ, Nφ); nwet = zeros(Int, Nλ, Nφ)
for j in 1:Ny, i in 1:Nx
    I, J = cellI(λm[i]), cellJ(φm[j]); ntot[I, J] += 1; wet[i, j] && (nwet[I, J] += 1)
end
# an altimetry cell is in the model domain when it holds (nearly) its full share of model cells, and wet when half of them are
full = floor(Int, (dλ / (λm[2] - λm[1])) * (dφ / (φm[2] - φm[1])) * 0.85)
mvalid = (ntot .>= full) .& (nwet .* 2 .>= ntot)
function binned(h)
    s = zeros(Nλ, Nφ); c = zeros(Int, Nλ, Nφ)
    for j in 1:Ny, i in 1:Nx
        (wet[i, j] && isfinite(h[i, j])) || continue
        I, J = cellI(λm[i]), cellJ(φm[j]); s[I, J] += h[i, j]; c[I, J] += 1
    end
    return [mvalid[I, J] && c[I, J] > 0 ? s[I, J] / c[I, J] : NaN for I in 1:Nλ, J in 1:Nφ]
end
dom_λ = extrema(λs[I] for I in 1:Nλ, J in 1:Nφ if mvalid[I, J])
dom_φ = extrema(φs[J] for I in 1:Nλ, J in 1:Nφ if mvalid[I, J])
# the wide satellite census: west of 52W, 32-44.5N
widevalid = [(λs[I] <= -52) && (32 <= φs[J] <= 44.5) for I in 1:Nλ, J in 1:Nφ]

# ---------------- ring detection ----------------
cellarea(J) = (Rₑ * deg2rad(dλ) * cosd(φs[J])) * (Rₑ * deg2rad(dφ))

function gradmag(h)
    g = fill(NaN, size(h))
    for J in 2:Nφ-1, I in 2:Nλ-1
        a, b, c, d = h[I+1, J], h[I-1, J], h[I, J+1], h[I, J-1]
        all(isfinite, (a, b, c, d)) || continue
        hx = (a - b) / (2Rₑ * deg2rad(dλ) * cosd(φs[J])); hy = (c - d) / (2Rₑ * deg2rad(dφ))
        g[I, J] = sqrt(hx^2 + hy^2)
    end
    return g
end

function gulf_stream_level(h)
    vals = filter(isfinite, h); g = gradmag(h)
    best, Lbest = -Inf, NaN
    for L in quantile(vals, 0.05):0.02:quantile(vals, 0.95)
        sel = [g[k] for k in eachindex(h) if isfinite(g[k]) && abs(h[k] - L) < 0.05]
        length(sel) < 20 && continue
        m = mean(sel)
        m > best && (best = m; Lbest = L)
    end
    return Lbest
end

const NBR = ((1, 0), (-1, 0), (0, 1), (0, -1))

# the region above L connected to the southernmost valid row (the Sargasso side of the Gulf Stream), or with south = false
# the region below L connected to the northernmost valid row (the Slope Sea side)
function side_region(h, L; south = true)
    inside(v) = isfinite(v) && (south ? v >= L : v < L)
    S = falses(size(h)); queue = Tuple{Int,Int}[]
    rows = south ? (1:Nφ) : (Nφ:-1:1)
    J0 = rows[findfirst(J -> any(isfinite, view(h, :, J)), rows)]
    for I in 1:Nλ
        inside(h[I, J0]) && (S[I, J0] = true; push!(queue, (I, J0)))
    end
    while !isempty(queue)
        I, J = popfirst!(queue)
        for (a, b) in NBR
            i, j = I + a, J + b
            (1 <= i <= Nλ && 1 <= j <= Nφ) || continue
            (!S[i, j] && inside(h[i, j])) || continue
            S[i, j] = true; push!(queue, (i, j))
        end
    end
    return S
end

# the Gulf Stream path: the northern edge of the Sargasso side at each longitude (NaN where the column has none)
function gulf_stream_path(h)
    S = side_region(h, gulf_stream_level(h))
    return [(J = findlast(view(S, I, :)); J === nothing ? NaN : φs[J]) for I in 1:Nλ]
end

# the connected region at or above `level` around (I, J), or nothing if it is not a closed contour of this maximum
function closed_region(h, I, J, level)
    top = h[I, J]
    seen = falses(size(h)); seen[I, J] = true
    queue = [(I, J)]; cells = Tuple{Int,Int}[]; area = 0.0
    while !isempty(queue)
        a, b = popfirst!(queue)
        push!(cells, (a, b)); area += cellarea(b)
        area > π * RMAX^2 && return nothing
        for (da, db) in NBR
            i, j = a + da, b + db
            (1 <= i <= Nλ && 1 <= j <= Nφ) || return nothing       # leaves the map
            isfinite(h[i, j]) || return nothing                    # meets land or the edge of the valid area
            (seen[i, j] || h[i, j] < level) && continue
            h[i, j] > top && return nothing                        # meets a higher maximum
            seen[i, j] = true; push!(queue, (i, j))
        end
    end
    return cells, area
end

# closed maxima of h (amplitude ≥ amin, equivalent radius ≥ RMIN) whose centres are not in `excluded`
function closed_maxima(h, excluded; amin)
    rings = NamedTuple[]
    for J in 2:Nφ-1, I in 2:Nλ-1
        c = h[I, J]
        (isfinite(c) && !excluded[I, J]) || continue
        any(isfinite(h[I+a, J+b]) && h[I+a, J+b] > c for a in -1:1, b in -1:1) && continue
        best = nothing; level = c - DH
        while true
            r = closed_region(h, I, J, level)
            r === nothing && break
            best = (r..., level); level -= DH
        end
        best === nothing && continue
        cells, area, lev = best
        amp = c - lev; rad = sqrt(area / π)
        (amp >= amin && rad >= RMIN) || continue
        push!(rings, (λ = λs[I], φ = φs[J], amp = amp, rad = rad))
    end
    return rings
end

# warm-core rings: closed highs off the Sargasso side; cold-core rings: closed lows off the Slope Sea side
function warm_core_rings(h; amin = AMIN)
    L = gulf_stream_level(h)
    return closed_maxima(h, side_region(h, L); amin), L
end
function cold_core_rings(h; amin = AMIN)
    L = gulf_stream_level(h)
    return closed_maxima(-h, side_region(h, L; south = false); amin), L
end

masked(h, valid) = [valid[I, J] && isfinite(h[I, J]) ? h[I, J] : NaN for I in 1:Nλ, J in 1:Nφ]
sat(d) = adt[:, :, findfirst(==(d), sdates)]

# ---------------- tracking ----------------
dist(p, q) = Rₑ * acos(clamp(sind(p.φ) * sind(q.φ) + cosd(p.φ) * cosd(q.φ) * cosd(p.λ - q.λ), -1, 1))

# links each day's rings to the tracks seen within the last 5 days (nearest first, at most 30 km plus 20 km per day apart)
function track(rings_by_day, dates)
    tracks = Vector{Vector{Tuple{Date,NamedTuple}}}()
    for (n, rings) in enumerate(rings_by_day)
        pairs = Tuple{Float64,Int,Int}[]
        for (k, r) in enumerate(rings), (t, tr) in enumerate(tracks)
            dlast, rlast = tr[end]
            gap = Dates.value(dates[n] - dlast)
            (1 <= gap <= 5) || continue
            d = dist(r, rlast)
            d <= 30e3 + 20e3 * gap && push!(pairs, (d, k, t))
        end
        sort!(pairs)
        usedk = falses(length(rings)); usedt = falses(length(tracks))
        for (d, k, t) in pairs
            (usedk[k] || usedt[t]) && continue
            push!(tracks[t], (dates[n], rings[k])); usedk[k] = true; usedt[t] = true
        end
        for k in eachindex(rings)
            usedk[k] || push!(tracks, [(dates[n], rings[k])])
        end
    end
    return tracks
end

# ---------------- census ----------------
# a NaN-aware box average over (2k + 1)² altimetry cells, to bring the model nearer the altimetry's effective resolution
function smoothed(h, k)
    k == 0 && return h
    out = fill(NaN, size(h))
    for J in 1:Nφ, I in 1:Nλ
        isfinite(h[I, J]) || continue
        s = 0.0; c = 0
        for b in max(1, J - k):min(Nφ, J + k), a in max(1, I - k):min(Nλ, I + k)
            isfinite(h[a, b]) && (s += h[a, b]; c += 1)
        end
        out[I, J] = s / c
    end
    return out
end

function census(name, dates, field, valid; amin = AMIN, smooth = 0, detector = warm_core_rings)
    rings_by_day = Vector{Vector{NamedTuple}}(undef, length(dates)); levels = zeros(length(dates))
    Threads.@threads for n in eachindex(dates)
        rings_by_day[n], levels[n] = detector(smoothed(masked(field(dates[n]), valid), smooth); amin)
    end
    tracks = track(rings_by_day, dates)
    @printf("%-22s %4d maps, %5d ring detections, %3d tracks (%d lasting >= 5 days), mean %.2f rings per map, GS level %.2f..%.2f m\n",
            name, length(dates), sum(length, rings_by_day), length(tracks), count(t -> length(t) >= 5, tracks),
            mean(length, rings_by_day), extrema(levels)...)
    return rings_by_day, tracks
end

# ---------------- GLORYS (daily mean SSH, 1/12°), on the altimetry grid ----------------
const GLORYS_ROOT = ENV["GLORYS_ROOT"]
function glorys_file(d)
    dir = joinpath(GLORYS_ROOT, Dates.format(d, "yyyy"), Dates.format(d, "mm"))
    isdir(dir) || return nothing
    tag = "_mean_" * Dates.format(d, "yyyymmdd") * "_"
    f = filter(x -> occursin(tag, x) && endswith(x, ".nc"), readdir(dir))
    return isempty(f) ? nothing : joinpath(dir, first(f))
end
gdates = filter(d -> glorys_file(d) !== nothing, sdates)
λg, φg, zos = NCD.Dataset(glorys_file(gdates[1])) do ds
    λa = Float64.(ds["longitude"][:]); φa = Float64.(ds["latitude"][:])
    gi = findall(x -> -80 <= x <= -50, λa); gj = findall(x -> 30 <= x <= 45, φa)
    λa[gi], φa[gj], fill(NaN, length(gi), length(gj), length(gdates))
end
for (n, d) in enumerate(gdates)     # the netCDF library is not thread-safe: one file at a time
    NCD.Dataset(glorys_file(d)) do ds
        λall = ds["longitude"][:]; φall = ds["latitude"][:]
        i1 = findfirst(==(Float32(λg[1])), λall); j1 = findfirst(==(Float32(φg[1])), φall)
        v = ds["zos"]
        NCD.dimnames(v)[1:2] == ("longitude", "latitude") || error("unexpected zos dimensions $(NCD.dimnames(v))")
        zos[:, :, n] = Float64.(coalesce.(v[i1:i1+length(λg)-1, j1:j1+length(φg)-1, 1], NaN))
    end
end
@printf("GLORYS: %d daily maps, %s to %s, %d x %d points\n", length(gdates), gdates[1], gdates[end], length(λg), length(φg))

# averages a field on its own regular grid onto the altimetry grid (cells at least half covered by finite values)
function binner(λsrc, φsrc)
    I = cellI.(λsrc); J = cellJ.(φsrc)
    return function (h)
        s = zeros(Nλ, Nφ); c = zeros(Int, Nλ, Nφ); n = zeros(Int, Nλ, Nφ)
        for j in eachindex(φsrc), i in eachindex(λsrc)
            (1 <= I[i] <= Nλ && 1 <= J[j] <= Nφ) || continue
            n[I[i], J[j]] += 1
            isfinite(h[i, j]) && (s[I[i], J[j]] += h[i, j]; c[I[i], J[j]] += 1)
        end
        return [c[a, b] > 0 && 2c[a, b] >= n[a, b] ? s[a, b] / c[a, b] : NaN for a in 1:Nλ, b in 1:Nφ]
    end
end
const glorys_binned = binner(λg, φg)
glorys(d) = glorys_binned(zos[:, :, findfirst(==(d), gdates)])

# ---------------- censuses ----------------
model_field(d) = binned(Hm[:, :, findfirst(==(d), mdates)])
cdays = [d for d in mdates if d in sdates && d in gdates]
println("sensitivity (comparison period):")
for (amin, smooth) in ((0.10, 2), (0.05, 0))
    census(@sprintf("model a≥%.2f s=%d", amin, smooth), cdays, model_field, mvalid; amin, smooth)
end
census("satellite box a≥0.05", cdays, sat, mvalid; amin = 0.05)
census("GLORYS box a≥0.05", cdays, glorys, mvalid; amin = 0.05)
println("reference:")
m_rings, m_tracks = census("model", cdays, model_field, mvalid)
gd_rings, gd_tracks = census("GLORYS, model box", cdays, glorys, mvalid)
sd_rings, sd_tracks = census("satellite, model box", cdays, sat, mvalid)
yday = [d for d in sdates if Date(2019, 4, 1) <= d <= Date(2020, 3, 31) && d in gdates]
gw_rings, gw_tracks = census("GLORYS, 80-52W", yday, glorys, widevalid)
sw_rings, sw_tracks = census("satellite, 80-52W", yday, sat, widevalid)

# the state the model started from
println("\nrings in the model box on 2019-04-01 (lon, lat, amplitude m, radius km):")
for (name, f) in (("satellite", sat), ("GLORYS", glorys))
    r, _ = warm_core_rings(masked(f(Date(2019, 4, 1)), mvalid))
    println(rpad("  " * name, 12), [(round(x.λ; digits = 2), round(x.φ; digits = 2), round(x.amp; digits = 2), round(Int, x.rad / 1e3)) for x in r])
end

# GLORYS against the satellite, day by day: rings of one with a ring of the other within 100 km
function matched(a, b)
    hits = 0; total = 0
    for n in eachindex(a), r in a[n]
        total += 1
        any(q -> dist(r, q) <= 100e3, b[n]) && (hits += 1)
    end
    return hits, total
end
@printf("GLORYS rings in the model box with a satellite ring within 100 km: %d of %d\n", matched(gd_rings, sd_rings)...)
@printf("satellite rings in the model box with a GLORYS ring within 100 km: %d of %d\n", matched(sd_rings, gd_rings)...)
@printf("satellite rings 80-52W (whole year) with a GLORYS ring within 100 km: %d of %d\n", matched(sw_rings, gw_rings)...)

println("\ncold-core rings (closed lows off the Slope Sea side, amplitude ≥ 10 cm):")
mc_rings, _ = census("model", cdays, model_field, mvalid; detector = cold_core_rings)
gc_rings, _ = census("GLORYS, model box", cdays, glorys, mvalid; detector = cold_core_rings)
sc_rings, _ = census("satellite, model box", cdays, sat, mvalid; detector = cold_core_rings)
census("GLORYS, 80-52W", yday, glorys, widevalid; detector = cold_core_rings)
census("satellite, 80-52W", yday, sat, widevalid; detector = cold_core_rings)

# how far the Gulf Stream meanders: the spread of its path latitude at each longitude over the comparison period
paths = Dict(name => reduce(hcat, [gulf_stream_path(masked(f(d), mvalid)) for d in cdays])
             for (name, f) in (("model", model_field), ("GLORYS", glorys), ("satellite", sat)))
inbox = [dom_λ[1] + 0.5 <= λs[I] <= dom_λ[2] - 0.5 for I in 1:Nλ]
pathstats(P) = ([mean(filter(isfinite, P[I, :])) for I in 1:Nλ], [std(filter(isfinite, P[I, :])) for I in 1:Nλ])
println("\nGulf Stream path latitude over the comparison period, by 2° of longitude (mean / standard deviation, °N):")
for name in ("model", "GLORYS", "satellite")
    m, s = pathstats(paths[name])
    bins = [(lo, findall(I -> inbox[I] && lo <= λs[I] < lo + 2, 1:Nλ)) for lo in -76:2:-66]
    println(rpad("  " * name, 12), join([@sprintf("%d..%dE %.2f/%.2f", lo, lo + 2, mean(m[ix]), mean(s[ix])) for (lo, ix) in bins if !isempty(ix)], "   "))
end

fig = Figure(size = (1100, 650), fontsize = 16)
axm = Axis(fig[1, 1], ylabel = "mean path latitude (°N)", title = "Gulf Stream path, $(cdays[1]) to $(cdays[end])")
axs = Axis(fig[2, 1], ylabel = "path latitude std (°)", xlabel = "longitude")
for (name, c) in (("model", :red), ("GLORYS", :blue), ("satellite", :black))
    m, s = pathstats(paths[name])
    lines!(axm, λs[inbox], m[inbox], color = c, linewidth = 2, label = name)
    lines!(axs, λs[inbox], s[inbox], color = c, linewidth = 2)
end
axislegend(axm, position = :lt)
save(joinpath(OUTDIR, "gulf_stream_path.png"), fig)

# ---------------- curvature along the Gulf Stream axis ----------------
# The axis is the longest contour at the Gulf Stream level, traced on the altimetry grid inside the model box, resampled every
# 10 km and smoothed with a 25 km Gaussian (the same for all three, so the grid scale does not count as curvature); points within
# 0.5° of the box edges (the sponge) are left out.
const CT = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "Contour")]   # loaded by Makie
const DS, SIGMA = 10e3, 25e3

# fills land and the outside of the valid area with the average of known neighbours, pass by pass, spreading outwards
function extended(h)
    f = copy(h); known = isfinite.(h)
    for pass in 1:60
        all(known) && break
        g = copy(f); newknown = copy(known)
        for J in 1:Nφ, I in 1:Nλ
            known[I, J] && continue
            s = 0.0; c = 0
            for (a, b) in NBR
                i, j = I + a, J + b
                (1 <= i <= Nλ && 1 <= j <= Nφ && known[i, j]) && (s += f[i, j]; c += 1)
            end
            c > 0 && (g[I, J] = s / c; newknown[I, J] = true)
        end
        f = g; known = newknown
    end
    return ifelse.(known, f, mean(filter(isfinite, h)))
end

# curvature κ along the axis, the positions, and the eastward component of its direction (negative where the axis folds back)
function axis_curvature(h; sigma = SIGMA)
    L = gulf_stream_level(h)
    z = extended(h)                      # the contour leaves the box instead of running along its edge
    best, bestlen = nothing, 0.0
    for line in CT.lines(CT.contour(λs, φs, z, L))
        x, y = CT.coordinates(line)
        len = sum(dist((λ = x[n], φ = y[n]), (λ = x[n+1], φ = y[n+1])) for n in 1:length(x)-1; init = 0.0)
        len > bestlen && (best = (x, y); bestlen = len)
    end
    best === nothing && return Float64[], Float64[], Float64[], Float64[]
    λl, φl = best
    λl[end] < λl[1] && (λl = reverse(λl); φl = reverse(φl))    # downstream = eastward overall
    φ0 = mean(φl)
    X = Rₑ * cosd(φ0) .* deg2rad.(λl); Y = Rₑ .* deg2rad.(φl)
    s = vcat(0.0, cumsum(hypot.(diff(X), diff(Y))))
    grid_s = 0.0:DS:s[end]
    interp(v) = [(k = clamp(searchsortedlast(s, q), 1, length(s) - 1); w = (q - s[k]) / max(s[k+1] - s[k], eps());
                  (1 - w) * v[k] + w * v[k+1]) for q in grid_s]
    Xs, Ys = interp(X), interp(Y)
    n = length(Xs); w = ceil(Int, 3sigma / DS)
    kernel = [exp(-(m * DS)^2 / (2sigma^2)) for m in -w:w]
    smooth(v) = [(r = max(1, i - w):min(n, i + w); sum(v[r] .* kernel[r .- i .+ w .+ 1]) / sum(kernel[r .- i .+ w .+ 1])) for i in 1:n]
    Xs, Ys = smooth(Xs), smooth(Ys)
    κ = Float64[]; λk = Float64[]; φk = Float64[]; east = Float64[]
    for i in 2:n-1
        x1, y1 = (Xs[i+1] - Xs[i-1]) / 2DS, (Ys[i+1] - Ys[i-1]) / 2DS
        x2, y2 = (Xs[i+1] - 2Xs[i] + Xs[i-1]) / DS^2, (Ys[i+1] - 2Ys[i] + Ys[i-1]) / DS^2
        λp = rad2deg(Xs[i] / (Rₑ * cosd(φ0))); φp = rad2deg(Ys[i] / Rₑ)
        (dom_λ[1] + 0.5 <= λp <= dom_λ[2] - 0.5 && dom_φ[1] + 0.5 <= φp <= dom_φ[2] - 0.5) || continue
        push!(κ, (x1 * y2 - y1 * x2) / (x1^2 + y1^2)^1.5); push!(λk, λp); push!(φk, φp); push!(east, x1 / hypot(x1, y1))
    end
    return κ, λk, φk, east
end

println("\ncurvature along the Gulf Stream axis (inside the box less 0.5°), radius of curvature R = 1/|κ|:")
curv = Dict{String,Vector{Float64}}()
for (name, f, sm) in ((MODEL_NAME, model_field, 0), (MODEL_NAME * " s=2", model_field, 2), ("GLORYS", glorys, 0), ("satellite", sat, 0))
    κs = reduce(vcat, [axis_curvature(smoothed(masked(f(d), mvalid), sm))[1] for d in cdays])
    curv[name] = κs
    a = abs.(κs) .* 1e3     # 1/km
    @printf("  %-16s %6d points: median |κ| %.4f /km, 90th pct %.4f /km; share of the axis with R < 50 km %4.1f%%, < 100 km %4.1f%%, < 200 km %4.1f%%\n",
            name, length(a), median(a), quantile(a, 0.9), 100 * mean(a .> 1 / 50), 100 * mean(a .> 1 / 100), 100 * mean(a .> 1 / 200))
end

# meander scale (60 km smoothing): curvature, and how much of the axis folds back (runs westward) or runs within 30° of north-south,
# per map: a crest or trough that overhangs is the step before a ring pinches off
println("\nmeander-scale axis shape (60 km smoothing), inside the box less 0.5°:")
for (name, f) in ((MODEL_NAME, model_field), ("GLORYS", glorys), ("satellite", sat))
    res = [axis_curvature(masked(f(d), mvalid); sigma = 60e3) for d in cdays]
    a = abs.(reduce(vcat, [r[1] for r in res])) .* 1e3
    fold = [isempty(r[4]) ? NaN : mean(r[4] .< 0) for r in res]
    steep = [isempty(r[4]) ? NaN : mean(abs.(r[4]) .< sind(30)) for r in res]
    @printf("  %-12s share of the axis with R < 100 km %4.1f%%, < 200 km %4.1f%%; folding back: %4.1f%% of the axis, on %d of %d maps; within 30° of north-south: %4.1f%%\n",
            name, 100 * mean(a .> 1 / 100), 100 * mean(a .> 1 / 200), 100 * mean(filter(isfinite, fold)),
            count(x -> isfinite(x) && x > 0, fold), length(fold), 100 * mean(filter(isfinite, steep)))
    # the axis points that fold back, for following those meanders in the 3-hourly output (foldback_events.jl); whether a
    # fold-back belongs to a crest or a trough is not given by the sign of κ (an overhang curves both ways), so they are all kept
    open(joinpath(OUTDIR, "foldback_points_$(replace(name, ' ' => '_')).csv"), "w") do io
        println(io, "date,lon,lat,curvature_per_km")
        for (d, r) in zip(cdays, res), k in eachindex(r[1])
            r[4][k] < 0 && @printf(io, "%s,%.3f,%.3f,%.5f\n", d, r[2][k], r[3][k], 1e3 * r[1][k])
        end
    end
end
fig = Figure(size = (900, 500), fontsize = 16)
ax = Axis(fig[1, 1], xlabel = "radius of curvature R (km)", ylabel = "share of the axis with a radius below R",
          xscale = log10, title = "Gulf Stream axis curvature, $(cdays[1]) to $(cdays[end])")
for (name, c) in ((MODEL_NAME, :red), (MODEL_NAME * " s=2", (:red, 0.4)), ("GLORYS", :blue), ("satellite", :black))
    R = sort(1 ./ abs.(curv[name]) ./ 1e3)
    lines!(ax, R, (1:length(R)) ./ length(R), color = c, linewidth = 2, label = name)
end
xlims!(ax, 20, 5000)
axislegend(ax, position = :lt)
save(joinpath(OUTDIR, "gulf_stream_curvature.png"), fig)

# ---------------- outputs ----------------
open(joinpath(OUTDIR, "ring_tracks.csv"), "w") do io
    println(io, "source,track,first_date,last_date,detections,first_lon,first_lat,last_lon,last_lat,max_amp_m,mean_radius_km")
    for (src, tracks) in (("model", m_tracks), ("glorys_model_box", gd_tracks), ("satellite_model_box", sd_tracks),
                          ("glorys_wide", gw_tracks), ("satellite_wide", sw_tracks))
        for (k, tr) in enumerate(tracks)
            f, l = tr[1][2], tr[end][2]
            @printf(io, "%s,%d,%s,%s,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.1f\n", src, k, tr[1][1], tr[end][1], length(tr),
                    f.λ, f.φ, l.λ, l.φ, maximum(x[2].amp for x in tr), mean(x[2].rad for x in tr) / 1e3)
        end
    end
end

# where tracks lasting >= 5 days first appear, and how many cross the model's east edge going west
function formation_summary(name, tracks)
    long = filter(t -> length(t) >= 5, tracks)
    println("\n$name 2019-04-01..2020-03-31, 80-52W, tracks lasting >= 5 days: ", length(long))
    println("  first seen by 4° of longitude from 80W: ", [count(t -> lo <= t[1][2].λ < lo + 4, long) for lo in -80:4:-56])
    east = dom_λ[2]
    @printf("  crossing %.2fE (the model's east edge) westward: %d\n", east,
            count(t -> any(x -> x[2].λ > east, t) && t[end][2].λ < east - 0.5, long))
    @printf("  first seen inside the model box: %d\n",
            count(t -> dom_λ[1] <= t[1][2].λ <= dom_λ[2] && dom_φ[1] <= t[1][2].φ <= dom_φ[2], long))
    @printf("  median lifetime %.0f days, median amplitude %.2f m\n", median(Dates.value(t[end][1] - t[1][1]) + 1 for t in long),
            median(maximum(x[2].amp for x in t) for t in long))
end
formation_summary("satellite", sw_tracks)
formation_summary("GLORYS", gw_tracks)
for (label, rb) in (("model", m_rings), ("GLORYS, model box", gd_rings), ("satellite, model box", sd_rings))
    @printf("%-22s ring-days in the comparison period: %d\n", label, sum(length, rb))
end

# daily counts
fig = Figure(size = (1100, 440), fontsize = 16)
ax = Axis(fig[1, 1], ylabel = "warm-core rings", xlabel = "days since 2019-04-01",
          title = "Warm-core rings north of the Gulf Stream (closed SSH highs, amplitude ≥ 10 cm)")
dn(ds) = Dates.value.(ds .- Date(2019, 4, 1))
lines!(ax, dn(yday), length.(sw_rings), color = (:gray, 0.6), label = "satellite, 80-52W")
lines!(ax, dn(yday), length.(gw_rings), color = (:dodgerblue, 0.5), label = "GLORYS, 80-52W")
lines!(ax, dn(cdays), length.(sd_rings), color = :black, linewidth = 2, label = "satellite, model box")
lines!(ax, dn(cdays), length.(gd_rings), color = :blue, linewidth = 2, label = "GLORYS, model box")
lines!(ax, dn(cdays), length.(m_rings), color = :red, linewidth = 2, label = "model")
axislegend(ax, position = :rt)
save(joinpath(OUTDIR, "ring_counts.png"), fig)

# side-by-side video (RING_VIDEO=false skips it): model SSH (full resolution), GLORYS and satellite ADT, each minus its mean
# over the model's ocean
if VIDEO
circle(r) = (θ = range(0, 2π; length = 40); [Point2f(r.λ + r.rad / (111e3 * cosd(r.φ)) * cos(t), r.φ + r.rad / 111e3 * sin(t)) for t in θ])
ringlines(rings) = isempty(rings) ? [Point2f(NaN, NaN)] : reduce(vcat, [vcat(circle(r), [Point2f(NaN, NaN)]) for r in rings])
k = Observable(1)
mframe = @lift begin
    h = [wet[i, j] ? Hm[i, j, findfirst(==(cdays[$k]), mdates)] : NaN for i in 1:Nx, j in 1:Ny]
    h .- mean(filter(isfinite, h))
end
demeaned(s) = s .- mean(filter(isfinite, masked(s, mvalid)))
gframe = @lift begin
    g = zos[:, :, findfirst(==(cdays[$k]), gdates)]
    g .- mean(filter(isfinite, masked(glorys_binned(g), mvalid)))
end
sframe = @lift demeaned(sat(cdays[$k]))
yk(k) = findfirst(==(cdays[k]), yday)
mrings = @lift ringlines(m_rings[$k])
grings = @lift ringlines(gd_rings[$k]); gwr = @lift ringlines(gw_rings[yk($k)])
srings = @lift ringlines(sd_rings[$k]); swr = @lift ringlines(sw_rings[yk($k)])
title = @lift @sprintf("De-tided SSH and warm-core rings  %s (day %d)", cdays[$k], Dates.value(cdays[$k] - Date(2019, 4, 1)))
fig = Figure(size = (2400, 700), fontsize = 18)
Label(fig[0, 1:4], title, fontsize = 22)
lims = ((-80, -52), (32, 44.5))
a1 = Axis(fig[1, 1], aspect = DataAspect(), title = "model (1/60°)", limits = lims)
a2 = Axis(fig[1, 2], aspect = DataAspect(), title = "GLORYS (1/12°)", limits = lims)
a3 = Axis(fig[1, 3], aspect = DataAspect(), title = "satellite altimetry (DUACS L4, 1/8°)", limits = lims)
hm = heatmap!(a1, λm, φm, mframe; colormap = :balance, colorrange = (-1.2, 1.2), nan_color = :gray80)
heatmap!(a2, λg, φg, gframe; colormap = :balance, colorrange = (-1.2, 1.2), nan_color = :gray80)
heatmap!(a3, λs, φs, sframe; colormap = :balance, colorrange = (-1.2, 1.2), nan_color = :gray80)
box = Point2f[(dom_λ[1], dom_φ[1]), (dom_λ[2], dom_φ[1]), (dom_λ[2], dom_φ[2]), (dom_λ[1], dom_φ[2]), (dom_λ[1], dom_φ[1])]
for a in (a1, a2, a3)
    lines!(a, box, color = :black, linewidth = 1)
end
lines!(a1, mrings, color = :black, linewidth = 2.5)
lines!(a2, gwr, color = (:black, 0.5), linewidth = 1.5); lines!(a2, grings, color = :black, linewidth = 2.5)
lines!(a3, swr, color = (:black, 0.5), linewidth = 1.5); lines!(a3, srings, color = :black, linewidth = 2.5)
Colorbar(fig[1, 4], hm; label = "SSH minus its mean over the model's ocean (m)")
resize_to_layout!(fig)
record(fig, joinpath(OUTDIR, "ssh_rings_model_glorys_satellite.mp4"), eachindex(cdays); framerate = 4) do i
    k[] = i
end
end
println("saved ", OUTDIR)
