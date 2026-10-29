# Closes the nitrogen budget of the shelf from a run with MAB_BGC_NBUDGET=true (bgc/nbudget.jl) and the multi-G sediment:
#
#   dW/dt = R − L + (return from the sediment) − F + (sponge) + (ScaleNegativeTracers) + X
#
#   W        N in the shelf water (hourly, in-run)            R   river input from the run's log (Q × (NO₃ + NH₄ + DOM))
#   L        lateral flux out through the shelf faces (in-run, every step, with the tides)
#   F        sinking flux of sPOM + bPOM into the sediment    return = λs Ns + λf Nf of the sediment (daily means)
#   sponge   N the BGC sponges add in the shelf columns       X is what is left, to be zero apart from the numerics of the
#   ScaleNegativeTracers: N added by the negative-tracer scaling        flux across the boundary and the sediment's own error
# Per day, kmol N/day, from the hourly in-run integrals (summed over the ranks). Usage:
#   BGC_TAG=/t0/.../bgc_feb19_woa_lim_mg30_nb BGC_LOG=/t0/.../bgc_feb19_woa_lim_mg30_nb_<job>.log julia --project=<env> bgc/nbudget_analysis.jl
using Printf, Statistics, Dates
using Oceananigans
using Oceananigans.Grids: znodes
using NumericalEarth.DataWrangling: jldopen

const TAG   = ENV["BGC_TAG"]
const LOG   = get(ENV, "BGC_LOG", "")
const START = DateTime(get(ENV, "MAB_START_DATE", "2019-02-01"))
const SHELF = parse(Float64, get(ENV, "BGC_SHELF_DEPTH", "100"))
const OUT   = get(ENV, "BGC_OUT", "/t0/workdir/enrique/runs/analysis/nbudget_" * basename(TAG))
const H, Nz = 7, parse(Int, get(ENV, "MAB_NZ", "50"))
const R_EARTH = 6.371e6
const CPD = parse(Int, get(ENV, "MAB_CELLS_PER_DEGREE", "12"))
const φ_south = 34 + 2 / 12
mkpath(OUT)
io = open(joinpath(OUT, "summary.txt"), "w")
function say(args...); println(args...); flush(stdout); println(io, args...); flush(io); end

# in-run files
files = sort(filter(f -> occursin(Regex("^" * basename(TAG) * "_nbudget_rank\\d+\\.txt\$"), basename(f)), readdir(dirname(TAG); join = true)))
isempty(files) && error("no $(basename(TAG))_nbudget_rank*.txt")
rows = map(files) do f
    lines = filter(l -> !startswith(l, "#") && !isempty(strip(l)), readlines(f))
    reduce(hcat, [parse.(Float64, split(l)) for l in lines])'
end
nt = minimum(size(r, 1) for r in rows)
t = rows[1][1:nt, 1]
all(all(r[1:nt, 1] .== t) for r in rows) || error("the ranks' times differ")
W, Lat, Spo, Neg = (sum(r[1:nt, c] for r in rows) for c in 2:5)
say("$(length(files)) ranks, $nt hourly lines, day $(round(t[1] / 86400, digits = 2)) to $(round(t[end] / 86400, digits = 2))")

# rivers from the log
function parse_river(log)
    out = Tuple{Float64, Float64}[]; tt = NaN
    for line in eachline(log)
        m = match(r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2})\s", line)
        m === nothing || (tt = Dates.value(DateTime(m.captures[1], dateformat"yyyy-mm-dd HH:MM") - START) / 1000)
        m2 = match(r"rivers: Q = ([\d.]+) m³/s.*DOM ([\d.]+), NH₄ ([\d.]+), NO₃ ([\d.]+)", line)
        m2 === nothing || isnan(tt) || push!(out, (tt, parse(Float64, m2.captures[1]) * (parse(Float64, m2.captures[2]) + parse(Float64, m2.captures[3]) + parse(Float64, m2.captures[4]))))
    end
    return out
end
river = isfile(LOG) ? parse_river(LOG) : Tuple{Float64, Float64}[]
river_rate(a, b) = (v = [q for (tt, q) in river if a ≤ tt ≤ b]; isempty(v) ? NaN : mean(v))     # mmol N/s

# sediment: daily means of the pools and fluxes over the shelf columns
bottom_files = ["$(TAG)_bottom_rank$(r).jld2" for r in 0:99 if isfile("$(TAG)_bottom_rank$(r).jld2")]
parts = map(bottom_files) do f; jldopen(f) do file; (I = file["I_OFF"], J = file["J_OFF"], kb = file["k_bottom"], dzb = file["Δz_bottom"]); end; end
iters(f, name) = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/$name"]))))
function exponential_scale(N, H, Δtop)
    top(h) = H * expm1(H / N / h) / expm1(H / h)
    lo, hi = H / 50, 100H
    for _ in 1:200; mid = (lo + hi) / 2; top(mid) > Δtop ? (hi = mid) : (lo = mid); end
    return (lo + hi) / 2
end
zg = RectilinearGrid(size = (1, 1, Nz), x = (0, 1), y = (0, 1), z = ExponentialDiscretization(Nz, -5500, 0; scale = exponential_scale(Nz, 5500, 2.0)))
Δz = diff(collect(znodes(zg, Face())))
Nx = maximum(p.I + size(p.kb, 1) for p in parts); Ny = maximum(p.J + size(p.kb, 2) for p in parts)
depth = zeros(Nx, Ny)
for p in parts, j in axes(p.kb, 2), i in axes(p.kb, 1)
    k = p.kb[i, j]; depth[p.I + i, p.J + j] = k > 0 ? sum(Δz[k+1:end]) + p.dzb[i, j] : 0.0
end
shelf = (depth .> 0) .& (depth .< SHELF)
area = [R_EARTH^2 * cosd(φ_south + (j - 0.5) / CPD) * deg2rad(1 / CPD)^2 for i in 1:Nx, j in 1:Ny]
function read2(name, n)
    A = zeros(Nx, Ny)
    for (p, r) in zip(parts, 0:length(parts)-1)
        jldopen("$(TAG)_bgc_sediment_rank$(r).jld2") do f
            it = iters(f, name)[n]; B = f["timeseries/$name/$it"]
            nx, ny = size(p.kb); hx = (size(B, 1) - nx) ÷ 2; hy = (size(B, 2) - ny) ÷ 2
            A[p.I+1:p.I+nx, p.J+1:p.J+ny] .= B[hx+1:hx+nx, hy+1:hy+ny, 1]
        end
    end
    return A
end
st = jldopen(f -> [f["timeseries/t/$it"] for it in iters(f, "Ns")], "$(TAG)_bgc_sediment_rank0.jld2")
F = Float64[]; Ret = Float64[]; Nsum = Float64[]
for n in eachindex(st)
    Ns, Nf, Nr = read2("Ns", n), read2("Nf", n), read2("Nr", n)
    fl = read2("tracked_sPOM", n) .+ read2("tracked_bPOM", n)
    push!(F, sum(fl[shelf] .* area[shelf]))                                       # mmol N/s
    push!(Ret, sum((0.2 / 86400 .* Ns .+ 2.0 / 86400 .* Nf)[shelf] .* area[shelf]))
    push!(Nsum, sum((Ns .+ Nf .+ Nr)[shelf] .* area[shelf]))
end
# sediment day n holds the mean over the day before st[n]
day(x) = x / 86400
tk = 1e-6
say("\nper day (from the sediment's daily means: day ending at the time; kmol N/day)")
say("  day ending  dW/dt      river R   lateral L  return    sinking F  sponge    ScaleNeg   X (remainder)")
Xs = Float64[]; rows_out = []
for n in eachindex(st)
    t1 = st[n]; t0 = t1 - 86400
    (t0 ≥ t[1] && t1 ≤ t[end]) || continue
    i0 = argmin(abs.(t .- t0)); i1 = argmin(abs.(t .- t1))
    dW = (W[i1] - W[i0]) * tk; Lt = (Lat[i1] - Lat[i0]) * tk; Sp = (Spo[i1] - Spo[i0]) * tk; Ng = (Neg[i1] - Neg[i0]) * tk
    R = river_rate(t0, t1) * 86400 * tk
    Fd = F[n] * 86400 * tk; Rt = Ret[n] * 86400 * tk
    X = dW - R + Lt - Rt + Fd - Sp - Ng
    push!(Xs, X); push!(rows_out, (day(t1), dW, R, Lt, Rt, Fd, Sp, Ng, X))
    say(@sprintf("  %9.1f %+10.3g %+10.3g %+10.3g %+10.3g %+10.3g %+10.3g %+10.3g %+10.3g", day(t1), dW, R, Lt, Rt, Fd, Sp, Ng, X))
end
FROM = parse(Float64, get(ENV, "BGC_BUDGET_FROM_DAY", string(round(day(t[1]) + 4))))
sel = [r for r in rows_out if r[1] ≥ FROM]
say(@sprintf("\nmeans from day %.0f on (pools spun up), %d days, kmol N/day:", FROM, length(sel)))
names = ("dW/dt (change of N in the shelf water)", "river input R", "lateral export L (with tides)", "return from the sediment", "sinking into the sediment F", "sponge N source", "ScaleNegativeTracers N change", "remainder X")
for (c, name) in enumerate(names)
    v = [r[c+1] for r in sel]
    say(@sprintf("  %-42s %+10.3g   (std %.2g)", name, mean(v), std(v)))
end
say(@sprintf("  mean N in the shelf water %.3g kmol; the remainder X is %.2f %% of it per day", mean(W) * tk, 100 * mean([r[9] for r in sel]) / (mean(W) * tk)))
say("\nsign convention: dW/dt = R − L + return − F + sponge + ScaleNeg + X; L outward; X = dW/dt − R + L − return + F − sponge − ScaleNeg")
