# Nitrogen budget of the shelf columns (water depth < BGC_SHELF_DEPTH, default 100 m) from the de-tided daily output of a run of
# bgc/mab_bgc.jl, to trace where the nitrogen goes.
#
#   N in the water  W = Σ dV (NO₃ + NH₄ + P + Z + DOM + sPOM + bPOM)                     (mmol N, shelf volume)
#   N in the sediment S = Σ A (Ns + Nf + Nr)                                                (multi-G runs only)
#   rivers          R = Q (NO₃ + NH₄ + DOM) of the discharge-weighted river chemistry in the run's log
#   export          E = net advective flux of the same N out through the shelf's lateral faces (the 100 m isobath and the open edges),
#                       from the daily u, v and the mean of the two cells' concentrations
#   non-conservative X = d(W + S)/dt − R + E : what is left, i.e. the sponge restoring at the open boundaries, ScaleNegativeTracers'
#                       zeroing of invalid cells, numerical diffusion and the tidal part of the flux that the daily de-tided fields miss
#
# LOBSTER and the sediment conserve nitrogen apart from the refractory pool (a store, counted in S), so a clearly negative X is a
# loss that no process of the BGC model explains. Daily frames, rates by centred differences over two frames. Run:
#   BGC_TAG=/t0/.../bgc_feb19_woa_lim_mg30 BGC_LOG=/t0/.../bgc_feb19_woa_lim_mg30_12243.log julia --project=<env with CairoMakie> bgc/nitrogen_budget.jl
using Printf, Statistics, Dates
using CairoMakie
using Oceananigans
using Oceananigans.Grids: znodes
using NumericalEarth.DataWrangling: jldopen

const TAG   = ENV["BGC_TAG"]
const LOG   = get(ENV, "BGC_LOG", "")
const START = DateTime(get(ENV, "MAB_START_DATE", "2019-02-01"))
const CPD   = parse(Int, get(ENV, "MAB_CELLS_PER_DEGREE", "12"))
const SHELF = parse(Float64, get(ENV, "BGC_SHELF_DEPTH", "100"))
const OUT   = get(ENV, "BGC_OUT", "/t0/workdir/enrique/runs/analysis/nitrogen_" * basename(TAG))
const H, Nz = 7, parse(Int, get(ENV, "MAB_NZ", "50"))
const R_EARTH = 6.371e6
const λ_west, φ_south = -76 + 2 / 12, 34 + 2 / 12
mkpath(OUT)
summary_io = open(joinpath(OUT, "summary.txt"), "w")
function say(args...)
    println(args...); flush(stdout)
    println(summary_io, args...); flush(summary_io)
end

function exponential_scale(N, H, Δtop)
    top(h) = H * expm1(H / N / h) / expm1(H / h)
    lo, hi = H / 50, 100H
    for _ in 1:200; mid = (lo + hi) / 2; top(mid) > Δtop ? (hi = mid) : (lo = mid); end
    return (lo + hi) / 2
end
Δz_top, Z_BOTTOM = parse(Float64, get(ENV, "MAB_DZ_TOP", "2")), parse(Float64, get(ENV, "MAB_ZBOTTOM", "5500"))
zgrid = RectilinearGrid(size = (1, 1, Nz), x = (0, 1), y = (0, 1),
                        z = ExponentialDiscretization(Nz, -Z_BOTTOM, 0; scale = exponential_scale(Nz, Z_BOTTOM, Δz_top)))
const Δz = diff(collect(znodes(zgrid, Face())))          # bottom first

iters(f, name) = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/$name"]))))
bottom_files = ["$(TAG)_bottom_rank$(r).jld2" for r in 0:99 if isfile("$(TAG)_bottom_rank$(r).jld2")]
const NR = length(bottom_files)
parts = map(bottom_files) do f
    jldopen(f) do file; (I = file["I_OFF"], J = file["J_OFF"], kb = file["k_bottom"], dzb = file["Δz_bottom"]); end
end
const Nx = maximum(p.I + size(p.kb, 1) for p in parts); const Ny = maximum(p.J + size(p.kb, 2) for p in parts)
kb = zeros(Int, Nx, Ny); depth = zeros(Nx, Ny); dzb = zeros(Nx, Ny)
for p in parts, j in axes(p.kb, 2), i in axes(p.kb, 1)
    k = p.kb[i, j]; kb[p.I + i, p.J + j] = k; dzb[p.I + i, p.J + j] = p.dzb[i, j]
    depth[p.I + i, p.J + j] = k > 0 ? sum(Δz[k+1:end]) + p.dzb[i, j] : 0.0
end
wet = kb .> 0
shelf = wet .& (depth .< SHELF)
φc(j) = φ_south + (j - 0.5) / CPD
const Δy = R_EARTH * deg2rad(1 / CPD)
Δx(j) = R_EARTH * cosd(φc(j)) * deg2rad(1 / CPD)
area = [Δx(j) * Δy for i in 1:Nx, j in 1:Ny]
say("shelf: $(count(shelf)) columns shallower than $SHELF m of $(count(wet)) wet columns, area $(round(sum(area[shelf]) / 1e9, digits = 1)) thousand km²")
thickness(i, j, k, s) = (k == kb[i, j] ? dzb[i, j] : Δz[k]) * s[i, j]      # level thickness with the z* factor

# ---------------- reading ----------------
function read3(tag, stem, name, n; faces = (0, 0))
    A = zeros(Float32, Nx + faces[1], Ny + faces[2], Nz)
    for (p, r) in zip(parts, 0:NR-1)
        jldopen("$(tag)_$(stem)_rank$(r).jld2") do f
            it = iters(f, name)[n]
            B = f["timeseries/$name/$it"]
            nx, ny = size(p.kb)
            ex = max(0, size(B, 1) - nx - 2H); ey = max(0, size(B, 2) - ny - 2H)
            B = B[H+1:H+nx+ex, H+1:H+ny+ey, H+1:H+Nz]
            nxk = (faces[1] == 1 && ex == 1) ? nx + 1 : nx
            nyk = (faces[2] == 1 && ey == 1) ? ny + 1 : ny
            A[p.I+1:p.I+nxk, p.J+1:p.J+nyk, :] .= B[1:nxk, 1:nyk, :]
        end
    end
    return A
end
function frame_times(tag, stem, name)
    jldopen("$(tag)_$(stem)_rank0.jld2") do f; [f["timeseries/t/$it"] for it in iters(f, name)]; end
end
function read_eta(tag, n)
    A = zeros(Float32, Nx, Ny)
    for (p, r) in zip(parts, 0:NR-1)
        jldopen("$(tag)_eta_daily_rank$(r).jld2") do f
            it = iters(f, "η")[n]; B = f["timeseries/η/$it"]
            nx, ny = size(p.kb); hx = (size(B, 1) - nx) ÷ 2; hy = (size(B, 2) - ny) ÷ 2
            A[p.I+1:p.I+nx, p.J+1:p.J+ny] .= B[hx+1:hx+nx, hy+1:hy+ny, 1]
        end
    end
    return A
end
function read2(tag, name, n)
    A = zeros(Float64, Nx, Ny)
    for (p, r) in zip(parts, 0:NR-1)
        jldopen("$(tag)_bgc_sediment_rank$(r).jld2") do f
            it = iters(f, name)[n]; B = f["timeseries/$name/$it"]
            nx, ny = size(p.kb); hx = (size(B, 1) - nx) ÷ 2; hy = (size(B, 2) - ny) ÷ 2
            A[p.I+1:p.I+nx, p.J+1:p.J+ny] .= B[hx+1:hx+nx, hy+1:hy+ny, 1]
        end
    end
    return A
end

const NAMES = ("NO₃", "NH₄", "P", "Z", "DOM", "sPOM", "bPOM")
vt = round.(Int, frame_times(TAG, "bgc_volume_daily", "P"))
et = round.(Int, frame_times(TAG, "eta_daily", "η"))
ut = round.(Int, frame_times(TAG, "volume_daily", "u"))
vt == et || say("note: the η frames differ from the BGC frames; matching by time")
has_sed = isfile("$(TAG)_bgc_sediment_rank0.jld2")
st = has_sed ? round.(Int, jldopen(f -> [f["timeseries/t/$it"] for it in iters(f, "Ns")], "$(TAG)_bgc_sediment_rank0.jld2")) : Int[]
nearest(ts, t) = (m = argmin(abs.(ts .- t)); abs(ts[m] - t) ≤ 43200 ? m : 0)
Nf = length(vt)
say("frames: $Nf daily, day $(round(first(vt) / 86400, digits = 1)) to $(round(last(vt) / 86400, digits = 1)); sediment fields: $has_sed")

# ---------------- rivers from the log ----------------
function parse_river(log)
    out = Tuple{Float64, Float64}[]      # (time in s, N input in mmol N/s)
    t = NaN
    for line in eachline(log)
        m = match(r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2})\s", line)
        m === nothing || (t = Dates.value(DateTime(m.captures[1], dateformat"yyyy-mm-dd HH:MM") - START) / 1000)
        m2 = match(r"rivers: Q = ([\d.]+) m³/s.*DOM ([\d.]+), NH₄ ([\d.]+), NO₃ ([\d.]+)", line)
        m2 === nothing || isnan(t) || push!(out, (t, parse(Float64, m2.captures[1]) * (parse(Float64, m2.captures[2]) + parse(Float64, m2.captures[3]) + parse(Float64, m2.captures[4]))))
    end
    return out
end
river = isfile(LOG) ? parse_river(LOG) : Tuple{Float64, Float64}[]
say(isfile(LOG) ? "river input: $(length(river)) samples from the log" : "no log given: river input not counted")
river_input(t) = (v = [q for (tt, q) in river if abs(tt - t) ≤ 43200]; isempty(v) ? NaN : mean(v))     # mmol N/s

# ---------------- one pass over the frames ----------------
Vol = fill(NaN, Nf); Ev = fill(NaN, Nf); W = fill(NaN, Nf); S = fill(NaN, Nf); Fsed = fill(NaN, Nf); Rel = fill(NaN, Nf); Nr = fill(NaN, Nf); E = fill(NaN, Nf); Rv = fill(NaN, Nf)
forms = Dict(n => fill(NaN, Nf) for n in NAMES)
for n in 1:Nf
    eta = read_eta(TAG, nearest(et, vt[n])); s = 1 .+ eta ./ max.(depth, 1)
    C = Dict(name => read3(TAG, "bgc_volume_daily", name, n) for name in NAMES)
    Ntot = sum(values(C))
    u = read3(TAG, "volume_daily", "u", nearest(ut, vt[n]); faces = (1, 0))
    v = read3(TAG, "volume_daily", "v", nearest(ut, vt[n]); faces = (0, 1))
    w_tot = 0.0; vol = 0.0; w_form = Dict(name => 0.0 for name in NAMES)
    for j in 1:Ny, i in 1:Nx
        shelf[i, j] || continue
        for k in kb[i, j]:Nz
            dV = area[i, j] * thickness(i, j, k, s)
            vol += dV; w_tot += dV * Ntot[i, j, k]
            for name in NAMES; w_form[name] += dV * C[name][i, j, k]; end
        end
    end
    Vol[n] = vol; W[n] = w_tot; for name in NAMES; forms[name][n] = w_form[name]; end
    # lateral export (mmol N/s): outward positive
    out = 0.0; outv = 0.0
    for j in 1:Ny, i in 1:Nx
        shelf[i, j] || continue
        for k in kb[i, j]:Nz
            dz = thickness(i, j, k, s)
            c0 = Ntot[i, j, k]
            # east, west, north, south neighbours: (neighbour inside?, face velocity with outward sign, neighbour concentration)
            if i == Nx || !shelf[i+1, j]
                nb = i < Nx && wet[i+1, j] && k ≥ kb[i+1, j] ? Ntot[i+1, j, k] : c0
                (i == Nx || (wet[i+1, j] && k ≥ kb[i+1, j])) && (out += u[i+1, j, k] * Δy * dz * 0.5 * (c0 + nb); outv += u[i+1, j, k] * Δy * dz)
            end
            if i == 1 || !shelf[i-1, j]
                nb = i > 1 && wet[i-1, j] && k ≥ kb[i-1, j] ? Ntot[i-1, j, k] : c0
                (i == 1 || (wet[i-1, j] && k ≥ kb[i-1, j])) && (out -= u[i, j, k] * Δy * dz * 0.5 * (c0 + nb); outv -= u[i, j, k] * Δy * dz)
            end
            if j == Ny || !shelf[i, j+1]
                nb = j < Ny && wet[i, j+1] && k ≥ kb[i, j+1] ? Ntot[i, j+1, k] : c0
                (j == Ny || (wet[i, j+1] && k ≥ kb[i, j+1])) && (out += v[i, j+1, k] * Δx(j) * dz * 0.5 * (c0 + nb); outv += v[i, j+1, k] * Δx(j) * dz)
            end
            if j == 1 || !shelf[i, j-1]
                nb = j > 1 && wet[i, j-1] && k ≥ kb[i, j-1] ? Ntot[i, j-1, k] : c0
                (j == 1 || (wet[i, j-1] && k ≥ kb[i, j-1])) && (out -= v[i, j, k] * Δx(j) * dz * 0.5 * (c0 + nb); outv -= v[i, j, k] * Δx(j) * dz)
            end
        end
    end
    E[n] = out * 86400                                      # mmol N/day
    Ev[n] = outv * 86400                                    # m³/day
    Rv[n] = river_input(vt[n]) * 86400
    if has_sed
        m = nearest(st, vt[n])
        if m > 0
            Ns, Nfast, Nref = read2(TAG, "Ns", m), read2(TAG, "Nf", m), read2(TAG, "Nr", m)
            fl = read2(TAG, "tracked_sPOM", m) .+ read2(TAG, "tracked_bPOM", m)
            S[n] = sum((Ns .+ Nfast .+ Nref)[shelf] .* area[shelf]); Nr[n] = sum(Nref[shelf] .* area[shelf])
            Fsed[n] = sum(fl[shelf] .* area[shelf]) * 86400   # mmol N/day
            # N the sediment gives back to the water: λs Ns + λf Nf (mmol N/m²/s; 2/day and 0.2/day, SimpleMultiG's defaults)
            Rel[n] = sum((0.2 / 86400 .* Ns .+ 2.0 / 86400 .* Nfast)[shelf] .* area[shelf]) * 86400
        end
    end
    (n % 5 == 0) && say("  frame $n of $Nf")
end

# ---------------- budget ----------------
# values below in kmol N (1e6 mmol) per day unless stated
tk = 1e-6        # mmol → kmol
d(x) = [i == 1 || i == length(x) ? NaN : (x[i+1] - x[i-1]) / ((vt[i+1] - vt[i-1]) / 86400) for i in eachindex(x)]
dW = d(W); dS = has_sed ? d(S) : zeros(Nf)
X = dW .+ dS .- Rv .+ E
say("\nnitrogen budget of the shelf columns, kmol N per day (mean over the frames with every term; positive = gain of the shelf water + sediment)")
ok = [i for i in 1:Nf if all(isfinite, (dW[i], dS[i], Rv[i], E[i]))]
mn(x) = isempty(ok) ? NaN : mean(x[ok]) * tk
say(@sprintf("  frames used: %d (day %.0f to %.0f)", length(ok), vt[first(ok)] / 86400, vt[last(ok)] / 86400))
say(@sprintf("  N in the water (W):             mean %10.3g kmol N, change dW/dt %+10.3g", mean(W[ok]) * tk, mn(dW)))
say(@sprintf("  N in the sediment (S):          mean %10.3g kmol N, change dS/dt %+10.3g   (refractory pool: %+10.3g)", has_sed ? mean(S[ok]) * tk : 0.0, mn(dS), has_sed ? mean(d(Nr)[ok]) * tk : 0.0))
say(@sprintf("  river input (R):                         %+10.3g", mn(Rv)))
say(@sprintf("  advective export through the shelf's faces (E):   %+10.3g", mn(E)))
say(@sprintf("  sinking into the sediment (flux of sPOM + bPOM): %+10.3g", has_sed ? mn(Fsed) : NaN))
say(@sprintf("  non-conservative X = dW/dt + dS/dt − R + E:      %+10.3g   (%.2f %% of the N in the water per day)", mn(X), 100 * mn(X) / (mean(W[ok]) * tk)))
if has_sed
    Xw = dW .- Rv .+ E .+ Fsed .- Rel
    say("\nwater only, with the sediment's own source and sink: dW/dt = R − E − (sinking into the sediment) + (return from the sediment) + Xw")
    say(@sprintf("  sinking into the sediment F: %+10.3g   return from the sediment (λs Ns + λf Nf): %+10.3g   net to the sediment F − return: %+10.3g (%.1f %% of F; 10 %% is the refractory fraction)", mn(Fsed), mn(Rel), mn(Fsed .- Rel), 100 * mean((Fsed .- Rel)[ok]) / mean(Fsed[ok])))
    say(@sprintf("  Xw = dW/dt − R + E + F − return: %+10.3g   (%.2f %% of the N in the water per day; standard error of the mean over the frames %.2g)", mn(Xw), 100 * mn(Xw) / (mean(W[ok]) * tk), std(Xw[ok]) * tk / sqrt(length(ok))))
    say(@sprintf("  refractory pool growth %+10.3g kmol N/day against 0.1 F = %+10.3g: ratio %.2f", mean(d(Nr)[ok]) * tk, 0.1 * mn(Fsed), mean(d(Nr)[ok]) / (0.1 * mean(Fsed[ok]))))
end
dV = d(Vol)
say(@sprintf("\nvolume closure of the shelf columns (flux calculation check): dV/dt = %+.3g m³/day, net volume flux out of the shelf faces %+.3g, their sum %+.3g (river water, precipitation and evaporation are about 1.6e8 m³/day); the mean of the volume in the shelf columns is %.3g m³ and the volume flux out of the shelf faces through a day is of order %.2g m³/day",
             mean(dV[ok]), mean(Ev[ok]), mean((dV .+ Ev)[ok]), mean(Vol[ok]), mean(abs.(Ev[ok]))))
say("\nN in the water by form, mean over the frames (kmol N): " * join([@sprintf("%s %.3g", n, mean(forms[n][ok]) * tk) for n in NAMES], ", "))
say("change of each form from the first to the last frame (kmol N): " * join([@sprintf("%s %+.3g", n, (forms[n][end] - forms[n][1]) * tk) for n in NAMES], ", "))
say("\nper frame (kmol N/day): day, W (kmol N), dW/dt, dS/dt, R, E, X")
for i in ok
    say(@sprintf("  %5.1f %10.3g %+10.3g %+10.3g %+10.3g %+10.3g %+10.3g", vt[i] / 86400, W[i] * tk, dW[i] * tk, dS[i] * tk, Rv[i] * tk, E[i] * tk, X[i] * tk))
end

days = vt ./ 86400
fig = Figure(size = (1100, 600))
ax = Axis(fig[1, 1], title = "Shelf nitrogen budget (kmol N per day)", xlabel = "run day", ylabel = "kmol N / day")
for (label, y, color) in (("dW/dt", dW, :black), ("dS/dt", dS, :brown), ("river input R", Rv, :royalblue), ("export E", -E, :seagreen), ("non-conservative X", X, :firebrick))
    lines!(ax, days, y .* tk; label = label, color = color, linewidth = label == "non-conservative X" ? 3 : 1.5)
end
axislegend(ax, position = :rb)
ax2 = Axis(fig[1, 2], title = "N in the shelf water by form (kmol N)", xlabel = "run day")
for (name, color) in zip(NAMES, (:royalblue, :seagreen, :goldenrod, :orange, :purple, :grey, :brown))
    lines!(ax2, days, forms[name] .* tk; label = name, color = color)
end
axislegend(ax2, position = :rt)
save(joinpath(OUT, "nitrogen_budget.png"), fig)
say("\nfigure and summary in $OUT")
