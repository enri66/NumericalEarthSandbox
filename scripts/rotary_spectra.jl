# Rotary frequency spectra of the currents at the OOI Pioneer moorings, ADCP against the runs: is the model's
# near-inertial energy missing, or at a shifted frequency (background vorticity moves the effective inertial frequency
# to f + ζ/2) outside the fixed 16-22 h band of ni_budget.jl? For each mooring and source the layer-mean velocity
# w = u + iv (hourly, gaps filled linearly) is split into 16-day segments overlapping by half, each Hann-windowed; the
# averaged periodograms give clockwise (negative frequency) and counterclockwise spectra. Reported:
#   the frequency of the clockwise peak between 0.85 f and 1.35 f, in units of f (the diurnal tides at 0.72-0.78 f
#   and the semidiurnal at about 1.5 f are kept out);
#   clockwise energy in a narrow band (0.85-1.15 f) and a broader one clear of the tides (0.85-1.38 f), model / ADCP,
#   and the fraction of the broad band in the narrow one (a spread or shifted peak lowers it);
#   for each run with daily output, the median ζ/f at the mooring cell (centred differences of the de-tided top-level
#   currents over the window) and the f + ζ/2 it implies.
# Usage:
#   MAB_TAGS=/t0/.../catke_fall/moor,/t0/.../res_test/r12,/t0/.../res_test/r24 OOI_DIR=/t0/workdir/enrique/Data/OOI/pioneer \
#       julia --project=. scripts/rotary_spectra.jl
using CairoMakie
include(joinpath(@__DIR__, "ooi_common.jl"))

const TAGS   = filter(!isempty, split(get(ENV, "MAB_TAGS", ""), ","))
const AXIS   = collect(DateTime(get(ENV, "SPEC_START", "2019-08-30")):Hour(1):DateTime(get(ENV, "SPEC_END", "2019-10-27")))
const SEG    = parse(Int, get(ENV, "SPEC_SEGMENT_HOURS", "384"))   # hours per segment (16 days: 0.05 f resolution)
const OUTDIR = get(ENV, "SPEC_OUTDIR", OOI_DIR)
const LAYERS = Dict("OSSM" => ((40, 100),), "PMUO" => ((40, 100),), "PMCO" => ((12, 40), (40, 80)), "CNSM" => ((12, 40), (40, 80)))
const f_cpd  = 2 * 7.292e-5 * sind(40.0) * 86400 / 2π         # cycles per day at 40°N (≈ 1.29)

function layer_series(src, (a, b))
    ti = Dict(t => n for (n, t) in enumerate(src.times))
    grid = collect(a:2.0:b)
    U = regrid(src.depths, src.U, grid); V = regrid(src.depths, src.V, grid)
    u = [haskey(ti, t) ? nanmean(U[:, ti[t]]) : NaN for t in AXIS]
    v = [haskey(ti, t) ? nanmean(V[:, ti[t]]) : NaN for t in AXIS]
    return fill_gaps(u) .+ im .* fill_gaps(v)
end

# Welch rotary spectrum; returns frequencies (cpd, ascending, negative = clockwise) and spectral density (m²/s²/cpd)
function rotary_spectrum(w)
    hann = [0.5 - 0.5cos(2π * (n - 1) / SEG) for n in 1:SEG]
    starts = 1:SEG÷2:length(w) - SEG + 1
    P = zeros(SEG)
    for s in starts
        x = w[s:s+SEG-1]; x = (x .- mean(x)) .* hann
        P .+= abs2.(FFTW.fft(x))
    end
    fr = FFTW.fftfreq(SEG, 24.0)                               # cycles per day (hourly samples)
    P ./= length(starts) * sum(abs2, hann)                     # density per cycle per hour (hourly samples)
    o = sortperm(fr)
    return fr[o], P[o] ./ 24.0                                 # density per cpd
end

band(fr, P, lo, hi) = sum(P[(fr .<= -lo * f_cpd) .& (fr .>= -hi * f_cpd)]) * (24.0 / SEG)

# The run's start date (its clock's zero), as stored in its mooring output
function run_start(prefix)
    f = first(filter(f -> occursin(r"_moorings_rank\d+\.jld2$", f) && startswith(basename(f), basename(prefix) * "_moorings"),
                     readdir(dirname(prefix); join = true)))
    return JLD2.jldopen(file -> DateTime(file["start_date"]), f)
end

# Median ζ/f at a cell over the window, from a run's de-tided daily top-level currents
function vorticity_ratio(prefix, cell)
    vol = prefix * "_volume_daily.jld2"
    (isfile(vol) || any(f -> startswith(basename(f), basename(prefix) * "_volume_daily_rank"), readdir(dirname(prefix)))) || return NaN
    U = open_series(vol, "u"; backend = OnDisk()); V = open_series(vol, "v"; backend = OnDisk())
    ug = U.grid.underlying_grid
    λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center()))
    i = argmin(abs.(λ .- cell[1])); j = argmin(abs.(φ .- cell[2]))
    R = 6.371e6; dx = R * cosd(φ[j]) * deg2rad(λ[2] - λ[1]); dy = R * deg2rad(φ[2] - φ[1])
    f = 2 * 7.292e-5 * sind(φ[j])
    t0 = run_start(prefix)
    ratios = Float64[]
    for n in eachindex(U.times)
        t = t0 + Second(round(Int, U.times[n]))
        AXIS[1] <= t <= AXIS[end] || continue
        u = interior(U[n]); v = interior(V[n]); k = size(u, 3)
        ζ = ((v[i+1, j, k] + v[i+1, j+1, k]) - (v[i-1, j, k] + v[i-1, j+1, k])) / (4dx) -
            ((u[i, j+1, k] + u[i+1, j+1, k]) - (u[i, j-1, k] + u[i+1, j-1, k])) / (4dy)
        push!(ratios, ζ / f)
    end
    return isempty(ratios) ? NaN : median(ratios)
end

positions = Dict("OSSM" => (-70.8869, 39.9375), "PMUO" => (-70.7702, 39.9393), "PMCO" => (-70.8792, 40.0968), "CNSM" => (-70.7783, 40.1333))
spectra = Dict{Tuple{String, Tuple{Int, Int}}, Any}()
@printf("f at 40°N = %.3f cpd (%.1f h); window %s to %s, %d-hour segments\n", f_cpd, 24 / f_cpd, AXIS[1], AXIS[end], SEG)
for name in ("OSSM", "PMUO", "PMCO", "CNSM")
    dataset, _, _ = MOORINGS[name]
    sources = [("ADCP", read_adcp(dataset)); [(basename(String(t)), read_model(String(t), name)) for t in TAGS]]
    sources = filter(s -> !isnothing(s[2]), sources)
    println("\n== $name")
    for tag in TAGS
        r = vorticity_ratio(String(tag), positions[name])
        isfinite(r) && @printf("   %-6s median ζ/f at the cell %+.3f → f + ζ/2 = %.3f f\n", basename(String(tag)), r, 1 + r / 2)
    end
    for layer in LAYERS[name]
        res = Dict(label => rotary_spectrum(layer_series(src, layer)) for (label, src) in sources)
        spectra[(name, layer)] = res
        fa, Pa = res["ADCP"]
        nA, wA = band(fa, Pa, 0.85, 1.15), band(fa, Pa, 0.85, 1.38)
        @printf("   %d-%d m   %-6s %12s %16s %16s %14s\n", layer..., "source", "CW peak (f)", "narrow / ADCP", "broad / ADCP", "narrow share")
        for (label, _) in sources
            fr, P = res[label]
            sel = findall(x -> -1.35f_cpd <= x <= -0.85f_cpd, fr)
            peak = -fr[sel[argmax(P[sel])]] / f_cpd
            n, w = band(fr, P, 0.85, 1.15), band(fr, P, 0.85, 1.38)
            @printf("             %-6s %12.2f %16.2f %16.2f %14.2f\n", label, peak, n / nA, w / wA, n / w)
        end
    end
end

# ---------------- figure ----------------
keys_sorted = [(n, l) for n in ("OSSM", "PMUO", "PMCO", "CNSM") for l in LAYERS[n]]
fig = Figure(size = (1500, 330 * cld(length(keys_sorted), 3) + 60), fontsize = 13)
Label(fig[0, 1:3], "rotary spectra of layer-mean currents (clockwise solid, counterclockwise dashed); dotted: f, K1, M2", fontsize = 16)
colors = Dict("ADCP" => :black)
palette = (:firebrick, :royalblue, :darkorange, :seagreen)
for (n, (name, layer)) in enumerate(keys_sorted)
    ax = Axis(fig[cld(n, 3), mod1(n, 3)], title = "$name $(layer[1])-$(layer[2]) m", xscale = log10, yscale = log10,
              xlabel = "frequency (cpd)", ylabel = mod1(n, 3) == 1 ? "PSD (m²/s²/cpd)" : "", limits = (0.2, 12, nothing, nothing))
    res = spectra[(name, layer)]
    for (c, label) in enumerate(sort(collect(keys(res)); by = l -> l == "ADCP" ? "" : l))
        fr, P = res[label]
        color = label == "ADCP" ? :black : palette[mod1(c - 1, length(palette))]
        cw = findall(<(0), fr); ccw = findall(>(0), fr)
        lines!(ax, -fr[cw], max.(P[cw], 1e-8); color, label)
        lines!(ax, fr[ccw], max.(P[ccw], 1e-8); color, linestyle = :dash)
    end
    vlines!(ax, [f_cpd, 1.0027, 1.9323]; color = :gray, linestyle = :dot)
    n == 1 && axislegend(ax; position = :lb)
end
out = joinpath(OUTDIR, get(ENV, "SPEC_FIGURE", "rotary_spectra.png"))
save(out, fig; px_per_unit = 1.2)
println("\nsaved ", out)
