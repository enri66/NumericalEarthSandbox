# ERA5's 10 m wind against the OOI Pioneer Central Surface Mooring (CP01CNSM, 40.13°N 70.78°W) in Aug-Oct 2019: is
# the forcing short of near-inertial energy? The buoy carries two bulk meteorology packages (METBK A and B, one-minute
# winds, earth-relative); their hourly means are checked against each other, averaged, and taken to 10 m with a
# neutral log profile from the anemometer height Z_ANEMOMETER (z₀ = 10⁻⁴ m). ERA5 is interpolated bilinearly to the
# buoy from the MAB-box hourly files. Both winds go through the same stress formula, τ = ρₐ C_D |U| U with the
# Large & Pond (1981) neutral drag coefficient, and are compared as
#   overall statistics of speed and stress;
#   clockwise and counterclockwise stress variance in frequency bands (near-inertial 16-22 h, f at 40.1°N: 18.6 h),
#       from the rotary spectrum of the hourly stress;
#   the near-inertial mixed-layer currents each forcing drives in a slab model (Pollard & Millard 1970),
#       dZ/dt + (r + i f) Z = T / (ρ₀ H),  Z = u + i v,  T = τˣ + i τʸ,  with H = SLAB_H and r = 1 / (4 days),
#       band-passed to 16-22 h.
# Usage:
#   OOI_DIR=/t0/workdir/enrique/Data/OOI/pioneer julia --project=. scripts/era5_vs_cnsm_wind.jl
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))
const FFTW = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "FFTW")]

const OOI_DIR = get(ENV, "OOI_DIR", joinpath(homedir(), "Data", "OOI", "pioneer"))
const ERA5_DIR = joinpath(GDIR, "era5")
const BOX = "_-76.0_-64.0_34.0_42.0.nc"
const Z_ANEMOMETER = parse(Float64, get(ENV, "Z_ANEMOMETER", "4.0"))   # m above the sea surface
const SLAB_H = parse(Float64, get(ENV, "SLAB_H", "40"))                # m
const T0, T1 = DateTime(2019, 8, 29), DateTime(2019, 11, 1)
const LON, LAT = -70.7783, 40.1333
const ρₐ, ρ₀ = 1.22, 1025.0
const f = 2 * 7.292e-5 * sind(LAT)
const STORMS = (("Dorian, 6-12 Sep", DateTime(2019, 9, 6), DateTime(2019, 9, 13)),
                ("17-23 Oct", DateTime(2019, 10, 17), DateTime(2019, 10, 24)))

# ---------------- buoy ----------------
function read_metbk(dataset)
    lines = readlines(joinpath(OOI_DIR, dataset * ".csv"))
    col = Dict(h => i for (i, h) in enumerate(split(lines[1], ",")))
    sums = Dict{DateTime, NTuple{3, Float64}}()
    for line in lines[3:end]
        v = split(line, ",")
        (v[col["eastward_wind_qc_agg"]] == "4" || v[col["northward_wind_qc_agg"]] == "4") && continue
        u = tryparse(Float64, v[col["eastward_wind"]]); w = tryparse(Float64, v[col["northward_wind"]])
        (isnothing(u) || isnothing(w) || !isfinite(u) || !isfinite(w) || hypot(u, w) > 60) && continue
        t = round(DateTime(v[col["time"]][1:19]), Hour)          # hourly means centred on the hour, like ERA5's instant
        s = get(sums, t, (0.0, 0.0, 0.0)); sums[t] = (s[1] + u, s[2] + w, s[3] + 1)
    end
    return Dict(t => (s[1] / s[3], s[2] / s[3]) for (t, s) in sums if s[3] >= 30)
end
A = read_metbk("ooi-cp01cnsm-sbd11-06-metbka000"); B = read_metbk("ooi-cp01cnsm-sbd12-06-metbka000")
hours = collect(T0:Hour(1):T1)
log_factor = log(10 / 1e-4) / log(Z_ANEMOMETER / 1e-4)
function buoy(t)
    a = get(A, t, nothing); b = get(B, t, nothing)
    isnothing(a) && isnothing(b) && return (NaN, NaN)
    isnothing(a) && return b .* log_factor
    isnothing(b) && return a .* log_factor
    return ((a[1] + b[1]) / 2, (a[2] + b[2]) / 2) .* log_factor
end
both = [t for t in hours if haskey(A, t) && haskey(B, t)]
sA = [hypot(A[t]...) for t in both]; sB = [hypot(B[t]...) for t in both]
@printf("METBK A vs B (%d common hours): mean speed %.2f vs %.2f m/s, corr %.3f, rms difference %.2f m/s\n",
        length(both), mean(sA), mean(sB), cor(sA, sB), sqrt(mean((sA .- sB) .^ 2)))
@printf("anemometer height %.1f m → 10 m factor %.3f (neutral log profile)\n", Z_ANEMOMETER, log_factor)

# ---------------- ERA5 ----------------
NCD5 = NCD
ds = NCD5.Dataset(joinpath(ERA5_DIR, "10m_u_component_of_wind_ERA5HourlySingleLevel_2019-09-01T00$BOX"))
elon = Float64.(ds["longitude"][:]); elat = Float64.(ds["latitude"][:]); close(ds)
order = sortperm(elat); elat_sorted = elat[order]
function era5(name, var, t)
    file = joinpath(ERA5_DIR, "$(name)_ERA5HourlySingleLevel_$(Dates.format(t, "yyyy-mm-ddTHH"))$BOX")
    isfile(file) || return NaN
    ds = NCD5.Dataset(file); a = Float64.(ds[var][:, :, 1])[:, order]; close(ds)
    return bilinear(a, elon, elat_sorted, LON, LAT)
end
E = [(era5("10m_u_component_of_wind", "u10", t), era5("10m_v_component_of_wind", "v10", t)) for t in hours]
O = [buoy(t) for t in hours]

# ---------------- stress and statistics ----------------
CD(U) = U < 11 ? 1.2e-3 : (0.49 + 0.065U) * 1e-3          # Large & Pond (1981), neutral, U at 10 m
stress((u, v)) = (U = hypot(u, v); ρₐ * CD(U) * U .* (u, v))
τE = stress.(E); τO = stress.(O)
ok = [all(isfinite, O[n]) && all(isfinite, E[n]) for n in eachindex(hours)]
speed(x) = hypot(x...)
@printf("\n%d hours with both ERA5 and the buoy (of %d)\n", count(ok), length(hours))
@printf("10 m speed: buoy %.2f, ERA5 %.2f m/s (ratio %.2f), corr %.3f; std buoy %.2f ERA5 %.2f\n", mean(speed.(O[ok])), mean(speed.(E[ok])),
        mean(speed.(E[ok])) / mean(speed.(O[ok])), cor(speed.(O[ok]), speed.(E[ok])), std(speed.(O[ok])), std(speed.(E[ok])))
@printf("|τ|: buoy %.4f, ERA5 %.4f N/m² (ratio %.2f)\n", mean(speed.(τO[ok])), mean(speed.(τE[ok])), mean(speed.(τE[ok])) / mean(speed.(τO[ok])))
for (label, a, b) in STORMS
    w = [ok[n] && a <= hours[n] < b for n in eachindex(hours)]
    @printf("  %-17s max speed buoy %.1f, ERA5 %.1f m/s; mean |τ| buoy %.3f, ERA5 %.3f N/m²\n", label, maximum(speed.(O[w])),
            maximum(speed.(E[w])), mean(speed.(τO[w])), mean(speed.(τE[w])))
end

# Gaps (only where the buoy is missing) filled linearly in time before the spectral steps
function filled(z)
    x = copy(z); good = findall(isfinite, x)
    for n in eachindex(x)
        isfinite(x[n]) && continue
        a = findlast(<(n), good); b = findfirst(>(n), good)
        x[n] = isnothing(a) ? x[good[b]] : isnothing(b) ? x[good[a]] : x[good[a]] + (x[good[b]] - x[good[a]]) * (n - good[a]) / (good[b] - good[a])
    end
    return x
end
complex_series(τ) = filled([t[1] for t in τ]) .+ im .* filled([t[2] for t in τ])
ZE = complex_series(τE); ZO = complex_series(τO)

# Rotary band variance: with Z = τˣ + i τʸ and exp(+iωt) counterclockwise, clockwise motion has negative frequency
freq = FFTW.fftfreq(length(ZE), 1.0)                         # cycles per hour
function band_variance(Z, lo, hi; clockwise)
    Y = FFTW.fft(Z .- mean(Z)); P = abs2.(Y) ./ length(Z)^2
    sel = [(clockwise ? -fr : fr) > 0 && 1 / hi <= abs(fr) <= 1 / lo for fr in freq]
    return sum(P[sel])
end
println("\nStress variance by band and rotation, ERA5 / buoy:")
for (label, lo, hi) in (("subinertial >30 h", 30.0, 1e5), ("near-inertial 16-22 h", 16.0, 22.0), ("semidiurnal-ish 10-16 h", 10.0, 16.0), ("high <10 h", 2.0, 10.0))
    for cw in (true, false)
        e = band_variance(ZE, lo, hi; clockwise = cw); o = band_variance(ZO, lo, hi; clockwise = cw)
        @printf("  %-24s %-17s buoy %.2e  ERA5 %.2e (N/m²)²  ratio %.2f\n", label, cw ? "clockwise" : "counterclockwise", o, e, e / o)
    end
end

# ---------------- slab model ----------------
function slab(T; H = SLAB_H, r = 1 / (4 * 86400), dt = 3600.0)
    Z = zeros(ComplexF64, length(T))
    for n in 2:length(T)                                     # exact for piecewise-constant forcing over each hour
        a = r + im * f
        F = (T[n-1] + T[n]) / 2 / (ρ₀ * H)
        Z[n] = Z[n-1] * exp(-a * dt) + F * (1 - exp(-a * dt)) / a
    end
    return Z
end
function bandpass(Z, lo, hi)
    Y = FFTW.fft(Z .- mean(Z)); Y[.!(1 / hi .<= abs.(freq) .<= 1 / lo)] .= 0
    return FFTW.ifft(Y)
end
uE = bandpass(slab(ZE), 16.0, 22.0); uO = bandpass(slab(ZO), 16.0, 22.0)
kE = 0.5 .* abs2.(uE); kO = 0.5 .* abs2.(uO)
@printf("\nSlab model (H = %g m, r = 1/4 days): near-inertial KE, ERA5 forcing / buoy forcing = %.2e / %.2e J/kg (ratio %.2f)\n",
        SLAB_H, mean(kE), mean(kO), mean(kE) / mean(kO))
for (label, a, b) in STORMS
    w = [a <= t < b for t in hours]
    @printf("  %-17s %.2e / %.2e (ratio %.2f)\n", label, mean(kE[w]), mean(kO[w]), mean(kE[w]) / mean(kO[w]))
end

# ---------------- figure ----------------
td = [Dates.value(t - T0) / 86_400_000 for t in hours]
fig = Figure(size = (1500, 1000), fontsize = 14)
Label(fig[0, 1:2], "ERA5 vs OOI Pioneer Central Surface Mooring (40.13°N, 70.78°W), Aug-Oct 2019", fontsize = 18)
ax = Axis(fig[1, 1:2], ylabel = "10 m wind speed (m/s)", xlabel = "days since 2019-08-29")
lines!(ax, td, speed.(O); color = :black, label = "buoy (METBK A/B, to 10 m)"); lines!(ax, td, speed.(E); color = :firebrick, label = "ERA5")
axislegend(ax; position = :lt)
ax = Axis(fig[2, 1:2], ylabel = "slab near-inertial KE (J/kg)", xlabel = "days since 2019-08-29", title = "slab model, H = $(Int(SLAB_H)) m")
lines!(ax, td, kO; color = :black, label = "buoy forcing"); lines!(ax, td, kE; color = :firebrick, label = "ERA5 forcing")
axislegend(ax; position = :lt)
ax = Axis(fig[3, 1], xscale = log10, yscale = log10, xlabel = "frequency (cycles/day)", ylabel = "clockwise stress spectrum",
          title = "rotary spectra of stress (clockwise solid, counterclockwise dashed)")
for (Z, color) in ((ZO, :black), (ZE, :firebrick))
    Y = FFTW.fft(Z .- mean(Z)); P = abs2.(Y) ./ length(Z)
    pos = findall(>(0), freq); neg = findall(<(0), freq)
    smooth(x) = [mean(x[max(1, n - 3):min(end, n + 3)]) for n in eachindex(x)]
    lines!(ax, 24 .* freq[pos], smooth(P[pos]); color, linestyle = :dash)
    o = sortperm(-freq[neg]); lines!(ax, -24 .* freq[neg][o], smooth(P[neg][o]); color)
end
vlines!(ax, 24 * f / 2π * 3600; color = :gray, linestyle = :dot)
ax = Axis(fig[3, 2], xlabel = "buoy speed (m/s)", ylabel = "ERA5 speed (m/s)", title = "hourly 10 m speed")
scatter!(ax, speed.(O[ok]), speed.(E[ok]); markersize = 3, color = (:black, 0.4)); lines!(ax, [0, 25], [0, 25]; color = :firebrick)
out = joinpath(OOI_DIR, "era5_vs_cnsm_wind.png")
save(out, fig; px_per_unit = 1.2)
println("saved ", out)
