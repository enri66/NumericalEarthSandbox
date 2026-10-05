# The tide of a script-05 run against observations: harmonic constants of sea level at the NOAA CO-OPS gauges
# (noaa_harcon_mab.jl) and of TPXO10 at the same points, and tidal currents at the OOI Pioneer moorings against the
# ADCPs. Least-squares fits use the same (ω, f, phase) the run's tidal forcing uses (earth_tidal_harmonics at the
# run's start date), so a complex amplitude C means η = Re(C f e^{i(ωt + φ)}) and a gauge's (A, G) is C = A e^{-iG}.
# Fitted: M2 S2 N2 K1 O1 (the record is about two months, too short to separate K2 from S2 or P1 from K1); the
# first TIDE_SKIP_DAYS are left out (the tidal ramp).
#   Sea level: model η (hourly output, nearest wet cell to each gauge) against the gauge and TPXO; amplitude ratio and
#   phase difference per constituent, and the complex error |C_model − C_obs|.
#   Currents: the depth-mean current over the common ADCP/model depth range at each mooring, fitted for each
#   constituent in u and v; tidal-ellipse semi-major axis and the model / ADCP ratio.
# Usage:
#   MAB_TAG=/t0/.../res_test/r12 MAB_START_DATE=2019-08-29 OOI_DIR=/t0/workdir/enrique/Data/OOI/pioneer \
#   TPXO_DIR=/t0/workdir/enrique/Data/TPXO10_atlas_v2_nc julia --project=. scripts/tide_check.jl
include(joinpath(@__DIR__, "ooi_common.jl"))
include(joinpath(@__DIR__, "noaa_harcon_mab.jl"))
using NumericalEarth: TPXO10Atlas, earth_tidal_harmonics
using Oceananigans: tidal_atlas_constants

const TAG       = get(ENV, "MAB_TAG", "r12")
const PREFIX    = run_prefix(TAG)
const TPXO_DIR  = get(ENV, "TPXO_DIR", joinpath(homedir(), "Data", "TPXO10_atlas_v2_nc"))
const SKIP_DAYS = parse(Float64, get(ENV, "TIDE_SKIP_DAYS", "2"))
const FIT       = (:M2, :S2, :N2, :K1, :O1)

harmonics = earth_tidal_harmonics(start_date; constituents = FIT, ramp_time = 0)

# Least squares for the complex amplitudes of FIT in a real series x(t), t in seconds since the run's start
function fit(t, x)
    ok = isfinite.(x) .& (t .>= SKIP_DAYS * 86400)
    t, x = t[ok], x[ok]
    cols = [ones(length(t))]
    for n in eachindex(FIT)
        θ = harmonics.frequencies[n] .* t .+ harmonics.phases[n]
        push!(cols, harmonics.nodal_factors[n] .* cos.(θ)); push!(cols, -harmonics.nodal_factors[n] .* sin.(θ))
    end
    β = reduce(hcat, cols) \ x
    return Dict(FIT[n] => complex(β[2n], β[2n+1]) for n in eachindex(FIT))   # Re(C e^{iθ}) = a cos θ − b sin θ
end
amp(C) = abs(C); pha(C) = mod(rad2deg(-angle(C)), 360)           # Greenwich lag of C = A e^{-iG}
dphase(a, b) = mod(pha(a) - pha(b) + 180, 360) - 180

# ---------------- sea level ----------------
η, ηt = surface_series(PREFIX * "_eta.jld2", "η")
vol = PREFIX * "_volume_daily.jld2"
grid = open_series(vol, "T"; backend = OnDisk()).grid; ug = grid.underlying_grid
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center()))
B = grid.immersed_boundary.bottom_height; bh = [B[i, j, 1] for i in 1:size(ug, 1), j in 1:size(ug, 2)]
wet = [(i, j) for i in axes(bh, 1), j in axes(bh, 2) if bh[i, j] < 0]
nearest_wet(lon, lat) = wet[argmin([hypot((λ[i] - lon) * cosd(lat), φ[j] - lat) for (i, j) in wet])]

stations = sort([(id, s) for (id, s) in NOAA_HARCON if λ[1] <= s.lon <= λ[end] && φ[1] <= s.lat <= φ[end]]; by = x -> x[2].lat)
nodes = [(s.lon, s.lat) for (_, s) in stations]
tpxo = Dict(c => getproperty(tidal_atlas_constants(TPXO10Atlas(), nodes, c; dir = TPXO_DIR), :sea_surface_height) for c in FIT)

println("Sea level, $(basename(PREFIX)): model / gauge amplitude and model − gauge phase (°); TPXO / gauge in brackets")
@printf("%-24s %-8s %8s", "station", "kind", "depth")
for c in FIT; @printf("  %20s", c); end
@printf("  %10s %10s\n", "|ΔC| M2", "TPXO |ΔC|")
for (k, (id, s)) in enumerate(stations)
    i, j = nearest_wet(s.lon, s.lat)
    C = fit(ηt, η[i, j, :])
    @printf("%-24s %-8s %6.0f m", first(s.name, 24), s.kind, -bh[i, j])
    for c in FIT
        Cg = s.con[c][1] * cis(-deg2rad(s.con[c][2])); Ct = tpxo[c][k]
        @printf("  %4.2f %+5.0f° (%4.2f)", amp(C[c]) / amp(Cg), dphase(C[c], Cg), amp(Ct) / amp(Cg))
    end
    Cg = s.con[:M2][1] * cis(-deg2rad(s.con[:M2][2]))
    @printf("  %8.3f m %8.3f m\n", abs(C[:M2] - Cg), abs(tpxo[:M2][k] - Cg))
end

# Model against TPXO over open water by depth class (M2 amplitude ratio and complex error)
println("\nM2 sea level against TPXO10 over the domain (every 3rd wet cell), by depth class")
sample = wet[1:3:end]
Ctp = getproperty(tidal_atlas_constants(TPXO10Atlas(), [(λ[i], φ[j]) for (i, j) in sample], :M2; dir = TPXO_DIR), :sea_surface_height)
for (label, lo, hi) in (("shelf < 50 m", 0, 50), ("shelf 50-200 m", 50, 200), ("slope 200-1000 m", 200, 1000), ("deep > 1000 m", 1000, 1e5))
    ks = [k for (k, (i, j)) in enumerate(sample) if lo <= -bh[i, j] < hi && abs(Ctp[k]) > 0]
    isempty(ks) && continue
    Cm = [fit(ηt, η[sample[k]..., :])[:M2] for k in ks]
    r = amp.(Cm) ./ amp.(Ctp[ks])
    @printf("  %-17s %5d cells  median amplitude ratio %.2f  cRMS %.3f m  cRMS / mean amplitude %.2f\n", label, length(ks),
            median(r), sqrt(mean(abs2.(Cm .- Ctp[ks]))), sqrt(mean(abs2.(Cm .- Ctp[ks]))) / mean(amp.(Ctp[ks])))
end

# ---------------- currents at the moorings ----------------
# Semi-major axis of the ellipse traced by (Re(Cu e^{iθ}), Re(Cv e^{iθ}))
semimajor(Cu, Cv) = (wp = abs(Cu + im * Cv) / 2; wm = abs(conj(Cu) + im * conj(Cv)) / 2; wp + wm)
println("\nDepth-mean tidal current semi-major axis (m/s), model / ADCP")
@printf("%-6s %-14s", "", "depths")
for c in FIT; @printf("  %22s", c); end
println()
for name in ("CNSM", "PMCO", "OSSM", "PMUO")
    dataset, Δz, ztop = MOORINGS[name]
    adcp = read_adcp(dataset); m = read_model(PREFIX, name)
    isnothing(m) && continue
    zbot = name in ("CNSM", "PMCO") ? 112.0 : 400.0
    grid_z = collect(ztop:Δz:zbot)
    depthmean(src, A) = [nanmean(regrid(src.depths, A[:, n:n], grid_z)) for n in axes(A, 2)]
    t_a = [Dates.value(t - start_date) / 1000 for t in adcp.times]; t_m = [Dates.value(t - start_date) / 1000 for t in m.times]
    t1 = min(t_a[end], t_m[end]); t0 = max(t_a[1], t_m[1])
    ka = findall(t -> t0 <= t <= t1, t_a); km = findall(t -> t0 <= t <= t1, t_m)
    Fa = (fit(t_a[ka], depthmean(adcp, adcp.U[:, ka])), fit(t_a[ka], depthmean(adcp, adcp.V[:, ka])))
    Fm = (fit(t_m[km], depthmean(m, m.U[:, km])), fit(t_m[km], depthmean(m, m.V[:, km])))
    @printf("%-6s %4.0f-%-4.0f m   ", name, grid_z[1], grid_z[end])
    for c in FIT
        sa = semimajor(Fa[1][c], Fa[2][c]); sm = semimajor(Fm[1][c], Fm[2][c])
        @printf("  %5.3f / %5.3f (%4.2f)", sm, sa, sm / sa)
    end
    println()
end
