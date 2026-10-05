# M2 tidal-current ellipses of a run's barotropic transport against TPXO10 over the shelf: in the band where the model's
# semidiurnal currents are too strong (southern New England / New York Bight outer shelf), is the excess in the
# cross-shelf or the along-shelf flow, and are the ellipses turned or shifted in phase? For each shelf cell (depth
# 20-200 m) in two regions — the band and the southern MAB shelf as a reference — the M2 transport (U, V) of the run
# (hourly MAB_BAROTROPIC_OUTPUT, fitted with the run's own nodal factors and phases) and of TPXO10 give
#   the semi-major and semi-minor axes and the inclination of the major axis (degrees from east, 0-180), found as the
#   direction of largest and smallest component amplitude, and the Greenwich phase of the major-axis component;
#   the amplitudes of the cross-shelf (along ∇h, the local bathymetric gradient) and along-shelf components.
# Reported per region: medians of model / TPXO for each, of the inclination difference and of the phase difference.
# Usage:
#   MAB_TAG=/t0/.../res_test/cd03 MAB_START_DATE=2019-08-29 TPXO_DIR=/t0/workdir/enrique/Data/TPXO10_atlas_v2_nc \
#       julia --project=. scripts/tidal_ellipses.jl
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))
using NumericalEarth: TPXO10Atlas, earth_tidal_harmonics
using Oceananigans: tidal_atlas_constants

const TAG       = get(ENV, "MAB_TAG", "cd03")
const PREFIX    = run_prefix(TAG)
const TPXO_DIR  = get(ENV, "TPXO_DIR", joinpath(homedir(), "Data", "TPXO10_atlas_v2_nc"))
const SKIP_DAYS = parse(Float64, get(ENV, "TIDE_SKIP_DAYS", "2"))
const FIT       = (:M2, :S2, :N2, :K1, :O1)
const REGIONS   = (("band 70-74W 39.5-41N", (-74.0, -70.0), (39.5, 41.0)), ("south 74-76W 35-39N", (-76.0, -74.0), (35.0, 39.0)))

harmonics = earth_tidal_harmonics(start_date; constituents = FIT, ramp_time = 0)
function fit_many(t, X)
    keep = t .>= SKIP_DAYS * 86400
    t = t[keep]; X = X[keep, :]
    cols = [ones(length(t))]
    for n in eachindex(FIT)
        θ = harmonics.frequencies[n] .* t .+ harmonics.phases[n]
        push!(cols, harmonics.nodal_factors[n] .* cos.(θ)); push!(cols, -harmonics.nodal_factors[n] .* sin.(θ))
    end
    β = reduce(hcat, cols) \ X
    return complex.(β[2, :], β[3, :])                          # M2 only (first constituent)
end

# Component of the ellipse along direction α (radians from east): complex amplitude Cu cos α + Cv sin α
component(Cu, Cv, α) = Cu * cos(α) + Cv * sin(α)
function ellipse(Cu, Cv)
    αs = range(0, π; length = 361)[1:end-1]
    a = [abs(component(Cu, Cv, α)) for α in αs]
    k = argmax(a)
    major = a[k]; minor = abs(component(Cu, Cv, αs[k] + π / 2))
    phase = mod(rad2deg(-angle(component(Cu, Cv, αs[k]))), 360)
    return (; major, minor, inclination = rad2deg(αs[k]), phase)
end
angdiff(a, b; period = 180) = mod(a - b + period / 2, period) - period / 2

grid = isfile(PREFIX * "_barotropic.jld2") ? JLD2.jldopen(f -> f["serialized/grid"], PREFIX * "_barotropic.jld2") :
                                             global_grid(PREFIX * "_barotropic.jld2")
ug = grid.underlying_grid
Nx, Ny = size(ug, 1), size(ug, 2)
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center()))
B = grid.immersed_boundary.bottom_height; h = [-B[i, j, 1] for i in 1:Nx, j in 1:Ny]
R = 6.371e6; dy = R * deg2rad(φ[2] - φ[1])

U, tU = surface_series(PREFIX * "_barotropic.jld2", "U")
V, _  = surface_series(PREFIX * "_barotropic.jld2", "V")

rows = []
for (label, λr, φr) in REGIONS
    cells = [(i, j) for i in 3:Nx-2, j in 3:Ny-2 if λr[1] <= λ[i] <= λr[2] && φr[1] <= φ[j] <= φr[2] && 20 <= h[i, j] <= 200 &&
             all(h[i+a, j+b] > 0 for a in -2:2, b in -2:2)]
    isempty(cells) && continue
    Cu = fit_many(collect(tU), reduce(hcat, [(U[i, j, :] .+ U[i+1, j, :]) ./ 2 for (i, j) in cells]))
    Cv = fit_many(collect(tU), reduce(hcat, [(V[i, j, :] .+ V[i, j+1, :]) ./ 2 for (i, j) in cells]))
    nodes = [(λ[i], φ[j]) for (i, j) in cells]
    Tu = getproperty(tidal_atlas_constants(TPXO10Atlas(), nodes, :M2; dir = TPXO_DIR), :eastward_transport)
    Tv = getproperty(tidal_atlas_constants(TPXO10Atlas(), nodes, :M2; dir = TPXO_DIR), :northward_transport)
    for (k, (i, j)) in enumerate(cells)
        # cross-shelf direction: the depth gradient, centred over two cells each way (smoother than one)
        dx = R * cosd(φ[j]) * deg2rad(λ[2] - λ[1])
        gx = (h[i+2, j] - h[i-2, j]) / (4dx); gy = (h[i, j+2] - h[i, j-2]) / (4dy)
        hypot(gx, gy) > 0 || continue
        α = atan(gy, gx)
        em, et = ellipse(Cu[k], Cv[k]), ellipse(Tu[k], Tv[k])
        et.major > 0 || continue
        push!(rows, (; region = label, λ = λ[i], φ = φ[j], h = h[i, j],
                       major = em.major / et.major, minor_m = em.minor, minor_t = et.minor,
                       dinc = angdiff(em.inclination, et.inclination), dphase = angdiff(em.phase, et.phase; period = 360),
                       cross = abs(component(Cu[k], Cv[k], α)) / abs(component(Tu[k], Tv[k], α)),
                       along = abs(component(Cu[k], Cv[k], α + π / 2)) / abs(component(Tu[k], Tv[k], α + π / 2)),
                       cross_share_t = abs(component(Tu[k], Tv[k], α)) / et.major,
                       inc_m = em.inclination, inc_t = et.inclination))
    end
end

col(rs, q) = [getproperty(r, q) for r in rs]
med(x) = median(filter(isfinite, x))
println("M2 barotropic transport ellipses, $(basename(PREFIX)) against TPXO10, shelf cells 20-200 m (medians)")
@printf("%-22s %6s %10s %10s %10s %12s %12s %16s\n", "region", "cells", "major", "cross", "along", "Δincl (°)", "Δphase (°)",
        "TPXO cross share")
for (label, _, _) in REGIONS
    rs = filter(r -> r.region == label, rows)
    isempty(rs) && continue
    @printf("%-22s %6d %10.2f %10.2f %10.2f %12.1f %12.1f %16.2f\n", label, length(rs), med(col(rs, :major)), med(col(rs, :cross)),
            med(col(rs, :along)), med(col(rs, :dinc)), med(col(rs, :dphase)), med(col(rs, :cross_share_t)))
end
println("(major, cross, along: model / TPXO amplitude; Δ: model − TPXO; TPXO cross share: |cross-shelf| / semi-major of TPXO)")

# ---------------- figure ----------------
fig = Figure(size = (1500, 520), fontsize = 13)
Label(fig[0, 1:3], "$(basename(PREFIX)) M2 transport ellipses against TPXO10, shelf 20-200 m", fontsize = 16)
for (c, (q, title, cr, cmap, f)) in enumerate(((:cross, "cross-shelf amplitude, model / TPXO", (-2, 2), :balance, x -> log2(x)),
                                               (:along, "along-shelf amplitude, model / TPXO", (-2, 2), :balance, x -> log2(x)),
                                               (:dinc, "inclination, model − TPXO (°)", (-45, 45), :balance, identity)))
    ax = Axis(fig[1, c], title = title, aspect = DataAspect())
    s = scatter!(ax, col(rows, :λ), col(rows, :φ); color = clamp.(f.(col(rows, q)), cr...), colorrange = cr, colormap = cmap,
                 markersize = 6, marker = :rect)
    Colorbar(fig[2, c], s; vertical = false, label = c < 3 ? "log₂ ratio" : "degrees")
end
out = PREFIX * "_tidal_ellipses.png"
save(out, fig; px_per_unit = 1.2)
println("saved ", out)
