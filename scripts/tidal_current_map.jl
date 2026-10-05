# Tidal currents over the whole domain against TPXO10: where are the model's barotropic tidal currents too strong?
# Reads one or more script-05 runs with MAB_BAROTROPIC_OUTPUT=true (hourly split-explicit transports U, V, m²/s),
# fits M2 S2 N2 K1 O1 at every STRIDE-th wet cell (cell-centred averages of U and V; the run's own nodal factors and
# phases, first TIDE_SKIP_DAYS left out) and compares the tidal-ellipse semi-major axis of the TRANSPORT with
# TPXO10's at the same points: median model / TPXO by depth class for each constituent, and maps of the M2 ratio.
# Usage:
#   MAB_TAGS=/t0/.../res_test/cd03,/t0/.../res_test/cd06 MAB_START_DATE=2019-08-29 \
#   TPXO_DIR=/t0/workdir/enrique/Data/TPXO10_atlas_v2_nc julia --project=. scripts/tidal_current_map.jl
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))
using NumericalEarth: TPXO10Atlas, earth_tidal_harmonics
using Oceananigans: tidal_atlas_constants

const TAGS      = filter(!isempty, split(get(ENV, "MAB_TAGS", ""), ","))
const TPXO_DIR  = get(ENV, "TPXO_DIR", joinpath(homedir(), "Data", "TPXO10_atlas_v2_nc"))
const SKIP_DAYS = parse(Float64, get(ENV, "TIDE_SKIP_DAYS", "2"))
const STRIDE    = parse(Int, get(ENV, "MAP_STRIDE", "2"))
const FIT       = (:M2, :S2, :N2, :K1, :O1)
const CLASSES   = (("shelf < 50 m", 0, 50), ("shelf 50-200 m", 50, 200), ("slope 200-1000 m", 200, 1000), ("deep > 1000 m", 1000, 1e5))

harmonics = earth_tidal_harmonics(start_date; constituents = FIT, ramp_time = 0)

# Complex amplitudes for all FIT constituents at once, for many series (columns of X) sharing one time axis
function fit_many(t, X)
    keep = t .>= SKIP_DAYS * 86400
    t = t[keep]; X = X[keep, :]
    cols = [ones(length(t))]
    for n in eachindex(FIT)
        θ = harmonics.frequencies[n] .* t .+ harmonics.phases[n]
        push!(cols, harmonics.nodal_factors[n] .* cos.(θ)); push!(cols, -harmonics.nodal_factors[n] .* sin.(θ))
    end
    A = reduce(hcat, cols)
    β = A \ X
    return Dict(FIT[n] => complex.(β[2n, :], β[2n+1, :]) for n in eachindex(FIT))
end
semimajor(Cu, Cv) = (abs(Cu + im * Cv) + abs(conj(Cu) + im * conj(Cv))) / 2

# Grid, wet cells and depths from the first run's barotropic output (its rank files carry the grid)
ref = run_prefix(String(first(TAGS)))
grid = isfile(ref * "_barotropic.jld2") ? JLD2.jldopen(f -> f["serialized/grid"], ref * "_barotropic.jld2") :
                                          global_grid(ref * "_barotropic.jld2")
ug = grid.underlying_grid
Nx, Ny = size(ug, 1), size(ug, 2)
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center()))
B = grid.immersed_boundary.bottom_height; bh = [B[i, j, 1] for i in 1:Nx, j in 1:Ny]
cells = [(i, j) for i in 2:STRIDE:Nx-1, j in 2:STRIDE:Ny-1 if bh[i, j] < 0 && bh[i-1, j] < 0 && bh[i+1, j] < 0 && bh[i, j-1] < 0 && bh[i, j+1] < 0]
depth = [-bh[i, j] for (i, j) in cells]
nodes = [(λ[i], φ[j]) for (i, j) in cells]
@printf("%d interior wet cells (every %d-th), TPXO10 at the same points\n", length(cells), STRIDE)
tpxo = Dict(c => (getproperty(tidal_atlas_constants(TPXO10Atlas(), nodes, c; dir = TPXO_DIR), :eastward_transport),
                  getproperty(tidal_atlas_constants(TPXO10Atlas(), nodes, c; dir = TPXO_DIR), :northward_transport)) for c in FIT)
tpxo_major = Dict(c => semimajor.(tpxo[c]...) for c in FIT)

results = Dict{String, Any}()
for tag in TAGS
    prefix = run_prefix(String(tag))
    U, tU = surface_series(prefix * "_barotropic.jld2", "U")
    V, _  = surface_series(prefix * "_barotropic.jld2", "V")
    Uc = reduce(hcat, [(U[i, j, :] .+ U[i+1, j, :]) ./ 2 for (i, j) in cells])     # (time, cell)
    Vc = reduce(hcat, [(V[i, j, :] .+ V[i, j+1, :]) ./ 2 for (i, j) in cells])
    Fu = fit_many(collect(tU), Uc); Fv = fit_many(collect(tU), Vc)
    major = Dict(c => semimajor.(Fu[c], Fv[c]) for c in FIT)
    results[basename(String(tag))] = major
    println("\n== $(basename(String(tag))): median model / TPXO transport semi-major axis by depth class (cells)")
    @printf("   %-17s %6s", "", "cells")
    for c in FIT; @printf(" %7s", c); end
    println()
    for (label, lo, hi) in CLASSES
        k = findall(d -> lo <= d < hi, depth)
        k = filter(m -> tpxo_major[:M2][m] > 0, k)
        isempty(k) && continue
        @printf("   %-17s %6d", label, length(k))
        for c in FIT
            r = major[c][k] ./ tpxo_major[c][k]
            @printf(" %7.2f", median(filter(isfinite, r)))
        end
        println()
    end
end

# ---------------- figure: M2 transport semi-major, TPXO and model / TPXO ----------------
labels = sort(collect(keys(results)))
fig = Figure(size = (520 * (1 + length(labels)), 520), fontsize = 13)
Label(fig[0, 1:1+length(labels)], "M2 barotropic transport semi-major axis (m²/s): TPXO10 and model / TPXO10", fontsize = 16)
ax = Axis(fig[1, 1], title = "TPXO10", aspect = DataAspect())
s = scatter!(ax, first.(nodes), last.(nodes); color = log10.(max.(tpxo_major[:M2], 1e-3)), colorrange = (-2, 1.5),
             colormap = :viridis, markersize = 5, marker = :rect)
Colorbar(fig[2, 1], s; vertical = false, label = "log₁₀ semi-major (m²/s)")
for (c, label) in enumerate(labels)
    ax = Axis(fig[1, 1 + c], title = "$label / TPXO", aspect = DataAspect())
    r = results[label][:M2] ./ tpxo_major[:M2]
    s = scatter!(ax, first.(nodes), last.(nodes); color = log2.(clamp.(r, 0.25, 4)), colorrange = (-2, 2), colormap = :balance,
                 markersize = 5, marker = :rect)
    Colorbar(fig[2, 1 + c], s; vertical = false, label = "log₂ ratio")
end
out = joinpath(dirname(ref), get(ENV, "MAP_FIGURE", "tidal_current_map.png"))
save(out, fig; px_per_unit = 1.2)
println("\nsaved ", out)
