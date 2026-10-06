# Does the model move too much tidal water in and out of Long Island Sound and Block Island Sound? Two estimates of
# the M2 volume flux for one or more script-05 runs (hourly η and MAB_BAROTROPIC_OUTPUT), each against TPXO10:
#   1. continuity: the flux into a basin is iω ∫ η dA, so the M2 elevation fitted at every wet cell of the basin
#      (the run's own nodal factors and phases), times the cell areas, gives the flux amplitude. TPXO10's elevations are
#      taken at the same cells, so both use the model's coastline and the comparison is of the tide itself;
#   2. a section: the M2 zonal transport summed across the eastern entrance of Long Island Sound (the column of u faces
#      nearest SECTION_LON between the coasts), model against TPXO10's zonal transport at the same faces.
# Basins: Long Island Sound (73.75-72.0°W, north of 40.95°N), and Long Island + Block Island Sounds (adding 72.0-71.4°W
# north of 41.0°N). Also the mean M2 amplitude over each basin, model against TPXO, and the gauges there.
# Usage:
#   MAB_TAGS=/t0/.../res_test/cd03,/t0/.../res_test/newdef MAB_START_DATE=2019-08-29 \
#   TPXO_DIR=/t0/workdir/enrique/Data/TPXO10_atlas_v2_nc julia --project=. scripts/lis_transport.jl
include(joinpath(@__DIR__, "mab_analysis_common.jl"))
using NumericalEarth: TPXO10Atlas, earth_tidal_harmonics
using Oceananigans: tidal_atlas_constants
using Oceananigans.Operators: Azᶜᶜᶜ

const TAGS        = filter(!isempty, split(get(ENV, "MAB_TAGS", ""), ","))
const TPXO_DIR    = get(ENV, "TPXO_DIR", joinpath(homedir(), "Data", "TPXO10_atlas_v2_nc"))
const SKIP_DAYS   = parse(Float64, get(ENV, "TIDE_SKIP_DAYS", "2"))
const SECTION_LON = parse(Float64, get(ENV, "SECTION_LON", "-72.05"))
const FIT         = (:M2, :S2, :N2, :K1, :O1)
const BASINS      = (("Long Island Sound", [((-73.75, -72.0), (40.95, 41.5))]),
                     ("Long Island + Block Island Sounds", [((-73.75, -72.0), (40.95, 41.5)), ((-72.0, -71.4), (41.0, 41.5))]))

harmonics = earth_tidal_harmonics(start_date; constituents = FIT, ramp_time = 0)
const ωM2 = harmonics.frequencies[1]
function fit_M2(t, X)
    keep = t .>= SKIP_DAYS * 86400
    t = t[keep]; X = X[keep, :]
    cols = [ones(length(t))]
    for n in eachindex(FIT)
        θ = harmonics.frequencies[n] .* t .+ harmonics.phases[n]
        push!(cols, harmonics.nodal_factors[n] .* cos.(θ)); push!(cols, -harmonics.nodal_factors[n] .* sin.(θ))
    end
    β = reduce(hcat, cols) \ X
    return complex.(β[2, :], β[3, :])
end

ref = run_prefix(String(first(TAGS)))
grid = global_grid(ref * "_barotropic.jld2"); ug = grid.underlying_grid
Nx, Ny = size(ug, 1), size(ug, 2)
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center())); λf = collect(λnodes(ug, Face()))
B = grid.immersed_boundary.bottom_height; wet = [B[i, j, 1] < 0 for i in 1:Nx, j in 1:Ny]
inbox((λr, φr), i, j) = λr[1] <= λ[i] <= λr[2] && φr[1] <= φ[j] <= φr[2]
basin_cells(boxes) = [(i, j) for i in 1:Nx, j in 1:Ny if wet[i, j] && any(b -> inbox(b, i, j), boxes)]

# Section: the u-face column nearest SECTION_LON, wet faces between 40.95°N and 41.5°N
is = argmin(abs.(λf .- SECTION_LON))
section = [j for j in 1:Ny if 40.95 <= φ[j] <= 41.5 && is > 1 && is <= Nx && wet[is-1, j] && wet[is, j]]
dy = [6.371e6 * deg2rad(φ[2] - φ[1]) for _ in section]
@printf("Section at %.3f°E: %d wet u faces between %.2f and %.2f°N\n", λf[is], length(section),
        isempty(section) ? NaN : φ[first(section)], isempty(section) ? NaN : φ[last(section)])

tpxo_η(cells) = getproperty(tidal_atlas_constants(TPXO10Atlas(), [(λ[i], φ[j]) for (i, j) in cells], :M2; dir = TPXO_DIR), :sea_surface_height)
tpxo_U_section = isempty(section) ? ComplexF64[] :
    getproperty(tidal_atlas_constants(TPXO10Atlas(), [(λf[is], φ[j]) for j in section], :M2; dir = TPXO_DIR), :eastward_transport)

for tag in TAGS
    prefix = run_prefix(String(tag))
    η, tη = surface_series(prefix * "_eta.jld2", "η")
    println("\n== $(basename(String(tag)))")
    for (name, boxes) in BASINS
        cells = basin_cells(boxes)
        A = [Azᶜᶜᶜ(i, j, 1, ug) for (i, j) in cells]
        Cm = fit_M2(collect(tη), reduce(hcat, [η[i, j, :] for (i, j) in cells]))
        Ct = tpxo_η(cells)
        ok = abs.(Ct) .> 0                                     # TPXO cells it counts as ocean
        Qm = ωM2 * abs(sum(Cm[ok] .* A[ok])); Qt = ωM2 * abs(sum(Ct[ok] .* A[ok]))
        @printf("   %-34s %4d cells (%4d in TPXO), area %6.0f km²: M2 flux model %7.0f, TPXO %7.0f m³/s (×%.2f); mean amplitude model %.2f, TPXO %.2f m\n",
                name, length(cells), count(ok), sum(A[ok]) / 1e6, Qm, Qt, Qm / Qt, mean(abs.(Cm[ok])), mean(abs.(Ct[ok])))
    end
    if !isempty(section)
        U, tU = surface_series(prefix * "_barotropic.jld2", "U")
        Cu = fit_M2(collect(tU), reduce(hcat, [U[is, j, :] for j in section]))
        Qm = abs(sum(Cu .* dy)); Qt = abs(sum(tpxo_U_section .* dy))
        @printf("   section at %.2f°E (eastern Long Island Sound): M2 flux model %7.0f, TPXO %7.0f m³/s (×%.2f)\n", λf[is], Qm, Qt, Qm / Qt)
    end
end
println("\nGauges (NOAA, M2 amplitude): New London 0.36 m, Montauk 0.28 m; TPXO at those points matches them (tide_check.jl)")
