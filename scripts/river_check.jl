# Where do the GloFAS river mouths in the MAB box go on the model grid? Builds the model's land mask exactly as script 05 does
# (ETOPO regrid, land fraction, optional cleanup), locates the river mouths in the GloFAS files, routes them onto the coastal
# wet cells, and prints the largest rivers with where each one lands, plus the total discharge and any mouth that sits on the
# edge of the downloaded window (a window edge is not a coast).
#
#   MAB_CELLS_PER_DEGREE=24 MAB_LAND_FRACTION=0.5 MAB_MASK_CLEANUP=true julia --project=. scripts/river_check.jl
using NumericalEarth, Oceananigans, Dates, Printf
using CDSAPI      # activates NumericalEarth's GloFAS backend
using NumericalEarth.DataWrangling: Metadata, Metadatum, BoundingBox
using NumericalEarth.Lands: coastal_outlet_indices
include(joinpath(@__DIR__, "land_fraction_mask.jl"))
include(joinpath(@__DIR__, "mask_cleanup.jl"))

const CPD  = parse(Int, get(ENV, "MAB_CELLS_PER_DEGREE", "12"))
const LAND_FRACTION = parse(Float64, get(ENV, "MAB_LAND_FRACTION", "0.5"))
const CLEANUP = get(ENV, "MAB_MASK_CLEANUP", "true") == "true"
const DATA_DIR = get(ENV, "MAB_DATA_DIR", joinpath(homedir(), "Data", "NumericalEarth"))
const date = DateTime(get(ENV, "RIVER_DATE", "2019-09-15"))
const data_λ, data_φ = (-76.0, -64.0), (34.0, 42.0)
const λ_bounds = (data_λ[1] + 2 / 12, data_λ[2] - 2 / 12); const φ_bounds = (data_φ[1] + 2 / 12, data_φ[2] - 2 / 12)
const Nλ = round(Int, (λ_bounds[2] - λ_bounds[1]) * CPD); const Nφ = round(Int, (φ_bounds[2] - φ_bounds[1]) * CPD)

grid = LatitudeLongitudeGrid(CPU(); size = (Nλ, Nφ, 2), longitude = λ_bounds, latitude = φ_bounds, z = (-5500, 0), halo = (7, 7, 2))
bottom = regrid_bathymetry(grid; dataset = ETOPO2022(), height_above_water = 1, cache = false,
                           minimum_depth = 10, major_basins = 1, interpolation_passes = 10)
LAND_FRACTION > 0 && apply_land_fraction!(bottom, grid, LAND_FRACTION; minimum_depth = 10)
if CLEANUP
    h = Array(interior(bottom))[:, :, 1]; wet0 = h .< 0
    wet = clean_mask(wet0)
    h[wet0 .& .!wet] .= 1.0
    set!(bottom, h); Oceananigans.BoundaryConditions.fill_halo_regions!(bottom)
end
ibg = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom))

region = BoundingBox(longitude = data_λ, latitude = data_φ)
meta = Metadata(:river_discharge; dataset = GloFASReanalysis(), start_date = date, end_date = date, dir = joinpath(DATA_DIR, "glofas"), region)
snapshot = Field(first(meta), CPU())
oi, oj, oλ, oφ = coastal_outlet_indices(snapshot)
Q = [snapshot[i, j, 1] for (i, j) in zip(oi, oj)]
@printf("%s: %d river mouths in the window, total discharge %.0f m³/s (finite cells %d)\n", Dates.format(date, "yyyy-mm-dd"), length(oi), sum(filter(isfinite, Q)), count(isfinite, interior(snapshot)))

order = sortperm(Q; rev = true)
wetcell(λ, φ) = (i = clamp(searchsortedlast(collect(λnodes(grid, Face())), λ), 1, Nλ); j = clamp(searchsortedlast(collect(φnodes(grid, Face())), φ), 1, Nφ); bottom[i, j, 1] < 0)
edge(λ, φ) = λ <= data_λ[1] + 0.1 || λ >= data_λ[2] - 0.1 || φ <= data_φ[1] + 0.1 || φ >= data_φ[2] - 0.1
println("largest mouths (discharge m³/s, position, on the window edge?, mouth cell is model ocean?):")
for n in order[1:min(15, end)]
    @printf("  %9.1f  %8.3f°E %7.3f°N  edge=%-5s  model-ocean=%s\n", Q[n], oλ[n], oφ[n], edge(oλ[n], oφ[n]), wetcell(oλ[n], oφ[n]))
end
println("mouths on the window edge: ", count(n -> edge(oλ[n], oφ[n]), eachindex(oi)), " of ", length(oi), ", carrying ",
        @sprintf("%.0f", sum(Q[n] for n in eachindex(oi) if edge(oλ[n], oφ[n]) && isfinite(Q[n]); init = 0.0)), " m³/s")

land = GloFASPrescribedLand(ibg; start_date = date, end_date = date + Day(1), dir = joinpath(DATA_DIR, "glofas"), region)
routing = land.river_routing.rivers
Nt = length(routing.target_i)
println("routing: ", length(routing.contribution_outlet_i), " mouth-to-cell contributions onto ", Nt, " ocean cells")
λc = collect(λnodes(grid, Center())); φc = collect(φnodes(grid, Center()))
weights = Float64.(Array(routing.contribution_weight))
println("ocean cells receiving freshwater, largest first (position of the cell centre):")
total = zeros(Nt)
for c in 1:Nt, k in routing.offsets[c]:routing.offsets[c+1]-1
    n = findfirst(m -> oi[m] == routing.contribution_outlet_i[k] && oj[m] == routing.contribution_outlet_j[k], eachindex(oi))
    isnothing(n) || isfinite(Q[n]) && (total[c] += Q[n] * 1000 / (1000 * 1.0))     # volume flux: weight is the density / cell area
end
for c in sortperm(total; rev = true)[1:min(12, Nt)]
    @printf("  %9.1f m³/s -> cell (%d, %d) at %.3f°E %.3f°N, depth %.0f m\n", total[c], routing.target_i[c], routing.target_j[c],
            λc[routing.target_i[c]], φc[routing.target_j[c]], -bottom[routing.target_i[c], routing.target_j[c], 1])
end
@printf("discharge delivered %.0f of %.0f m³/s\n", sum(total), sum(filter(isfinite, Q)))
