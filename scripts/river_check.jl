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
include(joinpath(@__DIR__, "glofas_land.jl"))

const CPD  = parse(Int, get(ENV, "MAB_CELLS_PER_DEGREE", "12"))
const LAND_FRACTION = parse(Float64, get(ENV, "MAB_LAND_FRACTION", "0.5"))
const CLEANUP = get(ENV, "MAB_MASK_CLEANUP", "true") == "true"
const DATA_DIR = get(ENV, "MAB_DATA_DIR", joinpath(homedir(), "Data", "NumericalEarth"))
const date = DateTime(get(ENV, "RIVER_DATE", "2019-09-15"))
const data_λ, data_φ = (-76.0, -64.0), (34.0, 42.0)
const λ_bounds = (data_λ[1] + 2 / 12, data_λ[2] - 2 / 12); const φ_bounds = (data_φ[1] + 2 / 12, data_φ[2] - 2 / 12)
const Nλ = round(Int, (λ_bounds[2] - λ_bounds[1]) * CPD); const Nφ = round(Int, (φ_bounds[2] - φ_bounds[1]) * CPD)

# 550 levels of 10 m: a cell is wet when its centre is above the bottom, so the surface layer must be thin (script 05 has 1-2 m)
grid = LatitudeLongitudeGrid(CPU(); size = (Nλ, Nφ, 550), longitude = λ_bounds, latitude = φ_bounds, z = (-5500, 0), halo = (7, 7, 7))
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
const EXTRA = get(ENV, "MAB_RIVER_EXTRA", "true") == "true"
if EXTRA                                  # the hand-placed mouths, found as glofas_land_with_mouths finds them
    gλ = collect(λnodes(snapshot.grid, Center())); gφ = collect(φnodes(snapshot.grid, Center()))
    oi = collect(oi); oj = collect(oj); oλ = collect(oλ); oφ = collect(oφ)
    for m in MAB_EXTRA_MOUTHS
        push!(oi, argmin(abs.(gλ .- m.glofas_λ))); push!(oj, argmin(abs.(gφ .- m.glofas_φ))); push!(oλ, m.mouth_λ); push!(oφ, m.mouth_φ)
    end
end
Q = [snapshot[i, j, 1] for (i, j) in zip(oi, oj)]
@printf("%s: %d river mouths in the window, total discharge %.0f m³/s (finite cells %d)\n", Dates.format(date, "yyyy-mm-dd"), length(oi), sum(filter(isfinite, Q)), count(isfinite, interior(snapshot)))

order = sortperm(Q; rev = true)
let g = snapshot.grid, vals = Array(interior(snapshot))[:, :, 1]
    gλ = collect(λnodes(g, Center())); gφ = collect(φnodes(g, Center()))
    big = sort([(vals[i, j], i, j) for i in axes(vals, 1), j in axes(vals, 2) if isfinite(vals[i, j])]; rev = true)
    println("largest GloFAS values anywhere in the window (m³/s, position, is it a mouth cell, ocean (NaN) neighbours):")
    for (q, i, j) in big[1:10]
        nb = count(d -> !isfinite(vals[clamp(i + d[1], 1, end), clamp(j + d[2], 1, end)]), ((1, 0), (-1, 0), (0, 1), (0, -1)))
        @printf("  %9.1f  %8.3f°E %7.3f°N  mouth=%-5s  NaN neighbours %d\n", q, gλ[i], gφ[j], (i, j) in zip(oi, oj), nb)
    end
end
wetcell(λ, φ) = (i = clamp(searchsortedlast(collect(λnodes(grid, Face())), λ), 1, Nλ); j = clamp(searchsortedlast(collect(φnodes(grid, Face())), φ), 1, Nφ); bottom[i, j, 1] < 0)
edge(λ, φ) = λ <= data_λ[1] + 0.1 || λ >= data_λ[2] - 0.1 || φ <= data_φ[1] + 0.1 || φ >= data_φ[2] - 0.1
println("largest mouths (discharge m³/s, position, on the window edge?, mouth cell is model ocean?):")
for n in order[1:min(15, end)]
    @printf("  %9.1f  %8.3f°E %7.3f°N  edge=%-5s  model-ocean=%s\n", Q[n], oλ[n], oφ[n], edge(oλ[n], oφ[n]), wetcell(oλ[n], oφ[n]))
end
println("mouths on the window edge: ", count(n -> edge(oλ[n], oφ[n]), eachindex(oi)), " of ", length(oi), ", carrying ",
        @sprintf("%.0f", sum(Q[n] for n in eachindex(oi) if edge(oλ[n], oφ[n]) && isfinite(Q[n]); init = 0.0)), " m³/s")

const SPREAD = parse(Int, get(ENV, "MAB_RIVER_SPREAD_CELLS", string(round(Int, 8 * CPD / 12))))
println("each river split over up to ", SPREAD, " cells")
land = glofas_land_with_mouths(ibg; extra_mouths = EXTRA ? MAB_EXTRA_MOUTHS : [], start_date = date, end_date = date + Day(1), dir = joinpath(DATA_DIR, "glofas"), region,
                               maximum_spread_cells = SPREAD, maximum_search_radius = round(Int, 5 * CPD / 12))
routing = land.river_routing.rivers
Nt = length(routing.target_i)
println("routing: ", length(routing.contribution_outlet_i), " mouth-to-cell contributions onto ", Nt, " ocean cells")
λc = collect(λnodes(grid, Center())); φc = collect(φnodes(grid, Center()))
using Oceananigans.Operators: Azᶜᶜᶜ
# a contribution k from mouth (oi, oj) onto ocean cell c deposits the volume flux  weight[k] * Q * area(c) / 1000  (m³/s);
# `weight` already includes the share of the mouth sent to that cell if build_river_routing divides it
weight = Float64.(Array(routing.contribution_weight)); coi = Array(routing.contribution_outlet_i); coj = Array(routing.contribution_outlet_j)
tcell = Array(routing.target_i), Array(routing.target_j); offs = Array(routing.offsets)
mouth_index = Dict((oi[m], oj[m]) => m for m in eachindex(oi))
delivered = zeros(Nt); per_mouth = Dict{Int, Vector{Tuple{Int, Float64}}}()
for c in 1:Nt, k in offs[c]:offs[c+1]-1
    m = mouth_index[(coi[k], coj[k])]
    vol = isfinite(Q[m]) ? weight[k] * Q[m] * Azᶜᶜᶜ(tcell[1][c], tcell[2][c], 1, ibg) / 1000 : 0.0
    delivered[c] += vol
    push!(get!(per_mouth, m, Tuple{Int, Float64}[]), (c, vol))
end
println("the three largest mouths and the cells receiving them (volume m³/s, cell centre, depth, distance from the mouth in km):")
for m in order[1:5]
    parts = sort(per_mouth[m]; by = last, rev = true)
    @printf("  mouth %.1f m³/s at %.3f°E %.3f°N -> %d cells, delivered %.1f m³/s\n", Q[m], oλ[m], oφ[m], length(parts), sum(last, parts))
    for (c, vol) in parts[1:min(4, end)]
        i, j = tcell[1][c], tcell[2][c]
        d = 111.0 * sqrt(((λc[i] - oλ[m]) * cosd(oφ[m]))^2 + (φc[j] - oφ[m])^2)
        @printf("      %7.1f m³/s -> cell (%d, %d) at %.3f°E %.3f°N, depth %.0f m, %.0f km from the mouth\n", vol, i, j, λc[i], φc[j], -bottom[i, j, 1], d)
    end
end
@printf("discharge delivered %.0f of %.0f m³/s into %d ocean cells\n", sum(delivered), sum(filter(isfinite, Q)), Nt)
far = [m for m in eachindex(oi) if haskey(per_mouth, m) && isfinite(Q[m]) && Q[m] > 1 && any(x -> 111.0 * sqrt(((λc[tcell[1][x[1]]] - oλ[m]) * cosd(oφ[m]))^2 + (φc[tcell[2][x[1]]] - oφ[m])^2) > 50, per_mouth[m])]
println("mouths over 1 m³/s with a receiving cell more than 50 km away: ", length(far))

# The distributed path: cut the whole-domain routing into a 4 x 2 block layout (as the test runs use) and check that the
# blocks together deliver the same discharge and every destination cell lands in exactly one block.
let nxb = Nλ ÷ 4, nyb = Nφ ÷ 2, total_blocks = 0.0, ncells = 0
    for bi in 0:3, bj in 0:1
        r = localize_routing(routing, bi * nxb, bj * nyb, nxb, nyb, CPU())
        ncells += length(r.target_i)
        toff = Array(r.offsets); w = Array(r.contribution_weight); coi = Array(r.contribution_outlet_i); coj = Array(r.contribution_outlet_j)
        for c in eachindex(r.target_i), k in toff[c]:toff[c+1]-1
            m = mouth_index[(coi[k], coj[k])]
            isfinite(Q[m]) && (total_blocks += w[k] * Q[m] * Azᶜᶜᶜ(r.target_i[c] + bi * nxb, r.target_j[c] + bj * nyb, 1, ibg) / 1000)
        end
    end
    @printf("4 x 2 blocks: %d of %d destination cells kept, discharge delivered %.0f of %.0f m³/s\n", ncells, Nt, total_blocks, sum(delivered))
end
