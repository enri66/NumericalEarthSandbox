# Build the MAB land mask exactly as script 05 does (ETOPO regrid + land fraction), clean it (mask_cleanup.jl), and write
# what a person should look at: a CSV of candidate cells and a figure with the cells the cleanup removed and the candidates
# that are left, for the main estuaries and sounds.
#
#   MAB_CELLS_PER_DEGREE=12 MAB_LAND_FRACTION=0.5 julia --project=. scripts/mask_review.jl
#
# MASK_OUT (default ~/Data/mab_analysis/mask) is the output directory; MASK_OVERRIDES is an overrides CSV (see
# mask_cleanup.jl) applied after the cleanup, e.g. the cells you flagged.
using NumericalEarth, Oceananigans, CairoMakie, Printf, DelimitedFiles
include(joinpath(@__DIR__, "land_fraction_mask.jl"))
include(joinpath(@__DIR__, "mask_cleanup.jl"))

const CPD       = parse(Int, get(ENV, "MAB_CELLS_PER_DEGREE", "12"))
const LAND_FRACTION = parse(Float64, get(ENV, "MAB_LAND_FRACTION", "0.5"))
const OUT       = get(ENV, "MASK_OUT", joinpath(homedir(), "Data", "mab_analysis", "mask"))
const n_pad, src_resolution = 2, 1 / 12                  # as script 05
const data_λ, data_φ = (-76.0, -64.0), (34.0, 42.0)
const λ_bounds = (data_λ[1] + n_pad * src_resolution, data_λ[2] - n_pad * src_resolution)
const φ_bounds = (data_φ[1] + n_pad * src_resolution, data_φ[2] - n_pad * src_resolution)
const Nλ = round(Int, (λ_bounds[2] - λ_bounds[1]) * CPD); const Nφ = round(Int, (φ_bounds[2] - φ_bounds[1]) * CPD)
mkpath(OUT)

grid = LatitudeLongitudeGrid(CPU(); size = (Nλ, Nφ, 2), longitude = λ_bounds, latitude = φ_bounds, z = (-5500, 0), halo = (7, 7, 2))
bottom = regrid_bathymetry(grid; dataset = ETOPO2022(), height_above_water = 1, cache = false,
                           minimum_depth = 10, major_basins = 1, interpolation_passes = 10)
LAND_FRACTION > 0 && apply_land_fraction!(bottom, grid, LAND_FRACTION; minimum_depth = 10)
h = Array(interior(bottom))[:, :, 1]
wet0 = h .< 0
λ = collect(Oceananigans.Grids.λnodes(grid, Oceananigans.Face()))[1:end-1] .+ (λ_bounds[2] - λ_bounds[1]) / Nλ / 2
φ = collect(Oceananigans.Grids.φnodes(grid, Oceananigans.Face()))[1:end-1] .+ (φ_bounds[2] - φ_bounds[1]) / Nφ / 2

overrides = []
if !isempty(get(ENV, "MASK_OVERRIDES", ""))
    for line in readlines(ENV["MASK_OVERRIDES"])
        (isempty(strip(line)) || startswith(line, "#")) && continue
        f = split(strip(line), ",")
        push!(overrides, length(f) == 3 ? (; action = f[1], i = parse(Int, f[2]), j = parse(Int, f[3])) :
                         (; action = f[1], lon0 = parse(Float64, f[2]), lon1 = parse(Float64, f[3]), lat0 = parse(Float64, f[4]), lat1 = parse(Float64, f[5])))
    end
end

wet = clean_mask(wet0; overrides, λ, φ)
connected = ocean_connected(wet0)
removed_isolated = wet0 .& .!connected
removed_dead = wet0 .& connected .& .!wet
candidates = mask_candidates(wet)
@printf("candidates left: %d thin-channel cells, %d one-cell islands, %d corner-only links\n",
        length(candidates.thin_channels), length(candidates.islands), length(candidates.corner_links))
write_candidates(joinpath(OUT, "candidates_$(CPD)cpd.csv"), candidates, λ, φ)

# the cleaned mask as a file the model script can read (1 = wet, 0 = land)
writedlm(joinpath(OUT, "wet_$(CPD)cpd.txt"), Int.(wet))

zooms = (("Chesapeake and Delaware Bays", (-77.4, -74.6), (36.6, 39.8)),
         ("New York, Long Island Sound", (-74.6, -71.4), (40.2, 41.6)),
         ("Cape Cod, Buzzards Bay", (-71.6, -69.6), (41.2, 42.2)),
         ("Pamlico Sound", (-77.0, -75.0), (34.8, 36.2)))
fig = Figure(size = (1400, 1100), fontsize = 13)
Label(fig[0, 1:2], @sprintf("MAB land mask, %d cells/degree, land fraction > %.2f: cleaned wet (blue), removed as isolated (red), removed dead ends (orange), left to review (black: thin channel, magenta: island, green: corner link)", CPD, LAND_FRACTION), fontsize = 14)
for (n, (title, lons, lats)) in enumerate(zooms)
    ax = Axis(fig[div(n - 1, 2) + 1, mod(n - 1, 2) + 1], title = title, aspect = DataAspect())
    ii = findall(x -> lons[1] <= x <= lons[2], λ); jj = findall(y -> lats[1] <= y <= lats[2], φ)
    (isempty(ii) || isempty(jj)) && continue
    code = zeros(length(ii), length(jj))
    for (a, i) in enumerate(ii), (b, j) in enumerate(jj)
        code[a, b] = wet[i, j] ? 1 : removed_isolated[i, j] ? 2 : removed_dead[i, j] ? 3 : 0
    end
    heatmap!(ax, λ[ii], φ[jj], code; colormap = [:tan, :lightskyblue, :red, :orange], colorrange = (0, 3))
    for (cells, color, marker) in ((candidates.thin_channels, :black, :xcross), (candidates.islands, :magenta, :circle), (candidates.corner_links, :green, :diamond))
        pts = [Point2f(λ[i], φ[j]) for (i, j) in cells if λ[i] in lons[1]..lons[2] && φ[j] in lats[1]..lats[2]]
        isempty(pts) || scatter!(ax, pts; color, marker, markersize = 9)
    end
end
save(joinpath(OUT, "mask_review_$(CPD)cpd.png"), fig; px_per_unit = 1.2)
println("wrote ", OUT)
