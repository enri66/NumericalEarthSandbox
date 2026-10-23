# Nitrate at the Pioneer surface moorings from four sources, Jan-Apr 2019 (mmol/m³):
#   SUNA      the nitrate sensor on the near-surface frame (7 m), monthly means of the hourly data (bgc/download_validation_data.sh)
#   raw       the Copernicus hindcast's monthly NO₃ (GLOBAL_MULTIYEAR_BGC_001_029, 1/4°) at the grid column nearest the mooring, from the file
#   pipeline  the same field as bgc_data.jl builds it for the model (inpainted, interpolated to the mooring), at the middle of each month
#   WOA       the World Ocean Atlas 2023 monthly climatology (µmol/kg × 1.025), same interpolation
# and the profiles of raw, pipeline and WOA in mid-February.
#   julia --project=<env of mpi05_bgc> bgc/nitrate_sources.jl
using NumericalEarth, Oceananigans, Printf, Statistics, Dates
using NumericalEarth.DataWrangling: Metadatum, BoundingBox, NearestNeighborInpainting
using NumericalEarth.WOA: WOAMonthly
using WorldOceanAtlasTools
using CopernicusMarine              # activates the Copernicus download backend (the files are on disk, it only skips them)
ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"
const NCD = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "NCDatasets")]

const DATA = get(ENV, "DATA_DIR", "/t0/workdir/enrique/Data")
const COP_DIR = joinpath(DATA, "NumericalEarth", "copernicus_bgc")
const OOI = joinpath(DATA, "OOI", "pioneer_bgc_2019")
const SITES = ("CNSM" => "cp01cnsm", "ISSM" => "cp03issm")
const COP_REGION = BoundingBox(longitude = (-77.0, -63.0), latitude = (33.0, 43.0))
const WOA_REGION = BoundingBox(longitude = (-78, -62), latitude = (32, 44))
interp(X, f) = Oceananigans.Fields.interpolate(X, f, Oceananigans.instantiated_location(f), f.grid)

# ---- SUNA
function suna(prefix)
    lines = readlines(joinpath(OOI, "$(prefix)_nutnr_all.csv")); col = Dict(h => i for (i, h) in enumerate(split(lines[1], ",")))
    v, q = col["mole_concentration_of_nitrate_in_sea_water_suna"], col["mole_concentration_of_nitrate_in_sea_water_suna_qc_agg"]
    out = Tuple{DateTime, Float64}[]
    for line in lines[3:end]
        f = split(line, ","); f[q] in ("1", "2") || continue
        x = tryparse(Float64, f[v]); (isnothing(x) || !(0 ≤ x ≤ 40)) && continue
        push!(out, (DateTime(f[1][1:19]), x))
    end
    lat = parse(Float64, split(lines[3], ",")[col["latitude"]]); lon = parse(Float64, split(lines[3], ",")[col["longitude"]])
    return out, lon, lat
end
mean_between(data, a, b) = (v = [x for (t, x) in data if a ≤ t < b]; isempty(v) ? NaN : mean(v))

# ---- raw Copernicus column: value at depth z (m, positive) in the grid column nearest (lon, lat)
function raw_column(month, lon, lat, depths)
    file = joinpath(COP_DIR, "no3_CopernicusBGCMonthly_2019-$(lpad(month, 2, '0'))-01T00-00-00_2019-$(lpad(month, 2, '0'))-01T00-00-00_-77.0_-63.0_33.0_43.0.nc")
    ds = NCD.NCDataset(file)
    names_ = Set(keys(ds)); lonv = Float64.(ds[first(filter(n -> n in names_, ("longitude", "lon")))][:]); latv = Float64.(ds[first(filter(n -> n in names_, ("latitude", "lat")))][:])
    z = Float64.(ds["depth"][:]); A = ds["no3"][:, :, :, 1]
    best = (Inf, 0, 0)
    for i in eachindex(lonv), j in eachindex(latv)
        d = hypot((lonv[i] - lon) * cosd(lat), latv[j] - lat); (d < best[1] && !ismissing(A[i, j, 1]) && isfinite(A[i, j, 1])) && (best = (d, i, j))
    end
    _, i, j = best; col = [ismissing(x) ? NaN : Float64(x) for x in A[i, j, :]]
    close(ds)
    vals = map(depths) do d
        k = searchsortedlast(z, d); (k < 1 || k ≥ length(z)) && return NaN
        w = (d - z[k]) / (z[k+1] - z[k]); (1 - w) * col[k] + w * col[k+1]
    end
    bottom = z[something(findlast(isfinite, col), 0)]
    return vals, (lonv[i], latv[j]), bottom
end

# ---- the model's pipeline and WOA
cop_field(month) = Field(Metadatum(:nitrate; dataset = CopernicusBGCMonthly(), region = COP_REGION, date = DateTime(2019, month, 1), dir = COP_DIR), CPU();
                         inpainting = NearestNeighborInpainting(Inf), cache_inpainted_data = false)
woa_field(month) = Field(Metadatum(:nitrate; dataset = WOAMonthly(), region = WOA_REGION, date = DateTime(2018, month, 1)), CPU(); cache_inpainted_data = false)
cop = Dict(m => cop_field(m) for m in 1:4); woa = Dict(m => woa_field(m) for m in 1:4)

for (name, prefix) in SITES
    data, lon, lat = suna(prefix)
    println("\n=== $name at $(round(lon, digits = 3))°E $(round(lat, digits = 3))°N; nitrate at 7 m (mmol/m³)")
    @printf("   SUNA:  Jan 25-31 %.2f   Feb %.2f   Mar %.2f   (first week of Feb %.2f)\n", mean_between(data, DateTime(2019, 1, 25), DateTime(2019, 2, 1)),
            mean_between(data, DateTime(2019, 2, 1), DateTime(2019, 3, 1)), mean_between(data, DateTime(2019, 3, 1), DateTime(2019, 4, 1)),
            mean_between(data, DateTime(2019, 2, 1), DateTime(2019, 2, 8)))
    println("   month     raw Copernicus (cell, bottom)          pipeline   WOA")
    for m in 1:4
        v, cell, bottom = raw_column(m, lon, lat, [7.0])
        @printf("   2019-%02d  %6.2f  (%.3f°E %.3f°N, %4.0f m)     %6.2f   %6.2f\n", m, v[1], cell..., bottom, interp((lon, lat, -7.0), cop[m]), 1.025 * interp((lon, lat, -7.0), woa[m]))
    end
    # what the model starts from: the pipeline linearly between the Jan and Feb frames at Feb 1 (17/31 of the way)
    w = 17 / 31; ic = (1 - w) * interp((lon, lat, -7.0), cop[1]) + w * interp((lon, lat, -7.0), cop[2])
    @printf("   model initial value at 7 m (pipeline, Feb 1): %.2f\n", ic)
    depths = [7.0, 25.0, 50.0, 100.0]
    v, _, _ = raw_column(2, lon, lat, depths)
    println("   mid-Feb profile, depths ", join(Int.(depths), ", "), " m:")
    println("      raw Copernicus  ", join([@sprintf("%6.2f", x) for x in v], " "))
    println("      pipeline        ", join([@sprintf("%6.2f", interp((lon, lat, -d), cop[2])) for d in depths], " "))
    println("      WOA             ", join([@sprintf("%6.2f", 1.025 * interp((lon, lat, -d), woa[2])) for d in depths], " "))
end
