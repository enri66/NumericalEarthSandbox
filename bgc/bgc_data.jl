# BGC data for mab_bgc.jl: the BGC tracers as FieldTimeSeries on the whole-domain source grid `src_grid` (the grid of
# the GLORYS series), for the initial state, the exterior values of the open boundaries and the sponge targets.
# Included by mab_bgc.jl after the GLORYS series are built; uses its `src_grid`, `start_date`, `stop_date`, `data_λ`,
# `data_φ`, `fts_S`, `DATA_DIR`, `comm` and `say`. MAB_BGC_DATA picks the source of NO₃ and O₂ (and P):
#
#   MAB_BGC_DATA=woa (default)
#     NO₃, O₂   WOA monthly climatology (product year 2023, 1°) at mid-month in the top 200 m, the WOA annual field below
#   MAB_BGC_DATA=copernicus
#     NO₃, O₂   Copernicus Marine global BGC hindcast (GLOBAL_MULTIYEAR_BGC_001_029, PISCES, 1/4°, NumericalEarth's
#               CopernicusBGCMonthly): the monthly mean of each month placed at mid-month, as the WOA months are
#     P         the same product's phytoplankton carbon (phyc, mmol C/m³) times the Redfield N:C 16/106
#               (LOBSTER's P is mmol N/m³)
#   both
#     Alk, DIC  placeholders: the Copernicus hindcast has no carbonate system (in v202406 dissic and talk are only in
#               the analysis-forecast product GLOBAL_ANALYSISFORECAST_BGC_001_028, from 2021-11). Alk = 2300 S / 35 and
#               DIC = 0.9 Alk (mmol/m³), from the GLORYS salinity of the frame closest to each mid-month time
#     the rest  uniform values (BGC_UNIFORM), used for the initial state and as the boundary exterior values
#
# Every rank builds the same series: the WOA and Copernicus fields are read for a box wider than the GLORYS box
# (2° for WOA, 1° for Copernicus) and inpainted without writing a cache, so the ranks neither wait on each other nor
# write the same files (see script 05's notes on the bathymetry and GLORYS caches). Missing Copernicus files are
# downloaded first by rank 0 (bgc/bgc_copernicus_download.jl fetches a whole year beforehand).

using NumericalEarth.DataWrangling: Metadatum, metadata_path, NearestNeighborInpainting
import Downloads
using NumericalEarth.WOA: WOAMonthly, WOAAnnual
using WorldOceanAtlasTools          # activates the WOA download backend

ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"   # the WOA files come through DataDeps, which otherwise asks before downloading

const WOA_REGION = BoundingBox(longitude = (data_λ[1] - 2, data_λ[2] + 2), latitude = (data_φ[1] - 2, data_φ[2] + 2))
const WOA_NAMES  = (NO₃ = :nitrate, O₂ = :dissolved_oxygen)
const μmol_per_kg_to_mmol_per_m³ = 1.025   # WOA reports µmol/kg; ρ = 1025 kg/m³

const BGC_DATA = get(ENV, "MAB_BGC_DATA", "woa")
BGC_DATA in ("woa", "copernicus") || error("MAB_BGC_DATA must be woa or copernicus, got $BGC_DATA")

const COPERNICUS_REGION = BoundingBox(longitude = (data_λ[1] - 1, data_λ[2] + 1), latitude = (data_φ[1] - 1, data_φ[2] + 1))
const COPERNICUS_DIR    = joinpath(DATA_DIR, "copernicus_bgc")
const COPERNICUS_NAMES  = (NO₃ = :nitrate, O₂ = :dissolved_oxygen, P = :phytoplankton_carbon)
const COPERNICUS_FACTOR = (NO₃ = 1.0, O₂ = 1.0, P = 16 / 106)   # the product's units are mmol/m³ (P: mmol C → mmol N)

# mmol N/m³ for the nitrogen pools, as LOBSTER's defaults; P comes from the data with MAB_BGC_DATA=copernicus
const BGC_UNIFORM = let uniform = (NH₄ = 0.05, P = 0.1, Z = 0.05, DOM = 0.5, sPOM = 0.01, bPOM = 0.01)
    BGC_DATA == "copernicus" ? Base.structdiff(uniform, NamedTuple{(:P,)}) : uniform
end

# Mid-month times covering the run, in seconds from start_date (the first may be negative)
function bgc_frame_dates(start_date, stop_date)
    first_month = firstdayofmonth(start_date - Day(15))
    last_month  = firstdayofmonth(stop_date + Day(45))
    return [DateTime(year(d), month(d), 15) for d in first_month:Month(1):last_month]
end

woa_field(name, dataset; date_kw...) =
    Field(Metadatum(name; dataset, region = WOA_REGION, date_kw...), CPU(); cache_inpainted_data = false)

# The monthly fields are used down to 200 m, blend linearly into the annual one between 200 and 300 m, and the annual
# field alone is used below: at 1° the monthly nutrient analyses are thinly sampled at depth (Feb nitrate reaches
# 35 mmol/m³ at 350 m at 61.5°W 42.5°N, against at most 25 anywhere in the annual field)
monthly_weight(z) = clamp((z + 300) / 100, 0, 1)

interp(X, f) = Oceananigans.Fields.interpolate(X, f, Oceananigans.instantiated_location(f), f.grid)

function woa_series(tracer, dates)
    name = WOA_NAMES[tracer]
    annual = woa_field(name, WOAAnnual())
    times = [Dates.value(d - start_date) / 1000 for d in dates]
    fts = FieldTimeSeries{Center, Center, Center}(src_grid, times)
    for (n, d) in enumerate(dates)
        monthly = woa_field(name, WOAMonthly(); date = DateTime(2018, month(d), 1))
        A = interior(fts[n])
        for k in axes(A, 3), j in axes(A, 2), i in axes(A, 1)
            X = Oceananigans.Grids.node(i, j, k, src_grid, Center(), Center(), Center())
            w = monthly_weight(X[3])
            value = w > 0 ? w * interp(X, monthly) + (1 - w) * interp(X, annual) : interp(X, annual)
            A[i, j, k] = μmol_per_kg_to_mmol_per_m³ * value
        end
    end
    say(@sprintf("  %s: WOA %s, %d monthly frames, range %.2f to %.2f mmol/m³", tracer, name, length(dates),
                 minimum(interior(fts[1])), maximum(interior(fts[1]))))
    return fts
end

copernicus_metadatum(name, date) =
    Metadatum(name; dataset = CopernicusBGCMonthly(), region = COPERNICUS_REGION, date = firstdayofmonth(date),
              dir = COPERNICUS_DIR)

# NumericalEarth downloads inside `@root` (rank 0 downloads, then every rank waits at a barrier), so every rank must
# call `download` for the same files: the list of missing files is made on every rank before any download starts.
function download_copernicus(dates)
    wanted  = [copernicus_metadatum(name, d) for name in COPERNICUS_NAMES for d in dates]
    absent = filter(md -> !isfile(metadata_path(md)), wanted)
    MPI.Barrier(comm)
    isempty(absent) || say("  downloading $(length(absent)) Copernicus BGC files into $COPERNICUS_DIR…")
    foreach(Downloads.download, absent)
    return nothing
end

function copernicus_series(tracer, dates)
    name = COPERNICUS_NAMES[tracer]
    times = [Dates.value(d - start_date) / 1000 for d in dates]
    fts = FieldTimeSeries{Center, Center, Center}(src_grid, times)
    for (n, d) in enumerate(dates)
        source = Field(copernicus_metadatum(name, d), CPU(); inpainting = NearestNeighborInpainting(Inf),
                       cache_inpainted_data = false)
        A = interior(fts[n])
        for k in axes(A, 3), j in axes(A, 2), i in axes(A, 1)
            X = Oceananigans.Grids.node(i, j, k, src_grid, Center(), Center(), Center())
            A[i, j, k] = COPERNICUS_FACTOR[tracer] * interp(X, source)
        end
    end
    say(@sprintf("  %s: Copernicus BGC monthly %s, %d frames, range %.2f to %.2f mmol/m³", tracer, name, length(dates),
                 minimum(interior(fts[1])), maximum(interior(fts[1]))))
    return fts
end

function carbonate_placeholder_series(dates)
    times = [Dates.value(d - start_date) / 1000 for d in dates]
    Alk = FieldTimeSeries{Center, Center, Center}(src_grid, times)
    DIC = FieldTimeSeries{Center, Center, Center}(src_grid, times)
    for n in eachindex(times)
        m = argmin(abs.(collect(fts_S.times) .- times[n]))
        S = interior(fts_S[m])
        interior(Alk[n]) .= 2300 .* S ./ 35
        interior(DIC[n]) .= 0.9 .* interior(Alk[n])
    end
    say("  Alk, DIC: placeholders from the GLORYS salinity (Alk = 2300 S / 35, DIC = 0.9 Alk)")
    return (; Alk, DIC)
end

function bgc_series()
    dates = bgc_frame_dates(start_date, stop_date)
    say("building the BGC series on the whole-domain source grid ($(length(dates)) frames, " *
        "$(Dates.format(first(dates), "yyyy-mm-dd")) → $(Dates.format(last(dates), "yyyy-mm-dd")))…")
    if BGC_DATA == "copernicus"
        download_copernicus(dates)
        data = NamedTuple(tracer => copernicus_series(tracer, dates) for tracer in keys(COPERNICUS_NAMES))
    else
        data = NamedTuple(tracer => woa_series(tracer, dates) for tracer in keys(WOA_NAMES))
    end
    return merge(data, carbonate_placeholder_series(dates))
end

# DIC (mmol/m³) whose total-scale pH, with alkalinity Alk (meq/m³), temperature T (°C), salinity S and pressure P (bar),
# is pHᵗ: bisection on DIC, since at fixed Alk the pH falls monotonically as DIC rises. For deriving DIC from the PISCES
# pH of the Copernicus reanalysis, which has no DIC.
using OceanBioME.Models.CarbonChemistryModel: CarbonChemistry

function dic_from_pH(pHᵗ, Alk, T, S; P = 0.0, carbon_chemistry = CarbonChemistry(), tolerance = 1e-4)
    pH(DIC) = carbon_chemistry(; DIC, Alk, T, S, P, output = Val(:pHᵗ))
    lo, hi = 0.5 * Alk, 1.2 * Alk
    (pH(lo) ≥ pHᵗ ≥ pH(hi)) || return NaN
    while hi - lo > tolerance
        mid = (lo + hi) / 2
        pH(mid) > pHᵗ ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

# The value of a series at time t on a field's nodes, interpolated in time between its two frames around t
function set_from_series!(field, fts, t)
    times = collect(fts.times)
    n₁ = clamp(searchsortedlast(times, t), 1, length(times))
    n₂ = min(n₁ + 1, length(times))
    w  = n₂ == n₁ ? 0.0 : clamp((t - times[n₁]) / (times[n₂] - times[n₁]), 0, 1)
    ℓ = Oceananigans.instantiated_location(field)
    A = interior(field)
    for k in axes(A, 3), j in axes(A, 2), i in axes(A, 1)
        X = Oceananigans.Grids.node(i, j, k, field.grid, ℓ...)
        A[i, j, k] = (1 - w) * interp(X, fts[n₁]) + w * interp(X, fts[n₂])
    end
    Oceananigans.BoundaryConditions.fill_halo_regions!(field)
    return field
end
