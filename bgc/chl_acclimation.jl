# Does a variable chlorophyll:nitrogen ratio bring the model's surface chlorophyll toward the satellite? A diagnostic on the output of
# a run of bgc/mab_bgc.jl with MAB_BGC_LIMITATION=true; the run itself keeps LOBSTER's fixed 1.31 mg Chl/mmol N. Matched cell-days
# (surface cell, model day with a valid satellite value) are binned by water depth as in validate_bgc.jl, and for each bin:
#   fixed       the model's chlorophyll with 1.31 mg Chl/mmol N
#   required    the ratio that would make the mean model chlorophyll equal the satellite mean: mean(satellite) / mean(P), mg Chl/mmol N
#   Cloern      the model's P times a photoacclimated ratio, Cloern et al. (1995, Limnol. Oceanogr. 40, 1313):
#                 Chl:C (g/g) = 0.003 + 0.0154 exp(0.050 T) exp(-0.059 I) f_N
#               T the temperature (°C) of the surface cell, I the daily irradiance (mol photons/m²/day) averaged over the mixed layer
#               (the saved levels down to the mixed-layer depth, PAR in W/m² × 4.57 µmol/J, the daily mean of the model's PAR),
#               f_N = min(1, 2 L_N) the nutrient-limited fraction of the growth rate (L_N = (L_NO₃ + L_NH₄) / 2 of bgc/bgc_limitation.jl,
#               whose maximum is 1/2); converted with the C:N of PhytoZoo, 6.56 mol/mol, and 12 mg C/mmol: Chl:N = 12 · 6.56 · Chl:C
#   Chl:N bounds  the range of the ratio for Chl:C between 0.005 and 0.06 g/g (about 0.4-4.7 mg Chl/mmol N), the usual extremes
# Neither the satellite's coastal bias (CDOM, sediment) nor the lag of the cells' acclimation to the mixing is represented; the
# Cloern ratio is a steady state for the mixed-layer light of the day.
#   julia -t 4 --project=<env with CairoMakie> bgc/chl_acclimation.jl
# Environment: BGC_TAG (run path prefix, needs *_bgc_limitation_rank*, *_bgc_volume_daily_rank*, *_bottom_rank*), MAB_START_DATE,
# MAB_CELLS_PER_DEGREE, MAB_NZ, BGC_SAT (satellite NetCDF), BGC_OUT (summary.txt).
using Printf, Statistics, Dates
using Oceananigans
using Oceananigans.Grids: znodes
using NumericalEarth.DataWrangling: jldopen
const NCD = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "NCDatasets")]

const TAG   = get(ENV, "BGC_TAG", "/t0/workdir/enrique/runs/bgc/bgc_feb19_woa_lim_mg30_fix")
const START = DateTime(get(ENV, "MAB_START_DATE", "2019-02-01"))
const CPD   = parse(Int, get(ENV, "MAB_CELLS_PER_DEGREE", "12"))
const SAT   = get(ENV, "BGC_SAT", "/t0/workdir/enrique/Data/ocean_colour/chl_l3_multi_4km_mab_2019.nc")
const OUT   = get(ENV, "BGC_OUT", "/t0/workdir/enrique/runs/analysis/acclimation_" * basename(TAG))
const H, Nz = 7, parse(Int, get(ENV, "MAB_NZ", "50"))
const CHL_PER_N = 1.31                                   # mg Chl per mmol N (OceanBioME's PhytoZoo)
const CARBON_RATIO = 6.56                                # mol C / mol N (PhytoZoo)
const PAR_TO_PHOTONS = 4.57e-6 * 86400                   # W/m² (PAR, daily mean) → mol photons/m²/day
const MLD_THRESHOLD, MLD_REFERENCE = 0.03, 10.0
const λ_west, φ_south = -76 + 2 / 12, 34 + 2 / 12        # script 05's box: GLORYS box edge plus two GLORYS cells
mkpath(OUT)
summary_io = open(joinpath(OUT, "summary.txt"), "w")
function say(args...)
    println(args...); flush(stdout)
    println(summary_io, args...); flush(summary_io)
end

# ---------------- the model's layout, depth and vertical grid (as in validate_bgc.jl) ----------------
function exponential_scale(N, H, Δtop)
    top(h) = H * expm1(H / N / h) / expm1(H / h)
    lo, hi = H / 50, 100H
    for _ in 1:200; mid = (lo + hi) / 2; top(mid) > Δtop ? (hi = mid) : (lo = mid); end
    return (lo + hi) / 2
end
Δz_top, Z_BOTTOM = parse(Float64, get(ENV, "MAB_DZ_TOP", "2")), parse(Float64, get(ENV, "MAB_ZBOTTOM", "5500"))
zgrid = RectilinearGrid(size = (1, 1, Nz), x = (0, 1), y = (0, 1),
                        z = ExponentialDiscretization(Nz, -Z_BOTTOM, 0; scale = exponential_scale(Nz, Z_BOTTOM, Δz_top)))
const zc = collect(znodes(zgrid, Center())); const zf = collect(znodes(zgrid, Face())); const Δz = diff(zf)   # bottom first

bottom_files = sort(filter(f -> occursin(Regex("^" * basename(TAG) * "_bottom_rank\\d+\\.jld2\$"), basename(f)), readdir(dirname(TAG); join = true)))
NR = length(bottom_files)
NR > 0 || error("no $(basename(TAG))_bottom_rank*.jld2 in $(dirname(TAG))")
parts = map(bottom_files) do f
    jldopen(f) do file
        (I = file["I_OFF"], J = file["J_OFF"], kb = file["k_bottom"], dzb = file["Δz_bottom"])
    end
end
const Nx = maximum(p.I + size(p.kb, 1) for p in parts); const Ny = maximum(p.J + size(p.kb, 2) for p in parts)
kb = zeros(Int, Nx, Ny); depth = zeros(Nx, Ny)
for p in parts, j in axes(p.kb, 2), i in axes(p.kb, 1)
    k = p.kb[i, j]
    kb[p.I + i, p.J + j] = k
    depth[p.I + i, p.J + j] = k > 0 ? sum(Δz[k+1:end]) + p.dzb[i, j] : 0.0
end
wet = kb .> 0
λm(i) = λ_west + (i - 0.5) / CPD
φm(j) = φ_south + (j - 0.5) / CPD
say("model: $(Nx)×$(Ny)×$(Nz) on $NR ranks, $(count(wet)) wet columns")

iters(f, name) = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/$name"]))))

# the volume frames (P, T, S): full columns; the limitation frames (PAR, nitrate and ammonium terms): the top L levels
function volume_frame(name, n; stem = "bgc_volume_daily")
    A = zeros(Float32, Nx, Ny, Nz)
    for (p, file) in zip(parts, ["$(TAG)_$(stem)_rank$(r).jld2" for r in 0:NR-1])
        jldopen(file) do f
            it = iters(f, name)[n]
            B = f["timeseries/$name/$it"][H+1:end-H, H+1:end-H, H+1:end-H]
            A[p.I+1:p.I+size(B, 1), p.J+1:p.J+size(B, 2), :] .= B
        end
    end
    return A
end
volume_files = ["$(TAG)_bgc_volume_daily_rank$(r).jld2" for r in 0:NR-1]
vol_times = jldopen(f -> [f["timeseries/t/$it"] for it in iters(f, "P")], volume_files[1])
vol_dates = [START + Millisecond(round(Int, 1000t)) for t in vol_times]

limitation_files = ["$(TAG)_bgc_limitation_rank$(r).jld2" for r in 0:NR-1]
lim_times = jldopen(f -> [f["timeseries/t/$it"] for it in iters(f, "PAR")], limitation_files[1])
lim_dates = [START + Millisecond(round(Int, 1000t)) for t in lim_times]
const L = jldopen(f -> size(f["timeseries/PAR/$(iters(f, "PAR")[1])"], 3), limitation_files[1])
const klev = Nz-L+1:Nz                                   # model levels of the saved ones, bottom first
function limitation_frame(name, n)
    A = fill(NaN32, Nx, Ny, L)
    for (p, file) in zip(parts, limitation_files)
        jldopen(file) do f
            B = f["timeseries/$name/$(iters(f, name)[n])"]
            nx, ny = size(p.kb)
            size(B, 1) == nx + 2H && (B = B[H+1:end-H, H+1:end-H, :])
            A[p.I+1:p.I+nx, p.J+1:p.J+ny, :] .= B
        end
    end
    return A
end
say("volume frames $(length(vol_dates)) ($(Date(first(vol_dates))) → $(Date(last(vol_dates)))), limitation frames $(length(lim_dates)) ($(Date(first(lim_dates))) → $(Date(last(lim_dates)))); limitation levels: top $L, down to $(round(-zf[Nz-L+1], digits = 1)) m")
lim_index = Dict(Date(d) => n for (n, d) in enumerate(lim_dates))

# mixed-layer depth from the daily T and S, as in bgc/limitation_analysis.jl: first depth at which the density exceeds that at 10 m by
# 0.03 kg/m³ (simplified equation of state of Roquet et al. 2015 at the surface pressure; mixed to the seafloor where it does not exceed it)
function density_anomaly(Θ, S)
    Θ′, S′ = Θ - 10, S - 35
    return -0.1655 * (1 + 0.5 * 0.05952 * Θ′) * Θ′ + 0.76554 * (1 - 0.5 * 0.00054914 * S′) * S′ - 0.0024341 * Θ′ * S′
end
function mixed_layer_depth(T, S)
    mld = fill(NaN, Nx, Ny)
    for j in 1:Ny, i in 1:Nx
        k0 = kb[i, j]; k0 > 0 || continue
        kref = Nz
        while kref > k0 && -zc[kref] < MLD_REFERENCE; kref -= 1; end
        σref = density_anomaly(T[i, j, kref], S[i, j, kref])
        crossed = false
        for k in kref-1:-1:k0
            σ = density_anomaly(T[i, j, k], S[i, j, k])
            if σ - σref ≥ MLD_THRESHOLD
                σ₊ = density_anomaly(T[i, j, k+1], S[i, j, k+1])
                w = σ == σ₊ ? 0.0 : (σref + MLD_THRESHOLD - σ₊) / (σ - σ₊)
                mld[i, j] = -zc[k+1] + w * (zc[k+1] - zc[k]); crossed = true; break
            end
        end
        crossed || (mld[i, j] = depth[i, j])
    end
    return mld
end
upper = -zf[klev .+ 1]; lower = -zf[klev]                # upper and lower face depth of each saved level (m, positive down)
# mean of a saved quantity A[i, j, l] over the mixed layer, weighted by the thickness of the levels (over the saved levels only
# where the mixed layer is deeper than them)
function mixed_layer_mean(A, mld)
    out = fill(NaN, Nx, Ny)
    for j in 1:Ny, i in 1:Nx
        (isfinite(mld[i, j]) && kb[i, j] > 0) || continue
        s = 0.0; w = 0.0
        for l in 1:L
            Δ = clamp(mld[i, j] - upper[l], 0.0, lower[l] - upper[l])
            isfinite(A[i, j, l]) && (s += Δ * A[i, j, l]; w += Δ)
        end
        w > 0 && (out[i, j] = s / w)
    end
    return out
end

# ---------------- satellite ----------------
ds = NCD.Dataset(SAT)
λs = Float64.(ds["longitude"][:]); φs = Float64.(ds["latitude"][:]); sat_times = ds["time"][:]
raw = ds["CHL"][:, :, :]
chl_sat = Array{Float32}(undef, size(raw)); for k in eachindex(raw); chl_sat[k] = ismissing(raw[k]) ? NaN32 : Float32(raw[k]); end
close(ds)
sat_day = Dict(Date(t) => n for (n, t) in enumerate(sat_times))
lon_idx = [findall(abs.(λs .- λm(i)) .< 0.5 / CPD) for i in 1:Nx]
lat_idx = [findall(abs.(φs .- φm(j)) .< 0.5 / CPD) for j in 1:Ny]
function sat_cell(n, i, j)
    v = [x for x in chl_sat[lon_idx[i], lat_idx[j], n] if isfinite(x) && x > 0]
    return length(v) ≥ 2 ? Float64(mean(v)) : NaN
end

# ---------------- one pass over the limitation frames ----------------
ts_files = ["$(TAG)_volume_daily_rank$(r).jld2" for r in 0:NR-1]
ts_times = jldopen(f -> [f["timeseries/t/$it"] for it in iters(f, "T")], ts_files[1])
nearest(times, t) = (m = argmin(abs.(times .- t)); abs(times[m] - t) < 43200 ? m : 0)    # the de-tided fields start a few days after a restart

sums = Dict(name => zeros(Nx, Ny) for name in (:sat, :P, :cloern, :ratio, :I, :fN, :T))
cnt = zeros(Int, Nx, Ny)
for n in eachindex(lim_dates)
    d = Date(lim_dates[n])
    haskey(sat_day, d) || continue
    mP, mT = nearest(vol_times, lim_times[n]), nearest(ts_times, lim_times[n])
    (mP > 0 && mT > 0) || continue
    P = volume_frame("P", mP)[:, :, Nz]
    T, S = volume_frame("T", mT; stem = "volume_daily"), volume_frame("S", mT; stem = "volume_daily")
    PAR, nitrate, ammonium = limitation_frame("PAR", n), limitation_frame("nitrate", n), limitation_frame("ammonium", n)
    I = PAR_TO_PHOTONS .* mixed_layer_mean(PAR, mixed_layer_depth(T, S))              # mol photons/m²/day over the mixed layer
    fN = min.(1.0, nitrate[:, :, L] .+ ammonium[:, :, L])                              # = min(1, 2 L_N)
    chl_C = @. 0.003 + 0.0154 * exp(0.050 * T[:, :, Nz]) * exp(-0.059 * I) * fN
    chl_N = 12 * CARBON_RATIO .* chl_C                                                # mg Chl/mmol N
    s = [wet[i, j] ? sat_cell(sat_day[d], i, j) : NaN for i in 1:Nx, j in 1:Ny]
    ok = isfinite.(s) .& (P .> 0) .& isfinite.(chl_N)
    sums[:sat][ok] .+= s[ok]; sums[:P][ok] .+= P[ok]; sums[:cloern][ok] .+= (chl_N .* P)[ok]; sums[:ratio][ok] .+= chl_N[ok]
    sums[:I][ok] .+= I[ok]; sums[:fN][ok] .+= fN[ok]; sums[:T][ok] .+= T[:, :, Nz][ok]; cnt[ok] .+= 1
    say(@sprintf("  %s: %d matched cells, mean mixed-layer light %.1f mol/m²/day, Cloern Chl:N %.2f", d, count(ok), mean(I[ok]), mean(chl_N[ok])))
end

say("\nsurface chlorophyll by water depth (mg/m³, matched cell-days; ratios in mg Chl/mmol N)")
say("  bin           cell-days  satellite  fixed(1.31)  Cloern  | fixed/sat  Cloern/sat | required ratio  Cloern ratio | mixed-layer light  f_N")
for (lo, hi) in ((0, 20), (20, 40), (40, 60), (60, 80), (80, 100), (100, 200), (200, 1000), (1000, Inf))
    sel = wet .& (depth .> lo) .& (depth .≤ hi) .& (cnt .> 0)
    count(sel) == 0 && continue
    n = sum(cnt[sel]); S(name) = sum(sums[name][sel])
    sat, fixed, cloern = S(:sat) / n, CHL_PER_N * S(:P) / n, S(:cloern) / n
    say(@sprintf("  %5.0f-%-5.0f m %9d  %9.3f  %11.3f  %6.3f  | %9.2f  %10.2f | %14.2f  %12.2f | %17.1f  %.2f",
                 lo, hi, n, sat, fixed, cloern, fixed / sat, cloern / sat, S(:sat) / S(:P), S(:ratio) / n, S(:I) / n, S(:fN) / n))
end
say("\nbounds: Chl:C 0.005-0.06 g/g is $(round(12 * CARBON_RATIO * 0.005, digits = 2))-$(round(12 * CARBON_RATIO * 0.06, digits = 2)) mg Chl/mmol N (the model's fixed 1.31 is Chl:C = $(round(CHL_PER_N / (12 * CARBON_RATIO), digits = 4)) g/g)")
