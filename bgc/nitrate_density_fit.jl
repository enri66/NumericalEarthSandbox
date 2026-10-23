# Nitrate as a function of potential density, fitted to WOA 2023 monthly climatology over the MAB box, and tested at the Pioneer moorings
# against the SUNA nitrate sensor, with the density taken from (a) GLORYS temperature and salinity at the mooring and (b) the
# temperature and salinity of the mooring's own oxygen optode (cp0xxxm_dosta_nsif.csv).
#
# σ = ρ(Θ, Sᴬ, 0) − 1000 (TEOS-10 at the surface, Θ ≈ the potential or in-situ temperature, Sᴬ = 1.004715 S): the predictor.
# The relation for each month (Jan-Apr) is the median WOA nitrate (µmol/kg × 1.025 → mmol/m³) in σ bins of 0.1 kg/m³ over the cells in the
# box (-76..-64°E, 34..42°N) shallower than 800 m (FIT_MAXDEPTH restricts it to WOA columns whose data end shallower than that), with the gaps
# interpolated and the result made non-decreasing in σ. It is written to nitrate_density_table.csv. Result (2026-10-23): the relation does not
# reproduce SUNA at CNSM and has no use for the model, see HANDOFF_BGC.md; no model option reads the table.
#   julia --project=<env of mpi05_bgc> bgc/nitrate_density_fit.jl [output directory for the csv]
using NumericalEarth, Oceananigans, Printf, Statistics, Dates
using Oceananigans.Grids: λnodes, φnodes, znodes
using NumericalEarth.DataWrangling: Metadatum, BoundingBox
using NumericalEarth.WOA: WOAMonthly
using WorldOceanAtlasTools
ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"
const NCD = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "NCDatasets")]
const SP  = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "SeawaterPolynomials")]
const eos = SP.TEOS10.TEOS10EquationOfState()
sigma(T, S) = SP.ρ(Float64(T), 1.004715 * Float64(S), 0.0, eos) - 1000

const OUTDIR = length(ARGS) ≥ 1 ? ARGS[1] : "."
const DATA = get(ENV, "DATA_DIR", "/t0/workdir/enrique/Data")
const GLORYS_DIR = joinpath(DATA, "NumericalEarth"); const OOI = joinpath(DATA, "OOI", "pioneer_bgc_2019")
const BOX = (λ = (-76.0, -64.0), φ = (34.0, 42.0), zmin = -800.0)
const MAXD = parse(Float64, get(ENV, "FIT_MAXDEPTH", "Inf"))      # only the WOA columns whose temperature data end shallower than this (shelf and upper slope)
const EDGES = 24.0:0.1:28.0; const CENTERS = collect(EDGES[1:end-1] .+ step(EDGES) / 2)
const WOA_REGION = BoundingBox(longitude = (-78, -62), latitude = (32, 44))
woa(var, m) = Field(Metadatum(var; dataset = WOAMonthly(), region = WOA_REGION, date = DateTime(2018, m, 1)), CPU(); inpainting = nothing, cache_inpainted_data = false)

# ---------------- the fit ----------------
function woa_samples(m)
    T, S, N = woa(:temperature, m), woa(:salinity, m), woa(:nitrate, m)
    λ, φ, z = collect(λnodes(N.grid, Center())), collect(φnodes(N.grid, Center())), collect(znodes(N.grid, Center()))
    At, As, An = Array(interior(T)), Array(interior(S)), Array(interior(N))
    column_bottom(i, j) = (k = findfirst(k -> isfinite(At[i, j, k]), eachindex(z)); isnothing(k) ? 0.0 : -z[k])      # z is bottom first
    σ = Float64[]; n = Float64[]
    for k in eachindex(z), j in eachindex(φ), i in eachindex(λ)
        (BOX.λ[1] < λ[i] < BOX.λ[2] && BOX.φ[1] < φ[j] < BOX.φ[2] && z[k] > BOX.zmin && column_bottom(i, j) ≤ MAXD) || continue
        t, s, x = At[i, j, k], As[i, j, k], An[i, j, k]
        (isfinite(t) && isfinite(s) && isfinite(x) && x ≥ 0) || continue
        push!(σ, sigma(t, s)); push!(n, 1.025 * x)
    end
    return σ, n
end
function lookup_table(σ, n)
    med = [(sel = (σ .≥ EDGES[b]) .& (σ .< EDGES[b+1]); count(sel) ≥ 4 ? median(n[sel]) : NaN) for b in eachindex(CENTERS)]
    ok = findall(isfinite, med)
    filled = [interp1(CENTERS[ok], med[ok], c) for c in CENTERS]
    return accumulate(max, filled)                       # non-decreasing in σ
end
function interp1(x, y, xq)
    xq ≤ x[1] && return y[1]; xq ≥ x[end] && return y[end]
    k = searchsortedlast(x, xq); w = (xq - x[k]) / (x[k+1] - x[k]); (1 - w) * y[k] + w * y[k+1]
end
nitrate_of(table, σ) = interp1(CENTERS, table, σ)

samples = Dict(m => woa_samples(m) for m in 1:4)
tables = Dict(m => lookup_table(samples[m]...) for m in 1:4)
println("fit restricted to WOA columns shallower than $MAXD m; WOA samples per month: ", join(["$m: $(length(samples[m][1]))" for m in 1:4], ", "))
println("\nnitrate (mmol/m³) as a function of σ: month tables at selected σ")
println("   σ      ", join(["Jan   ", "Feb   ", "Mar   ", "Apr   "], " "), "   WOA scatter in the Feb bin (p25–p75)")
for σ in 24.5:0.5:27.5
    b = clamp(searchsortedlast(EDGES, σ), 1, length(CENTERS))
    s = samples[2][1]; sel = (s .≥ EDGES[b]) .& (s .< EDGES[b+1]); q = count(sel) ≥ 4 ? quantile(samples[2][2][sel], [0.25, 0.75]) : [NaN, NaN]
    @printf("   %.1f   %s    n=%d  %.1f–%.1f\n", σ, join([@sprintf("%5.1f ", nitrate_of(tables[m], σ)) for m in 1:4], " "), count(sel), q...)
end
println("\nin-sample skill on the WOA cells (relation vs WOA nitrate): ", join([let (σ, n) = samples[m]; p = [nitrate_of(tables[m], x) for x in σ]; @sprintf("month %d: rms %.2f, r² %.2f", m, sqrt(mean((p .- n) .^ 2)), cor(p, n)^2) end for m in 1:4], "; "))
open(joinpath(OUTDIR, "nitrate_density_table.csv"), "w") do io
    println(io, "# nitrate (mmol/m3) vs sigma = rho(Theta, 1.004715 S, 0 dbar) - 1000 (kg/m3), WOA 2023 monthly medians over -76..-64E 34..42N shallower than 800 m, non-decreasing")
    println(io, "sigma,jan,feb,mar,apr")
    for (b, c) in enumerate(CENTERS); println(io, @sprintf("%.2f,%.3f,%.3f,%.3f,%.3f", c, (tables[m][b] for m in 1:4)...)); end
end

# ---------------- the test at the moorings ----------------
table_at(date) = (f = (Dates.day(date) - 1 + Dates.hour(date) / 24) / Dates.daysinmonth(date); m = Dates.month(date); w = f < 0.5 ? 0.5 - f : f - 0.5; m₂ = clamp(f < 0.5 ? m - 1 : m + 1, 1, 4);
                  (σ -> (1 - w) * nitrate_of(tables[clamp(m, 1, 4)], σ) + w * nitrate_of(tables[m₂], σ)))
function read_csv(path, columns)
    lines = readlines(path); col = Dict(h => i for (i, h) in enumerate(split(lines[1], ",")))
    out = Vector{Vector{Float64}}(); t = DateTime[]
    for line in lines[3:end]
        f = split(line, ","); vals = [tryparse(Float64, f[col[c]]) for c in columns]; any(isnothing, vals) && continue
        push!(t, DateTime(f[1][1:19])); push!(out, Float64.(vals))
    end
    lat = parse(Float64, split(lines[3], ",")[col["latitude"]]); lon = parse(Float64, split(lines[3], ",")[col["longitude"]])
    return t, out, lon, lat
end
window_mean(t, v, T; half = Hour(36), n = 24) = (s = [x for (τ, x) in zip(t, v) if abs(τ - T) ≤ half]; length(s) ≥ n ? mean(s) : NaN)
function glorys_at(lon, lat, date, depth)
    day = Dates.format(date, "yyyy-mm-dd")
    vals = map(("thetao", "so")) do v
        ds = NCD.NCDataset(joinpath(GLORYS_DIR, "$(v)_GLORYSDaily_$(day)T00-00-00_$(day)T00-00-00_-76.0_-64.0_34.0_42.0.nc"))
        lons, lats, zs = Float64.(ds["longitude"][:]), Float64.(ds["latitude"][:]), Float64.(ds["depth"][:])
        A = ds[v][:, :, :, 1]; best = (Inf, 1, 1)
        for i in eachindex(lons), j in eachindex(lats); d = hypot((lons[i] - lon) * cosd(lat), lats[j] - lat); (d < best[1] && !ismissing(A[i, j, 1])) && (best = (d, i, j)); end
        col = [ismissing(x) ? NaN : Float64(x) for x in A[best[2], best[3], :]]; close(ds)
        k = searchsortedlast(zs, depth); w = (depth - zs[k]) / (zs[k+1] - zs[k]); (1 - w) * col[k] + w * col[k+1]
    end
    return vals[1], vals[2]
end
days = [DateTime(2019, 2, 1) + Day(d) for d in 0:58]
for (name, prefix) in ("CNSM" => "cp01cnsm", "ISSM" => "cp03issm")
    tn, vn, lon, lat = read_csv(joinpath(OOI, "$(prefix)_nutnr_all.csv"), ["mole_concentration_of_nitrate_in_sea_water_suna"])
    keep = [0 ≤ x[1] ≤ 40 for x in vn]; tn, nit = tn[keep], [x[1] for x in vn[keep]]
    td, vd, _, _ = read_csv(joinpath(OOI, "$(prefix)_dosta_nsif.csv"), ["sea_water_temperature", "sea_water_practical_salinity"])
    ok = [0 < x[1] < 30 && 25 < x[2] < 38 for x in vd]; td, vd = td[ok], vd[ok]
    suna = [window_mean(tn, nit, T) for T in days]
    pred_dosta = [(f = table_at(T); window_mean(td, [f(sigma(x[1], x[2])) for x in vd], T)) for T in days]
    sig_dosta = [window_mean(td, [sigma(x[1], x[2]) for x in vd], T) for T in days]
    pred_glorys = Float64[]; sig_glorys = Float64[]
    for T in days
        θ, S = glorys_at(lon, lat, T, 7.0); push!(sig_glorys, sigma(θ, S)); push!(pred_glorys, table_at(T)(sig_glorys[end]))
    end
    println("\n=== $name ($(round(lon, digits = 3))°E $(round(lat, digits = 3))°N), 7 m, days 2019-02-01 → 03-31")
    for (label, rng) in (("all", 1:59), ("Feb", 1:28), ("Mar", 29:59))
        s, d, g = suna[rng], pred_dosta[rng], pred_glorys[rng]
        @printf("   %-3s mean SUNA %5.2f | from the optode's T,S: %5.2f | from GLORYS T,S: %5.2f\n", label, mean(filter(isfinite, s)), mean(filter(isfinite, d)), mean(g))
    end
    for (label, p) in (("from the optode's T,S", pred_dosta), ("from GLORYS T,S", pred_glorys))
        ok = isfinite.(suna) .& isfinite.(p)
        @printf("   %-22s vs SUNA over %d days: bias %+.2f, rms %.2f, r %+.2f, std ratio %.2f\n", label, count(ok), mean(p[ok] .- suna[ok]), sqrt(mean((p[ok] .- suna[ok]) .^ 2)), cor(p[ok], suna[ok]), std(p[ok]) / std(suna[ok]))
    end
    @printf("   σ (kg/m³): SUNA site (optode) %.2f ± %.2f; GLORYS %.2f ± %.2f; correlation of the two σ series %.2f\n", mean(filter(isfinite, sig_dosta)), std(filter(isfinite, sig_dosta)), mean(sig_glorys), std(sig_glorys),
            cor(sig_glorys[isfinite.(sig_dosta)], sig_dosta[isfinite.(sig_dosta)]))
end
