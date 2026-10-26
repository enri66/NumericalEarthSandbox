# Growth limitation of LOBSTER's phytoplankton, mixed-layer depth and PAR, from a run of bgc/mab_bgc.jl with
# MAB_BGC_LIMITATION=true (bgc/bgc_limitation.jl). It reads the daily means of hourly samples in
# <TAG>_bgc_limitation_rank<r>.jld2 (PAR, the light, nitrate and ammonium terms, the growth rate in 1/day; the top
# MAB_BGC_LIMITATION_LEVELS levels) and the de-tided daily T and S in <TAG>_volume_daily_rank<r>.jld2, from which the
# mixed-layer depth is computed here (first depth at which the density exceeds that at 10 m by 0.03 kg/m³; the simplified
# equation of state of Roquet et al. 2015 at the surface pressure; mixed to the seafloor where it does not exceed it).
# The nutrient term is L_N = (L_NO₃ + L_NH₄) / 2 and the growth rate is μ₀ L_light L_N (the Q₁₀ factor is 1 in LOBSTER).
# Written: summary.txt and figures (timeseries, maps, profiles) in BGC_OUT.
#
#   julia -t 4 --project=<env with CairoMakie> bgc/limitation_analysis.jl
# Environment: BGC_TAG (run path prefix), MAB_START_DATE, MAB_CELLS_PER_DEGREE, MAB_NZ, BGC_OUT.
using Printf, Statistics, Dates
using CairoMakie
using Oceananigans
using Oceananigans.Grids: znodes
using NumericalEarth.DataWrangling: jldopen

const TAG   = get(ENV, "BGC_TAG", "/t0/workdir/enrique/runs/bgc/bgc_feb19_woa_lim")
const START = DateTime(get(ENV, "MAB_START_DATE", "2019-02-01"))
const CPD   = parse(Int, get(ENV, "MAB_CELLS_PER_DEGREE", "12"))
const OUT   = get(ENV, "BGC_OUT", "/t0/workdir/enrique/runs/analysis/limitation_" * basename(TAG))
const H, Nz = 7, parse(Int, get(ENV, "MAB_NZ", "50"))
const λ_west, φ_south = -76 + 2 / 12, 34 + 2 / 12        # script 05's box: GLORYS box edge plus two GLORYS cells
const MLD_THRESHOLD, MLD_REFERENCE = 0.03, 10.0
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
regions = (shelf = (depth .> 0) .& (depth .< 100), slope = (depth .≥ 100) .& (depth .< 1000), deep = depth .≥ 1000)
say("model: $(Nx)×$(Ny)×$(Nz) on $NR ranks, $(count(wet)) wet columns; shelf $(count(regions.shelf)), slope $(count(regions.slope)), deep $(count(regions.deep)) columns")

# ---------------- the limitation output ----------------
iters(f, name) = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/$name"]))))
limitation_files = ["$(TAG)_bgc_limitation_rank$(r).jld2" for r in 0:NR-1]
# runs from before 2026-10-26 lack the concentrations and the daylight fraction
const FIELDS = Tuple(name for name in (:PAR, :light, :nitrate, :ammonium, :temperature, :growth, :NO₃, :NH₄, :daylight)
                     if jldopen(f -> haskey(f["timeseries"], String(name)), limitation_files[1]))
const HAS_DAYLIGHT = :daylight in FIELDS
lim_times = jldopen(f -> [f["timeseries/t/$it"] for it in iters(f, "PAR")], limitation_files[1])
Nd = length(lim_times)
const L = jldopen(f -> size(f["timeseries/PAR/$(iters(f, "PAR")[1])"], 3), limitation_files[1])
const klev = Nz-L+1:Nz                                   # model levels of the saved ones, bottom first
lim = Dict(name => fill(NaN32, Nx, Ny, L, Nd) for name in FIELDS)
for (p, file) in zip(parts, limitation_files), name in FIELDS
    jldopen(file) do f
        for (n, it) in enumerate(iters(f, String(name)))
            B = f["timeseries/$name/$it"]
            nx, ny = size(p.kb)
            size(B, 1) == nx + 2H && (B = B[H+1:end-H, H+1:end-H, :])
            lim[name][p.I+1:p.I+nx, p.J+1:p.J+ny, :, n] .= B
        end
    end
end
dates = [START + Millisecond(round(Int, 1000t)) for t in lim_times]
say("limitation frames: $Nd daily means, $(Dates.format(first(dates), "yyyy-mm-dd")) → $(Dates.format(last(dates), "yyyy-mm-dd")); top $L levels, down to $(round(-zf[Nz-L+1], digits = 1)) m")
const depth_c = -zc[klev]                                # level-centre depths (m), bottom first within the saved levels
upper = -zf[klev .+ 1]; lower = -zf[klev]                # upper and lower face depth of each saved level

# ---------------- mixed-layer depth from the daily T and S ----------------
# simplified equation of state, Roquet et al. 2015: density anomaly (kg/m³) of potential temperature Θ (°C) and salinity S at the surface
function density_anomaly(Θ, S)
    Θ′, S′ = Θ - 10, S - 35
    return -0.1655 * (1 + 0.5 * 0.05952 * Θ′) * Θ′ + 0.76554 * (1 - 0.5 * 0.00054914 * S′) * S′ - 0.0024341 * Θ′ * S′
end
volume_files = ["$(TAG)_volume_daily_rank$(r).jld2" for r in 0:NR-1]
vol_times = jldopen(f -> [f["timeseries/t/$it"] for it in iters(f, "T")], volume_files[1])
function volume_frame(name, n)
    A = zeros(Float32, Nx, Ny, Nz)
    for (p, file) in zip(parts, volume_files)
        jldopen(file) do f
            it = iters(f, name)[n]
            B = f["timeseries/$name/$it"][H+1:end-H, H+1:end-H, H+1:end-H]
            A[p.I+1:p.I+size(B, 1), p.J+1:p.J+size(B, 2), :] .= B
        end
    end
    return A
end
mld = fill(NaN, Nx, Ny, Nd); to_bottom = falses(Nx, Ny, Nd)
function mixed_layer_depth!(mld, to_bottom, n, T, S)
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
                mld[i, j, n] = -zc[k+1] + w * (zc[k+1] - zc[k]); crossed = true; break
            end
        end
        crossed || (mld[i, j, n] = depth[i, j]; to_bottom[i, j, n] = true)
    end
end
for n in 1:Nd
    isempty(vol_times) && break
    m = argmin(abs.(vol_times .- lim_times[n]))
    abs(vol_times[m] - lim_times[n]) < 43200 || continue          # the de-tided T, S start a few days after the restart
    mixed_layer_depth!(mld, to_bottom, n, volume_frame("T", m), volume_frame("S", m))
end
say("mixed-layer depth computed on $(count(n -> any(!isnan, mld[:, :, n]), 1:Nd)) of the $Nd days")

# ---------------- derived quantities ----------------
L_N = (lim[:nitrate] .+ lim[:ammonium]) ./ 2                      # the nutrient term
surface(name) = lim[name][:, :, L, :]
# the light term of the daytime: the daily mean of the term includes the night (the fraction of the samples with PAR > 1 W/m² is saved)
daytime_light = HAS_DAYLIGHT ? ifelse.(lim[:daylight] .> 0.02f0, lim[:light] ./ max.(lim[:daylight], 0.02f0), NaN32) : lim[:light]
HAS_DAYLIGHT || say("no daylight fraction in this run: the light term below is the daily mean, which includes the night")
PAR0 = surface(:PAR)
# euphotic depth (1% of the surface PAR) from the saved levels; NaN where PAR stays above 1% to the lowest saved level
zeu = fill(NaN, Nx, Ny, Nd)
for n in 1:Nd, j in 1:Ny, i in 1:Nx
    (wet[i, j] && isfinite(PAR0[i, j, n]) && PAR0[i, j, n] > 0) || continue
    for l in L-1:-1:1
        if lim[:PAR][i, j, l, n] < 0.01 * PAR0[i, j, n]
            zeu[i, j, n] = depth_c[l+1] + (depth_c[l] - depth_c[l+1]) * ((0.01 * PAR0[i, j, n]) - lim[:PAR][i, j, l+1, n]) / (lim[:PAR][i, j, l, n] - lim[:PAR][i, j, l+1, n])
            break
        end
    end
end
# mean of a saved quantity over the mixed layer, weighted by the thickness of the levels (over the saved levels only where the
# mixed layer is deeper than them)
function mixed_layer_mean(A)
    out = fill(NaN, Nx, Ny, Nd)
    for n in 1:Nd, j in 1:Ny, i in 1:Nx
        (isfinite(mld[i, j, n]) && kb[i, j] > 0) || continue
        s = 0.0; w = 0.0
        for l in 1:L
            Δ = clamp(mld[i, j, n] - upper[l], 0.0, lower[l] - upper[l])
            isfinite(A[i, j, l, n]) && (s += Δ * A[i, j, l, n]; w += Δ)
        end
        w > 0 && (out[i, j, n] = s / w)
    end
    return out
end
region_values(A, mask) = reshape(A, Nx * Ny, size(A, 3))[vec(mask), :]

# ---------------- moorings ----------------
const MOORINGS = ("CNSM" => (-70.88, 40.10, 147.0), "ISSM" => (-70.88, 40.36, 92.0))        # lon, lat, nominal water depth (m)
cells = Dict{String, Tuple{Int, Int}}()
for (name, (lon, lat, nominal)) in MOORINGS
    cand = [(i, j) for i in 1:Nx, j in 1:Ny if wet[i, j] && abs(λm(i) - lon) * cosd(lat) < 0.13 && abs(φm(j) - lat) < 0.13]
    best = argmin([abs(depth[i, j] - nominal) / nominal + 0.5 * hypot((λm(i) - lon) * cosd(lat), φm(j) - lat) / 0.13 for (i, j) in cand])
    cells[name] = cand[best]
    say(@sprintf("%s: %.3f°E %.3f°N, nominal %.0f m → cell (%d, %d) at %.3f°E %.3f°N, model depth %.0f m", name, lon, lat, nominal, cand[best]..., λm(cand[best][1]), φm(cand[best][2]), depth[cand[best]...]))
end

# series of a field for a mooring cell or a region mean (over wet cells with finite values on that day)
point_series(A, (i, j)) = [A[i, j, n] for n in 1:Nd]
region_series(A, mask) = [(v = A[:, :, n][mask .& isfinite.(A[:, :, n])]; isempty(v) ? NaN : mean(v)) for n in 1:Nd]
series_sets = (CNSM = A -> point_series(A, cells["CNSM"]), ISSM = A -> point_series(A, cells["ISSM"]),
               shelf = A -> region_series(A, regions.shelf), slope = A -> region_series(A, regions.slope), deep = A -> region_series(A, regions.deep))
ML_light = mixed_layer_mean(lim[:light]); ML_growth = mixed_layer_mean(lim[:growth]); ML_N = mixed_layer_mean(L_N)
quantities = (; PAR = PAR0, light = surface(:light), daytime = daytime_light[:, :, L, :],
                NO₃ = HAS_DAYLIGHT ? surface(:NO₃) : fill(NaN32, Nx, Ny, Nd), NH₄ = HAS_DAYLIGHT ? surface(:NH₄) : fill(NaN32, Nx, Ny, Nd), nitrate = surface(:nitrate), ammonium = surface(:ammonium), nutrient = L_N[:, :, L, :],
                growth = surface(:growth), MLD = mld, zeu = zeu, ML_light, ML_nutrient = ML_N, ML_growth)

say("\nmeans over the $Nd days at the surface level (top cell, ~$(round(depth_c[L], digits = 1)) m): PAR W/m², terms 0-1, growth 1/day; mixed-layer depth and euphotic depth in m")
say(@sprintf("%-6s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s", "", "PAR", "L_light", "L_NO₃", "L_NH₄", "L_N", "μ/P", "MLD", "z_eu", "ML:Ll", "ML:LN", "ML:μ/P"))
means = Dict{Symbol, Dict{Symbol, Float64}}()
for (key, f) in pairs(series_sets)
    s = Dict(q => f(A) for (q, A) in pairs(quantities))
    means[key] = Dict(q => (v = filter(isfinite, x); isempty(v) ? NaN : mean(v)) for (q, x) in s)
    m = means[key]
    say(@sprintf("%-6s %7.1f %7.3f %7.3f %7.3f %7.3f %7.3f %7.1f %7.1f %7.3f %7.3f %7.3f", String(key), m[:PAR], m[:light], m[:nitrate], m[:ammonium], m[:nutrient], m[:growth], m[:MLD], m[:zeu], m[:ML_light], m[:ML_nutrient], m[:ML_growth]))
end

months = unique(Dates.yearmonth.(dates))
say("\nmonthly means of the surface level (daytime light term is the mean of the term over the hours with PAR > 1 W/m²; NO₃, NH₄ in mmol/m³)")
say(@sprintf("%-8s %-8s %5s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s", "month", "", "days", "PAR", "L_light", "L_day", "L_N", "μ/P", "MLD", "z_eu", "NO₃", "NH₄", "LD<LN %"))
for (y, m) in months, (key, f) in pairs(series_sets)
    sel = findall(d -> Dates.yearmonth(d) == (y, m), dates)
    mean_over(q) = (v = filter(isfinite, f(quantities[q])[sel]); isempty(v) ? NaN : mean(v))
    a = f(quantities.daytime)[sel]; b = f(quantities.nutrient)[sel]; ok = isfinite.(a) .& isfinite.(b)
    say(@sprintf("%04d-%02d  %-8s %5d %7.1f %7.3f %7.3f %7.3f %7.3f %7.1f %7.1f %7.3f %7.4f %7.0f", y, m, String(key), length(sel), mean_over(:PAR), mean_over(:light), mean_over(:daytime),
                 mean_over(:nutrient), mean_over(:growth), mean_over(:MLD), mean_over(:zeu), mean_over(:NO₃), mean_over(:NH₄), any(ok) ? 100 * mean(a[ok] .< b[ok]) : NaN))
end

say("\nwhich factor limits growth at the surface: share of cell-days with the daytime L_light < L_N (light-limited) and L_N < L_light (nutrient-limited), and of the nutrient term from nitrate (L_NO₃/2 over L_N)")
for (key, mask) in (("shelf", regions.shelf), ("slope", regions.slope), ("deep", regions.deep))
    a = region_values(quantities.daytime, mask); b = region_values(quantities.nutrient, mask); c = region_values(quantities.nitrate, mask)
    ok = isfinite.(a) .& isfinite.(b) .& (b .> 0)
    say(@sprintf("  %-6s light-limited %5.1f %%   nutrient-limited %5.1f %%   nitrate share of L_N %5.1f %%   L_N < 0.1: %5.1f %%   L_light < 0.1: %5.1f %%",
                 key, 100 * mean(a[ok] .< b[ok]), 100 * mean(b[ok] .< a[ok]), 100 * mean(c[ok] ./ 2 ./ b[ok]), 100 * mean(b[ok] .< 0.1), 100 * mean(a[ok] .< 0.1)))
end
say("\nfraction of columns with the mixed layer reaching the seafloor, and mean mixed-layer depth / water depth on the shelf")
let n_ok = [count(isfinite, mld[:, :, n][regions.shelf]) for n in 1:Nd]
    tb = [count(to_bottom[:, :, n][regions.shelf]) for n in 1:Nd]
    say(@sprintf("  shelf: %.1f %% of column-days mixed to the seafloor (over %d days with MLD)", 100 * sum(tb) / max(1, sum(n_ok)), count(>(0), n_ok)))
end

# ---------------- figures ----------------
days = [Dates.value(d - first(dates)) / 86400000 for d in dates]
colors = (light = :goldenrod, nitrate = :royalblue, ammonium = :seagreen, nutrient = :black, growth = :firebrick)
fig = Figure(size = (1700, 800))
for (c, key) in enumerate(keys(series_sets))
    f = series_sets[key]
    ax = Axis(fig[1, c], title = String(key) * " — surface terms", xlabel = "days from $(Date(first(dates)))", ylabel = "limitation (0-1)", limits = (nothing, (0, 1)))
    for q in (:light, :nitrate, :ammonium, :nutrient); lines!(ax, days, f(quantities[q]); color = colors[q], label = String(q), linewidth = q == :nutrient ? 3 : 2); end
    ax2 = Axis(fig[1, c], yaxisposition = :right, ylabel = "growth (1/day)", limits = (nothing, (0, 0.6)))
    hidespines!(ax2); hidexdecorations!(ax2)
    lines!(ax2, days, f(quantities.growth); color = colors.growth, linestyle = :dash, linewidth = 2, label = "μ/P")
    c == 1 && axislegend(ax, position = :lt, nbanks = 2)
    ax3 = Axis(fig[2, c], title = String(key) * " — mixed layer and PAR", xlabel = "days from $(Date(first(dates)))", ylabel = "depth (m)", yreversed = true)
    lines!(ax3, days, f(quantities.MLD); color = :black, linewidth = 2, label = "MLD")
    lines!(ax3, days, f(quantities.zeu); color = :orange, linewidth = 2, label = "1% PAR depth")
    ax4 = Axis(fig[2, c], yaxisposition = :right, ylabel = "surface PAR (W/m²)")
    hidespines!(ax4); hidexdecorations!(ax4)
    lines!(ax4, days, f(quantities.PAR); color = :grey, linewidth = 1.5, label = "PAR")
    c == 1 && axislegend(ax3, position = :rb)
end
save(joinpath(OUT, "limitation_timeseries.png"), fig)

land = .!wet
function mean_map(A)
    M = fill(NaN, Nx, Ny)
    for j in 1:Ny, i in 1:Nx
        v = filter(isfinite, view(A, i, j, :)); wet[i, j] && !isempty(v) && (M[i, j] = mean(v))
    end
    return M
end
figm = Figure(size = (1700, 1000))
for (n, (q, label, range)) in enumerate(((:PAR, "surface PAR (W/m²)", (0, 80)), (:light, "L_light", (0, 0.4)), (:nitrate, "L_NO₃", (0, 1)), (:ammonium, "L_NH₄", (0, 1)),
                                         (:nutrient, "L_N", (0, 1)), (:growth, "growth μ/P (1/day)", (0, 0.6)), (:MLD, "mixed-layer depth (m)", (0, 150)), (:zeu, "1% PAR depth (m)", (0, 40))))
    ax = Axis(figm[fld1(n, 4), mod1(n, 4)], title = label * ", mean over the run", xlabel = "lon", ylabel = "lat", aspect = DataAspect())
    hm = heatmap!(ax, [λm(i) for i in 1:Nx], [φm(j) for j in 1:Ny], mean_map(quantities[q]); colorrange = range, colormap = :viridis)
    Colorbar(figm[fld1(n, 4), mod1(n, 4)][1, 2], hm)
    for (name, (i, j)) in cells; scatter!(ax, [λm(i)], [φm(j)]; color = :red, marker = :star5, markersize = 12); end
end
save(joinpath(OUT, "limitation_maps.png"), figm)

figp = Figure(size = (1500, 500))
for (c, key) in enumerate((:CNSM, :ISSM, :shelf, :slope))
    ax = Axis(figp[1, c], title = String(key) * " — mean profile, top $L levels", xlabel = "term (0-1) / growth (1/day)", ylabel = "depth (m)", yreversed = true, limits = ((0, 1), nothing))
    prof(name) = key in (:CNSM, :ISSM) ? [mean(filter(isfinite, lim[name][cells[String(key)]..., l, :])) for l in 1:L] :
                 [mean(filter(isfinite, region_values(lim[name][:, :, l, :], key == :shelf ? regions.shelf : regions.slope))) for l in 1:L]
    for name in (:light, :nitrate, :ammonium); lines!(ax, prof(name), depth_c; color = colors[name], label = String(name), linewidth = 2); end
    lines!(ax, prof(:growth), depth_c; color = colors.growth, linestyle = :dash, linewidth = 2, label = "μ/P (1/day)")
    c == 1 && axislegend(ax, position = :rb)
end
save(joinpath(OUT, "limitation_profiles.png"), figp)
say("\nfigures and summary in $OUT")
