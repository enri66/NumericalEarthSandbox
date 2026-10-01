# The surface fluxes of a script-04 run (bulk formulas applied to ERA5 and the model's own SST and currents) against
# ERA5's own surface fluxes over the same days: wind stress, sensible and latent heat. Daily means; ERA5 is
# interpolated to the model's ocean cells from its ocean points only (land-sea mask 0), so coastal values are not
# contaminated by land. Needs the run's `<tag>_fluxes_daily.jld2` (MAB_FLUX_OUTPUT, the default) and ERA5's
# mean_surface_* fields (scripts/download_era5_year.jl with ERA5_NAMES). Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/flux_test/fluxes julia --project=. scripts/fluxes_vs_era5.jl
#
# Sign conventions are reported, not assumed: each model flux is compared with ERA5 as is and with its sign reversed,
# and the table says which matches.
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const TAG    = get(ENV, "MAB_TAG", "fluxes")
const PREFIX = run_prefix(TAG)
const ERA5_DIR = joinpath(GDIR, "era5")
const BOX = "_-76.0_-64.0_34.0_42.0.nc"

# ---------------- model ----------------
path = PREFIX * "_fluxes_daily.jld2"
names = (:τx, :τy, :sensible_heat, :latent_heat)
M = Dict(n => FieldTimeSeries(path, String(n)) for n in (names..., :friction_velocity))
grid = M[:τx].grid; ug = grid.underlying_grid
Nx, Ny = ug.Nx, ug.Ny
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center()))
B = grid.immersed_boundary.bottom_height
wet = [B[i, j, 1] < 0 for i in 1:Nx, j in 1:Ny]
centred(A) = A[1:Nx, 1:Ny]   # all these fields are cell-centred on the exchange grid
times = M[:τx].times
days = [floor(Int, t / 86400 - 0.5) for t in times]   # an AveragedTimeInterval frame at t averages the day before t

# ---------------- ERA5 ----------------
era5_names = Dict(:τx => ("mean_eastward_turbulent_surface_stress", "avg_iews"),
                  :τy => ("mean_northward_turbulent_surface_stress", "avg_inss"),
                  :sensible_heat => ("mean_surface_sensible_heat_flux", "avg_ishf"),
                  :latent_heat => ("mean_surface_latent_heat_flux", "avg_slhtf"))
era5_file(cds, date) = joinpath(ERA5_DIR, "$(cds)_ERA5HourlySingleLevel_$(Dates.format(date, "yyyy-mm-ddTHH"))$BOX")

lsm_ds = NCD.Dataset(joinpath(ERA5_DIR, "land_sea_mask_ERA5HourlySingleLevel_2019-01-01T00$BOX"))
elon = Float64.(lsm_ds["longitude"][:]); elat = Float64.(lsm_ds["latitude"][:])
eland = Float64.(lsm_ds["lsm"][:, :, 1]) .> 0
close(lsm_ds)
# ERA5 latitudes are stored north to south; bilinear() wants increasing coordinates
order = sortperm(elat); elat = elat[order]; eland = eland[:, order]

function era5_daily(name, day)
    cds, var = era5_names[name]
    acc = zeros(length(elon), length(elat)); n = 0
    for h in 1:24   # ERA5's mean rates are stamped at the end of their hour
        f = era5_file(cds, start_date + Day(day) + Hour(h))
        isfile(f) || continue
        ds = NCD.Dataset(f); a = Float64.(ds[var][:, :, 1]); close(ds)
        acc .+= a[:, order]; n += 1
    end
    n == 24 || return nothing
    A = acc ./ n
    A[eland] .= NaN
    return [wet[i, j] ? bilinear(A, elon, elat, λ[i], φ[j]) : NaN for i in 1:Nx, j in 1:Ny]
end

# ---------------- comparison ----------------
stats = Dict{Symbol, Vector{NTuple{4, Float64}}}()   # per day: model mean, ERA5 mean, rms difference, correlation
mean_maps = Dict{Symbol, Any}()
for name in names
    rows = NTuple{4, Float64}[]; msum = zeros(Nx, Ny); esum = zeros(Nx, Ny); nd = 0
    for (n, d) in enumerate(days)
        d < 0 && continue
        e = era5_daily(name, d); isnothing(e) && continue
        m = centred(Array(interior(M[name][n]))[:, :, 1])
        ok = wet .& isfinite.(e)
        x = m[ok]; y = e[ok]
        push!(rows, (mean(x), mean(y), sqrt(mean((x .- y) .^ 2)), cor(x, y)))
        msum .+= ifelse.(ok, m, 0.0); esum .+= ifelse.(ok, e, 0.0); nd += 1
    end
    stats[name] = rows
    mean_maps[name] = (msum ./ nd, esum ./ nd)
end

println("Daily-mean domain averages over ocean cells, $(basename(PREFIX)) vs ERA5 (days with complete ERA5):")
for name in names
    r = stats[name]
    m = mean(first.(r)); e = mean(getindex.(r, 2)); c = mean(last.(r))
    sign_note = c < 0 ? "  (opposite sign convention: model = −ERA5)" : ""
    @printf("  %-14s model %+9.4f  ERA5 %+9.4f  |model|/|ERA5| %.2f  daily rms diff %.4f  mean daily pattern corr %+.2f%s\n",
            name, m, e, abs(m) / abs(e), mean(getindex.(r, 3)), c, sign_note)
end

# wind-stress magnitude, which does not depend on the sign convention
mx, ex = mean_maps[:τx]; my, ey = mean_maps[:τy]
mag_m = sqrt.(mx .^ 2 .+ my .^ 2); mag_e = sqrt.(ex .^ 2 .+ ey .^ 2)
ok = wet .& isfinite.(mag_e) .& (mag_e .> 0)
@printf("  |τ| of the period-mean stress: model %.4f, ERA5 %.4f N/m² (ratio %.2f)\n", mean(mag_m[ok]), mean(mag_e[ok]), mean(mag_m[ok]) / mean(mag_e[ok]))
u★ = mean(mean(centred(Array(interior(M[:friction_velocity][n]))[:, :, 1])[wet]) for n in eachindex(days) if days[n] >= 0)
@printf("  model friction velocity, domain and period mean: %.4f m/s\n", u★)

println("\nDaily domain means (model / ERA5):")
@printf("%5s %22s %22s %22s %22s\n", "day", "τx", "τy", "sensible", "latent")
for k in eachindex(stats[:τx])
    @printf("%5d", k)
    for name in names
        @printf("   %+8.3f / %+8.3f", stats[name][k][1], stats[name][k][2])
    end
    println()
end

# ---------------- maps of the period means ----------------
fig = Figure(size = (1800, 1250), fontsize = 15)
Label(fig[0, 1:3], "$(basename(PREFIX)): period-mean surface fluxes, model (bulk formulas) vs ERA5", fontsize = 20)
mapped = ((:latent_heat, "latent heat (W/m²)"), (:sensible_heat, "sensible heat (W/m²)"))
mask(A) = [wet[i, j] ? A[i, j] : NaN for i in 1:Nx, j in 1:Ny]
for (row, (name, label)) in enumerate(mapped)
    m, e = mean_maps[name]
    s = mean(last.(stats[name])) < 0 ? -1 : 1      # put the model on ERA5's sign convention
    cr = extrema(filter(isfinite, mask(e)))
    for (col, (A, t)) in enumerate(((s .* m, "model"), (e, "ERA5"), (s .* m .- e, "model − ERA5")))
        ax = Axis(fig[row, col], title = "$(label): $(t)", aspect = DataAspect())
        lim = col == 3 ? maximum(abs, filter(isfinite, mask(A))) : 0.0
        hm = heatmap!(ax, λ, φ, mask(A); colormap = col == 3 ? :balance : :viridis, colorrange = col == 3 ? (-lim, lim) : cr,
                      nan_color = :gray85)
        Colorbar(fig[row, col][1, 2], hm)
    end
end
for (col, (A, t)) in enumerate(((mag_m, "model"), (mag_e, "ERA5"), (mag_m .- mag_e, "model − ERA5")))
    ax = Axis(fig[3, col], title = "|period-mean τ| (N/m²): $(t)", aspect = DataAspect())
    lim = col == 3 ? maximum(abs, filter(isfinite, mask(A))) : maximum(filter(isfinite, mask(mag_e)))
    hm = heatmap!(ax, λ, φ, mask(A); colormap = col == 3 ? :balance : :viridis, colorrange = col == 3 ? (-lim, lim) : (0, lim),
                  nan_color = :gray85)
    Colorbar(fig[3, col][1, 2], hm)
end
out = PREFIX * "_fluxes_vs_era5.png"
save(out, fig; px_per_unit = 1.2)
println("saved ", out)
