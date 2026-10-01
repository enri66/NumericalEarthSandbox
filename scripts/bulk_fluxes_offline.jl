# NumericalEarth's bulk turbulent fluxes computed offline from ERA5's own surface state and compared with ERA5's own
# fluxes, which separates the bulk formulas from a run's SST and currents. This builds script 04's coupled model (the
# same flux code the runs use, with flooded ERA5) without time stepping it; every hour it sets the clocks, puts
# ERA5's SST (ERA5_SURFACE=sea_surface_temperature, the default) or skin temperature (ERA5_SURFACE=skin_temperature)
# into the ocean's top cells with zero surface currents, recomputes the fluxes with `update_state!`, and compares them
# with ERA5's mean fluxes over that hour, interpolated to the model's cells from ERA5 points with no land. Results are
# overall and binned by air-sea temperature difference and by wind speed. Needs, in the ERA5 cache, the forcing fields,
# ERA5's mean_surface_* fluxes, sea_surface_temperature or skin_temperature, and land_sea_mask. Usage:
#   ERA5_SURFACE=skin_temperature ERA5_DAYS=30 julia --project=. scripts/bulk_fluxes_offline.jl
const NDAYS   = parse(Int, get(ENV, "ERA5_DAYS", "30"))
const SURFACE = get(ENV, "ERA5_SURFACE", "sea_surface_temperature")

# Script 04's model, built but not run (its output writers are created but never written to)
ENV["MAB_DAYS"] = string(NDAYS)
ENV["MAB_TAG"]  = joinpath(tempdir(), "bulk_fluxes_offline")
script04 = read(joinpath(@__DIR__, "04_mab_glorys_tides_reservoirs.jl"), String)
include_string(Main, script04[1:findfirst("run!(simulation", script04).start - 1],
               joinpath(@__DIR__, "04_mab_glorys_tides_reservoirs.jl"))

using Oceananigans.TimeSteppers: update_state!

const NCD5 = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "NCDatasets")]
const ERA5_DIR = joinpath(DATA_DIR, "era5")
const BOX = "_-76.0_-64.0_34.0_42.0.nc"
era5_file(cds, date) = joinpath(ERA5_DIR, "$(cds)_ERA5HourlySingleLevel_$(Dates.format(date, "yyyy-mm-ddTHH"))$BOX")

# ERA5 stores latitude north to south; everything here runs south to north
function read2d(cds, var, date)
    ds = NCD5.Dataset(era5_file(cds, date))
    a = [ismissing(x) ? NaN : Float64(x) for x in ds[var][:, :, 1]]
    close(ds)
    return a[:, end:-1:1]
end
lsm_ds = NCD5.Dataset(joinpath(ERA5_DIR, "land_sea_mask_ERA5HourlySingleLevel_2019-01-01T00$BOX"))
elon = Float64.(lsm_ds["longitude"][:]); elat = reverse(Float64.(lsm_ds["latitude"][:]))
eland = Float64.(lsm_ds["lsm"][:, end:-1:1, 1]) .> 0
close(lsm_ds)

# Bilinear interpolation from ERA5 points with no land; NaN unless all four surrounding points are ocean
function to_model(A)
    B = copy(A); B[eland] .= NaN
    out = fill(NaN, Nλ, Nφ)
    for i in 1:Nλ, j in 1:Nφ
        bh[i, j] < 0 || continue
        x, y = λc[i], φc[j]
        a = clamp(searchsortedlast(elon, x), 1, length(elon) - 1); b = clamp(searchsortedlast(elat, y), 1, length(elat) - 1)
        wx = (x - elon[a]) / (elon[a+1] - elon[a]); wy = (y - elat[b]) / (elat[b+1] - elat[b])
        w = ((1 - wx) * (1 - wy), wx * (1 - wy), (1 - wx) * wy, wx * wy); v = (B[a, b], B[a+1, b], B[a, b+1], B[a+1, b+1])
        all(isfinite, v) && (out[i, j] = sum(w .* v))
    end
    return out
end

ocean_model = ocean.model
T, u, v = ocean_model.tracers.T, ocean_model.velocities.u, ocean_model.velocities.v
ao = model.interfaces.atmosphere_ocean_interface.fluxes
surface_var = SURFACE == "skin_temperature" ? "skt" : "sst"

names = (:sensible, :latent, :τx, :τy)
era5 = (sensible = ("mean_surface_sensible_heat_flux", "avg_ishf"), latent = ("mean_surface_latent_heat_flux", "avg_slhtf"),
        τx = ("mean_eastward_turbulent_surface_stress", "avg_iews"), τy = ("mean_northward_turbulent_surface_stress", "avg_inss"))
model_field = (sensible = ao.sensible_heat, latent = ao.latent_heat, τx = ao.x_momentum, τy = ao.y_momentum)
# One record per ocean cell and hour where ERA5 and the model are both defined for all four fluxes
rec = Dict(n => (Float64[], Float64[]) for n in names)
ΔT = Float64[]; U = Float64[]

for h in 1:24NDAYS
    t = 3600.0 * h
    date = start_date + Hour(h)
    # The prescribed atmosphere and radiation advance (and reload their data windows) through their own time steps;
    # the ocean is not stepped, only its clock and the coupled model's are moved along
    time_step!(atmosphere, 3600); time_step!(radiation, 3600)
    for c in (model.clock, ocean_model.clock)
        c.time = t; c.iteration = h
    end
    ts = to_model(read2d(SURFACE, surface_var, date) .- 273.15)
    ok = isfinite.(ts)
    Tdata = Array(interior(T)); top = view(Tdata, :, :, size(Tdata, 3)); top[ok] .= ts[ok]
    set!(T, Tdata)
    set!(u, 0); set!(v, 0)                        # no currents (only the top cell enters the fluxes)
    update_state!(model)

    m = Dict(n => -Array(interior(model_field[n]))[1:Nλ, 1:Nφ, 1] for n in names)   # NumericalEarth: positive upward
    e = Dict(n => to_model(read2d(era5[n]..., date)) for n in names)               # ERA5: mean over the hour ending at `date`
    k = ok .& reduce(.&, [isfinite.(m[n]) .& isfinite.(e[n]) for n in names])
    for n in names
        append!(rec[n][1], m[n][k]); append!(rec[n][2], e[n][k])
    end
    ta = to_model(read2d("2m_temperature", "t2m", date) .- 273.15)
    U10 = to_model(hypot.(read2d("10m_u_component_of_wind", "u10", date), read2d("10m_v_component_of_wind", "v10", date)))
    append!(ΔT, (ta .- ts)[k]); append!(U, U10[k])
    h % 24 == 0 && @printf("day %d done\n", h ÷ 24)
end

# ---------------- report ----------------
println("\nNumericalEarth bulk fluxes (script 04's coupled model) from ERA5's $(SURFACE), no currents, vs ERA5's own fluxes, ",
        "hourly at ocean cells ($(length(ΔT)) cell-hours; ERA5 sign convention, positive into the ocean)")
for n in names
    mm, ee = rec[n]
    @printf("  %-9s model %+9.4f  ERA5 %+9.4f  model/ERA5 %.3f  rms diff %.4f  corr %.3f\n", n, mean(mm), mean(ee), mean(mm) / mean(ee),
            sqrt(mean((mm .- ee) .^ 2)), cor(mm, ee))
end
mτ = hypot.(rec[:τx][1], rec[:τy][1]); eτ = hypot.(rec[:τx][2], rec[:τy][2])
@printf("  |τ|       model %.4f  ERA5 %.4f  model/ERA5 %.3f\n", mean(mτ), mean(eτ), mean(mτ) / mean(eτ))

S, L = rec[:sensible], rec[:latent]
println("\nBy air - sea temperature difference (stable when positive):")
@printf("  %-12s %8s %22s %22s\n", "ΔT (°C)", "count", "sensible model / ERA5", "latent model / ERA5")
for (lo, hi) in ((-Inf, -4), (-4, -2), (-2, -1), (-1, 0), (0, 1), (1, 2), (2, 4), (4, Inf))
    k = findall(x -> lo <= x < hi, ΔT)
    length(k) > 50 || continue
    @printf("  %5.0f..%-5.0f %8d %10.1f / %-10.1f %10.1f / %-10.1f\n", lo, hi, length(k),
            mean(S[1][k]), mean(S[2][k]), mean(L[1][k]), mean(L[2][k]))
end

println("\nBy 10 m wind speed:")
@printf("  %-12s %8s %24s %26s\n", "U10 (m/s)", "count", "|τ| model / ERA5 (N/m²)", "latent model / ERA5 (W/m²)")
for (lo, hi) in ((0, 3), (3, 6), (6, 9), (9, 12), (12, 15), (15, Inf))
    k = findall(x -> lo <= x < hi, U)
    length(k) > 50 || continue
    @printf("  %5.0f..%-5.0f %8d %11.4f / %-10.4f %12.1f / %-10.1f\n", lo, hi, length(k), mean(mτ[k]), mean(eτ[k]),
            mean(L[1][k]), mean(L[2][k]))
end
