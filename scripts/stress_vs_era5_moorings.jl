# The surface stress the coupled model applies at the OOI Pioneer mooring cells, against simple ERA5-based estimates,
# from a script-05 run whose mooring output carries the surface wind and fluxes (MAB_MOORINGS=pioneer, d2398c8 or
# later). At each mooring cell, hourly:
#   model:     the stress NumericalEarth applies to the ocean, -(x_momentum, y_momentum) (N/m², into the ocean);
#   absolute:  ρₐ C_D |U| U from the model's own atmospheric wind U (Large & Pond 1981), as in column_cnsm.jl;
#   relative:  the same with U minus the model's surface current.
# Compared: the wind the model uses against ERA5 interpolated to the cell; mean |τ|; stress variance in the
# near-inertial band (16-22 h), clockwise and counterclockwise; and the near-inertial wind work on the surface current,
# W = ⟨τ_NI · u_NI⟩ / ρ₀ (W/kg·m), with u the model's top-cell current, for each stress.
# Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/res_test/r12s MAB_START_DATE=2019-08-29 julia --project=. scripts/stress_vs_era5_moorings.jl
include(joinpath(@__DIR__, "mab_analysis_common.jl"))
const FFTW = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "FFTW")]

const TAG      = get(ENV, "MAB_TAG", "r12s")
const PREFIX   = run_prefix(TAG)
const ERA5_DIR = joinpath(GDIR, "era5")
const BOX      = "_-76.0_-64.0_34.0_42.0.nc"
const BAND     = (16.0, 22.0)
const ρₐ, ρ₀   = 1.22, 1026.0

CD(U) = U < 11 ? 1.2e-3 : (0.49 + 0.065U) * 1e-3
lp(u, v) = (U = hypot(u, v); (ρₐ * CD(U) * U * u, ρₐ * CD(U) * U * v))

ds = NCD.Dataset(joinpath(ERA5_DIR, "10m_u_component_of_wind_ERA5HourlySingleLevel_2019-09-01T00$BOX"))
elon = Float64.(ds["longitude"][:]); elat = Float64.(ds["latitude"][:]); close(ds)
order = sortperm(elat); elat_sorted = elat[order]
function era5(name, var, t, λ, φ)
    f = joinpath(ERA5_DIR, "$(name)_ERA5HourlySingleLevel_$(Dates.format(t, "yyyy-mm-ddTHH"))$BOX")
    isfile(f) || return NaN
    ds = NCD.Dataset(f); a = Float64.(ds[var][:, :, 1])[:, order]; close(ds)
    return bilinear(a, elon, elat_sorted, λ, φ)
end

freq(n) = FFTW.fftfreq(n, 1.0)
function band(Z)                                              # complex series, near-inertial band (both rotations)
    Y = FFTW.fft(Z .- sum(Z) / length(Z)); fr = freq(length(Z))
    Y[.!(1 / BAND[2] .<= abs.(fr) .<= 1 / BAND[1])] .= 0
    return FFTW.ifft(Y)
end
function rotary(Z; clockwise)                                 # band variance of one rotation (clockwise: negative frequencies)
    Y = FFTW.fft(Z .- sum(Z) / length(Z)); P = abs2.(Y) ./ length(Z)^2; fr = freq(length(Z))
    return sum(P[[(clockwise ? -f : f) > 0 && 1 / BAND[2] <= abs(f) <= 1 / BAND[1] for f in fr]])
end

files = filter(f -> occursin(r"_moorings_rank\d+\.jld2$", f) && startswith(basename(f), basename(PREFIX) * "_moorings"),
               readdir(dirname(PREFIX); join = true))
for f in files, name in JLD2.jldopen(file -> filter(k -> file[k] isa JLD2.Group && haskey(file[k], "τx"), collect(keys(file))), f)
    m = JLD2.jldopen(f) do file
        t0 = DateTime(file["start_date"])
        (; t = [round(t0 + Second(round(Int, s)), Hour) for s in file["$name/time"]], cell = file["$name/cell"],
           ua = file["$name/ua"], va = file["$name/va"], τx = file["$name/τx"], τy = file["$name/τy"],
           us = file["$name/u"][end, :], vs = file["$name/v"][end, :])
    end
    keep = [true; diff(m.t) .> Hour(0)]
    t = m.t[keep]; ua = m.ua[keep]; va = m.va[keep]; us = m.us[keep]; vs = m.vs[keep]
    τm = -(m.τx[keep] .+ im .* m.τy[keep])                    # into the ocean
    τa = [complex(lp(ua[n], va[n])...) for n in eachindex(t)]
    τr = [complex(lp(ua[n] - us[n], va[n] - vs[n])...) for n in eachindex(t)]
    λ, φ = m.cell[3], m.cell[4]
    eu = [era5("10m_u_component_of_wind", "u10", s, λ, φ) for s in t]; ev = [era5("10m_v_component_of_wind", "v10", s, λ, φ) for s in t]
    ok = isfinite.(eu) .& isfinite.(ev)
    @printf("\n== %s: cell (%.3f, %.3f), %d hours %s - %s\n", name, λ, φ, length(t), t[1], t[end])
    @printf("   model wind vs ERA5 at the cell: mean speed %.2f vs %.2f m/s, rms vector difference %.2f m/s\n",
            mean(hypot.(ua[ok], va[ok])), mean(hypot.(eu[ok], ev[ok])), sqrt(mean((ua[ok] .- eu[ok]) .^ 2 .+ (va[ok] .- ev[ok]) .^ 2)))
    @printf("   mean surface current %.2f m/s\n", mean(hypot.(us, vs)))
    U = us .+ im .* vs; Ub = band(U)
    @printf("   %-10s %10s %22s %22s %18s\n", "stress", "mean |τ|", "NI clockwise var", "NI counterclockwise", "NI wind work")
    for (label, τ) in (("model", τm), ("absolute", τa), ("relative", τr))
        W = mean(real.(band(τ) .* conj.(Ub))) / ρ₀
        @printf("   %-10s %10.4f %22.3e %22.3e %18.3e\n", label, mean(abs.(τ)), rotary(τ; clockwise = true), rotary(τ; clockwise = false), W)
    end
    @printf("   model / absolute: mean |τ| %.2f, NI clockwise variance %.2f;  model / relative: %.2f, %.2f\n",
            mean(abs.(τm)) / mean(abs.(τa)), rotary(τm; clockwise = true) / rotary(τa; clockwise = true),
            mean(abs.(τm)) / mean(abs.(τr)), rotary(τm; clockwise = true) / rotary(τr; clockwise = true))
end
