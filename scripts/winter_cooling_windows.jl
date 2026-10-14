# Which stretch of a winter has the strongest surface cooling over the MAB box? A bulk estimate of the daily mean net
# surface heat loss over ocean points, from hourly ERA5 (10 m wind, 2 m temperature and dewpoint, surface pressure,
# downward shortwave and longwave) and GLORYS daily surface temperature on the ERA5 grid (nearest point):
#   sensible ρₐ cₚ C U (SST − T₂), latent ρₐ L C U (0.98 qₛₐₜ(SST) − q(T_d)), C = 1.2×10⁻³,
#   net longwave ε σ SST⁴ − ε LW↓ (ε = 0.97), shortwave −(1 − 0.06) SW↓.
# Prints the daily values and the windows of MAB_WINDOW days (default 14) with the largest mean loss, for the whole box
# and for deep water (GLORYS depth > 1000 m). Usage:
#   MAB_FROM=2019-12-01 MAB_TO=2020-03-31 julia --project=/t0/workdir/enrique/mpi05_ib scripts/winter_cooling_windows.jl
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const FROM = Date(get(ENV, "MAB_FROM", "2019-12-01")); const TO = Date(get(ENV, "MAB_TO", "2020-03-31"))
const WINDOW = parse(Int, get(ENV, "MAB_WINDOW", "14"))
const EDIR = joinpath(GDIR, "era5")
const ρₐ, cₚ, L, C, ε, σ = 1.25, 1004.0, 2.5e6, 1.2e-3, 0.97, 5.67e-8

efile(var, t) = joinpath(EDIR, "$(var)_ERA5HourlySingleLevel_$(Dates.format(t, "yyyy-mm-ddTHH"))_-76.0_-64.0_34.0_42.0.nc")
function era5(var, t)
    NCD.NCDataset(efile(var, t)) do ds
        name = only(filter(k -> !(k in ("longitude", "latitude", "valid_time", "time", "number", "expver")), keys(ds)))
        A = Array(ds[name])
        return Float64.(coalesce.(ndims(A) == 3 ? A[:, :, 1] : A, NaN))
    end
end
elon, elat = NCD.NCDataset(efile("2m_temperature", DateTime(FROM)), "r") do ds
    Array(ds["longitude"]), Array(ds["latitude"])
end
qsat(T, p) = (e = 611.2 * exp(17.67 * (T - 273.15) / (T - 29.65)); 0.622e / (p - 0.378e))

# GLORYS surface temperature and depth on the ERA5 grid
glon, glat, _, _ = read_glorys("thetao", FROM)
H = NCD.NCDataset(only(filter(f -> startswith(f, "deptho_GLORYSStatic"), readdir(GDIR))) |> f -> joinpath(GDIR, f)) do ds
    Float64.(coalesce.(Array(ds["deptho"]), NaN))
end
nearest(x, xs) = argmin(abs.(xs .- x))
gi = [nearest(x, glon) for x in elon]; gj = [nearest(y, glat) for y in elat]
deep = [H[gi[i], gj[j]] > 1000 for i in eachindex(elon), j in eachindex(elat)]

days = FROM:Day(1):TO
loss = Float64[]; loss_deep = Float64[]
for d in days
    _, _, _, GT = read_glorys("thetao", d)
    sst = [GT[gi[i], gj[j], 1] + 273.15 for i in eachindex(elon), j in eachindex(elat)]
    Q = zeros(size(sst))
    for h in 0:23
        t = DateTime(d) + Hour(h)
        U = hypot.(era5("10m_u_component_of_wind", t), era5("10m_v_component_of_wind", t))
        T2 = era5("2m_temperature", t); Td = era5("2m_dewpoint_temperature", t); p = era5("surface_pressure", t)
        SW = era5("surface_solar_radiation_downwards", t) ./ 3600; LW = era5("surface_thermal_radiation_downwards", t) ./ 3600
        Q .+= (ρₐ * cₚ * C .* U .* (sst .- T2) .+ ρₐ * L * C .* U .* (0.98 .* qsat.(sst, p) .- qsat.(Td, p)) .+
               ε * σ .* sst .^ 4 .- ε .* LW .- 0.94 .* SW) ./ 24
    end
    ok = isfinite.(Q)
    push!(loss, mean(Q[ok])); push!(loss_deep, mean(Q[ok .& deep]))
    @printf("%s  net loss %6.0f W/m² (deep water %6.0f)\n", d, loss[end], loss_deep[end]); flush(stdout)
end

for (label, x) in (("whole box", loss), ("deep water", loss_deep))
    m = [mean(x[k:k+WINDOW-1]) for k in 1:length(x)-WINDOW+1]
    order = sortperm(m; rev = true)
    println("\n$label: $(WINDOW)-day windows with the largest mean loss (non-overlapping)")
    taken = Int[]
    for k in order
        any(abs(k - t) < WINDOW for t in taken) && continue
        push!(taken, k)
        @printf("  %s to %s  %6.0f W/m²\n", days[k], days[k+WINDOW-1], m[k])
        length(taken) == 5 && break
    end
end
