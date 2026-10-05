# A single CATKE column at the OOI Pioneer central surface mooring (CP01CNSM, 40.13°N 70.78°W), 29 Aug - 31 Oct 2019:
# how much near-inertial energy does CATKE's vertical mixing alone leave in the upper ocean, against the 3D model, a
# slab mixed layer and the ADCP? The column has the 3D run's vertical grid down to its depth at the mooring, starts
# from that run's T and S there (at rest), and is forced by
#   the ERA5 10 m wind at the buoy, through τ = ρₐ C_D |ΔU| ΔU (Large & Pond 1981), with ΔU the wind itself
#       ("absolute") or the wind minus the surface current ("relative");
#   the buoy's own net heat flux (METBK A and B, positive upward), applied at the surface (no salinity flux).
# Variants (COLUMN_VARIANTS, comma separated): catke, catke_relative, combo (the combo run's coefficients),
# combo_relative; each column's hourly profiles are saved to <tag>_column_cnsm_<variant>.jld2. The near-inertial kinetic energy (16-22 h band) is averaged over 12-72 m, the window of the
# mooring comparison, and printed beside the 3D run's own at the mooring cell (MAB_TAG, MAB_MOORINGS=pioneer).
# Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/res_test/r12 MAB_START_DATE=2019-08-29 OOI_DIR=/t0/workdir/enrique/Data/OOI/pioneer \
#       julia --project=. scripts/column_cnsm.jl
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEMixingLength, CATKEEquation
const FFTW = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "FFTW")]

const TAG      = get(ENV, "MAB_TAG", "r12")
const PREFIX   = run_prefix(TAG)
const OOI_DIR  = get(ENV, "OOI_DIR", joinpath(homedir(), "Data", "OOI", "pioneer"))
const ERA5_DIR = joinpath(GDIR, "era5")
const BOX      = "_-76.0_-64.0_34.0_42.0.nc"
const VARIANTS = split(get(ENV, "COLUMN_VARIANTS", "catke,catke_relative,combo,combo_relative"), ",")
const LON, LAT = -70.7783, 40.1333
const ρₐ, ρ₀, cₚ = 1.22, 1026.0, 3991.9
const DAYS = 63
const T0 = start_date
const BAND = (16.0, 22.0)
const WINDOW = (12.0, 72.0)                                   # m, as in moorings_vs_ooi.jl for the shelf moorings

# ---------------- the 3D run's column at the mooring ----------------
moor = let f = first(filter(f -> occursin(r"_moorings_rank\d+\.jld2$", f) && startswith(basename(f), basename(PREFIX) * "_moorings"),
                             readdir(dirname(PREFIX); join = true)))
    JLD2.jldopen(f) do file
        (; zc = file["z_center"], zf = file["z_face"], cell = file["CNSM/cell"],
           t = [round(DateTime(file["start_date"]) + Second(round(Int, s)), Hour) for s in file["CNSM/time"]],
           u = file["CNSM/u"], v = file["CNSM/v"], T = file["CNSM/T"], S = file["CNSM/S"])
    end
end
depth = moor.cell[5]
kbot = findfirst(z -> z > -depth, moor.zf)                    # deepest face above the model bottom
zfaces = moor.zf[kbot:end]; zfaces[1] = -depth
kc = kbot:length(moor.zc)
@printf("3D column at (%.3f, %.3f), %.0f m deep: %d levels, top cell %.1f m\n", moor.cell[3], moor.cell[4], depth, length(kc),
        zfaces[end] - zfaces[end-1])

# ---------------- forcing ----------------
NCD5 = NCD
ds = NCD5.Dataset(joinpath(ERA5_DIR, "10m_u_component_of_wind_ERA5HourlySingleLevel_2019-09-01T00$BOX"))
elon = Float64.(ds["longitude"][:]); elat = Float64.(ds["latitude"][:]); close(ds)
order = sortperm(elat); elat_sorted = elat[order]
function era5(name, var, t)
    ds = NCD5.Dataset(joinpath(ERA5_DIR, "$(name)_ERA5HourlySingleLevel_$(Dates.format(t, "yyyy-mm-ddTHH"))$BOX"))
    a = Float64.(ds[var][:, :, 1])[:, order]; close(ds)
    return bilinear(a, elon, elat_sorted, LON, LAT)
end
hours = collect(T0:Hour(1):T0 + Day(DAYS) + Hour(1))
Ua = [era5("10m_u_component_of_wind", "u10", t) for t in hours]
Va = [era5("10m_v_component_of_wind", "v10", t) for t in hours]

function buoy_heat_flux()
    sums = Dict{DateTime, NTuple{2, Float64}}()
    for d in ("ooi-cp01cnsm-sbd11-06-metbka000", "ooi-cp01cnsm-sbd12-06-metbka000")
        lines = readlines(joinpath(OOI_DIR, d * "_flux.csv"))
        for line in lines[3:end]
            v = split(line, ",")
            q = tryparse(Float64, v[2]); (isnothing(q) || !isfinite(q) || v[3] == "4") && continue
            t = round(DateTime(v[1][1:19]), Hour); s = get(sums, t, (0.0, 0.0)); sums[t] = (s[1] + q, s[2] + 1)
        end
    end
    Q = [haskey(sums, t) ? sums[t][1] / sums[t][2] : NaN for t in hours]
    good = findall(isfinite, Q)
    for n in eachindex(Q)                                     # short gaps: linear in time
        isfinite(Q[n]) && continue
        a = findlast(<(n), good); b = findfirst(>(n), good)
        Q[n] = isnothing(a) ? Q[good[b]] : isnothing(b) ? Q[good[a]] : Q[good[a]] + (Q[good[b]] - Q[good[a]]) * (n - good[a]) / (good[b] - good[a])
    end
    return Q
end
Qup = buoy_heat_flux()
@printf("forcing: ERA5 mean 10 m speed %.2f m/s; buoy net heat flux mean %.1f W/m² (positive upward)\n",
        mean(hypot.(Ua, Va)), mean(Qup))

# Hourly series → value at time t (s since T0), linear
@inline function at(series, t)
    x = t / 3600 + 1; n = clamp(floor(Int, x), 1, length(series) - 1); w = clamp(x - n, 0, 1)
    return @inbounds (1 - w) * series[n] + w * series[n+1]
end
CD(U) = U < 11 ? 1.2e-3 : (0.49 + 0.065U) * 1e-3

# ---------------- columns ----------------
function catke(variant)
    kw = startswith(variant, "combo") ? (; Cᵉc = 0.0, Cˡᵒc = 0.1845, Cʰⁱc = 0.049) : (;)
    tke = startswith(variant, "combo") ? (; Cᵂu★ = 1.59) : (;)
    return CATKEVerticalDiffusivity(VerticallyImplicitTimeDiscretization(); mixing_length = CATKEMixingLength(; kw...),
                                    turbulent_kinetic_energy_equation = CATKEEquation(; Cᵂϵ = 1.0, tke...))
end

function run_column(variant)
    grid = RectilinearGrid(size = length(kc), z = zfaces, topology = (Flat, Flat, Bounded))
    relative = endswith(variant, "relative")
    p = (; Ua, Va, relative)
    @inline Δ(t, u, v, p) = (at(p.Ua, t) - p.relative * u, at(p.Va, t) - p.relative * v)
    @inline function τx(t, u, v, p); (a, b) = Δ(t, u, v, p); U = hypot(a, b); -ρₐ * CD(U) * U * a / ρ₀; end
    @inline function τy(t, u, v, p); (a, b) = Δ(t, u, v, p); U = hypot(a, b); -ρₐ * CD(U) * U * b / ρ₀; end
    u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(τx; field_dependencies = (:u, :v), parameters = p))
    v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(τy; field_dependencies = (:u, :v), parameters = p))
    T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition((t, p) -> at(p, t) / (ρ₀ * cₚ); parameters = Qup))
    model = HydrostaticFreeSurfaceModel(grid; coriolis = FPlane(latitude = LAT), closure = catke(variant),
                                        buoyancy = SeawaterBuoyancy(equation_of_state = SWP.TEOS10.TEOS10EquationOfState()),
                                        tracers = (:T, :S), boundary_conditions = (u = u_bcs, v = v_bcs, T = T_bcs))
    set!(model, T = reshape(moor.T[kc, 1], 1, 1, :), S = reshape(moor.S[kc, 1], 1, 1, :))
    sim = Simulation(model; Δt = 300.0, stop_time = DAYS * 86400.0)
    out = (t = Float64[], u = Vector{Float64}[], v = Vector{Float64}[], T = Vector{Float64}[], S = Vector{Float64}[],
           κu = Vector{Float64}[], e = Vector{Float64}[])
    function record(sim)
        push!(out.t, sim.model.clock.time)
        push!(out.u, Array(interior(sim.model.velocities.u))[1, 1, :]); push!(out.v, Array(interior(sim.model.velocities.v))[1, 1, :])
        push!(out.T, Array(interior(sim.model.tracers.T))[1, 1, :]); push!(out.S, Array(interior(sim.model.tracers.S))[1, 1, :])
        push!(out.κu, Array(interior(sim.model.closure_fields.κu))[1, 1, :]); push!(out.e, Array(interior(sim.model.tracers.e))[1, 1, :])
    end
    add_callback!(sim, record, TimeInterval(3600.0))
    run!(sim)
    return out
end

# ---------------- near-inertial energy over the window ----------------
grid4 = collect(WINDOW[1]:4.0:WINDOW[2])
function on_grid(depths, U)                                    # depths increasing; U (levels × times)
    [begin k = clamp(searchsortedlast(depths, z), 1, length(depths) - 1); w = (z - depths[k]) / (depths[k+1] - depths[k])
         (1 - w) * U[k, n] + w * U[k+1, n] end for z in grid4, n in axes(U, 2)]
end
function bandpass(x)
    Y = FFTW.fft(x .- mean(x)); fr = FFTW.fftfreq(length(x), 1.0)
    Y[.!(1 / BAND[2] .<= abs.(fr) .<= 1 / BAND[1])] .= 0
    return real(FFTW.ifft(Y))
end
function near_inertial_KE(depths, U, V)
    u = on_grid(depths, U); v = on_grid(depths, V)
    ub = reduce(hcat, [bandpass(u[m, :]) for m in axes(u, 1)])'; vb = reduce(hcat, [bandpass(v[m, :]) for m in axes(v, 1)])'
    return vec(mean(0.5 .* (ub .^ 2 .+ vb .^ 2); dims = 1))
end
storm(times, a, b) = [a <= t < b for t in times]
STORMS = (("Dorian, 6-12 Sep", DateTime(2019, 9, 6), DateTime(2019, 9, 13)), ("17-23 Oct", DateTime(2019, 10, 17), DateTime(2019, 10, 24)))

depth_c = reverse(-moor.zc[kc])
keep = [true; diff(moor.t) .> Hour(0)]
t3 = moor.t[keep]
ke3 = near_inertial_KE(depth_c, reverse(moor.u[kc, keep]; dims = 1), reverse(moor.v[kc, keep]; dims = 1))
mld(T, S) = first(column_depths(depth_c, reverse(T), reverse(S)))
mld3 = [mld(moor.T[kc, n], moor.S[kc, n]) for n in findall(keep)]

results = Dict{String, Any}()
for variant in VARIANTS
    @info "column: $variant"
    out = run_column(String(variant))
    times = [T0 + Second(round(Int, s)) for s in out.t]
    # hourly profiles, bottom to top, for comparison with the 3D run's mooring column (z_face: the column's faces)
    JLD2.jldsave(PREFIX * "_column_cnsm_$(variant).jld2"; start_date = string(T0), time = out.t, z_face = zfaces,
                 u = reduce(hcat, out.u), v = reduce(hcat, out.v), T = reduce(hcat, out.T), S = reduce(hcat, out.S),
                 κu = reduce(hcat, out.κu), e = reduce(hcat, out.e))
    ke = near_inertial_KE(depth_c, reverse(reduce(hcat, out.u); dims = 1), reverse(reduce(hcat, out.v); dims = 1))
    results[variant] = (; times, ke, mld = [mld(out.T[n], out.S[n]) for n in eachindex(out.t)])
end

println("\nNear-inertial KE (16-22 h), mean over $(Int(WINDOW[1]))-$(Int(WINDOW[2])) m (J/kg): record / Dorian / 17-23 Oct;  median MLD Sep / Oct (m)")
function report(label, times, ke, h)
    sep = [Dates.month(t) == 9 for t in times]; oct = [Dates.month(t) == 10 for t in times]
    @printf("  %-22s %.2e / %.2e / %.2e    %5.1f / %5.1f\n", label, mean(ke), [mean(ke[storm(times, a, b)]) for (_, a, b) in STORMS]...,
            median(filter(isfinite, h[sep])), median(filter(isfinite, h[oct])))
end
report("3D model ($(basename(PREFIX)))", t3, ke3, mld3)
for variant in VARIANTS
    r = results[variant]; report("column $variant", r.times, r.ke, r.mld)
end
println("  (ADCP over the same window and record: 4.04e-03 J/kg; slab with ERA5 stress, H = 40 m: 3.7e-03)")

fig = Figure(size = (1400, 800), fontsize = 14)
td(ts) = [Dates.value(t - T0) / 86_400_000 for t in ts]
ax = Axis(fig[1, 1], ylabel = "near-inertial KE, $(Int(WINDOW[1]))-$(Int(WINDOW[2])) m (J/kg)", title = "CNSM: CATKE column vs 3D model")
lines!(ax, td(t3), ke3; color = :black, linewidth = 2, label = "3D model")
for (variant, color) in zip(VARIANTS, (:royalblue, :deepskyblue, :firebrick, :orange))
    lines!(ax, td(results[variant].times), results[variant].ke; color, label = "column $variant")
end
axislegend(ax; position = :lt)
ax = Axis(fig[2, 1], ylabel = "mixed-layer depth (m)", xlabel = "days since $(Dates.format(T0, "yyyy-mm-dd"))", yreversed = true)
lines!(ax, td(t3), mld3; color = :black, linewidth = 2)
for (variant, color) in zip(VARIANTS, (:royalblue, :deepskyblue, :firebrick, :orange))
    lines!(ax, td(results[variant].times), results[variant].mld; color)
end
out = PREFIX * "_column_cnsm.png"
save(out, fig; px_per_unit = 1.2)
println("saved ", out)
