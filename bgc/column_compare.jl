# Single-column comparison of LOBSTER and PISCES (OceanBioME): does either give more surface chlorophyll in shallow water than in deep water
# under the same forcing? A column of each water depth (MAB_COLUMN_DEPTHS), forced from 1 February by
#   - the daily-mean PAR at MAB_COLUMN_LATITUDE (top-of-atmosphere insolation × transmissivity MAB_COLUMN_TRANSMISSIVITY × 0.43 PAR fraction),
#   - an idealised mixed layer that shoals from 40 m (1 February) to 10 m (e^-1 time 45 days), the turbulent diffusivity 1e-2 m²/s in it
#     and 1e-5 below, so a column shallower than the mixed layer is mixed to the bottom (the shelf in winter),
#   - temperature warming from 5 to 12 °C at the surface (relaxed to the profile in 1 hour), salinity 33,
#   - a lateral supply of the same inorganic nutrients per unit volume at every depth: relaxation of NO₃ (and, in PISCES, PO₄, Si, Fe,
#     DIC and Alk) to the winter offshore values with the time scale MAB_COLUMN_SUPPLY_DAYS (0 = none).
# Both models have no sediment (sinking particles leave through the bottom, open_bottom), so the water column is the only difference between
# depths. The output (stdout and BGC_OUT/column_<model>.txt) is the surface chlorophyll (top cell), the water-column-mean chlorophyll,
# surface nitrate and, for PISCES, the diatom share and Chl:C, every MAB_COLUMN_EVERY days and a table of the means over days 30-90.
#   julia --project=<env with OceanBioME> bgc/column_compare.jl
# Environment: MAB_COLUMN_MODEL (lobster | pisces), MAB_COLUMN_DEPTHS ("10,20,40,60,100,200"), MAB_COLUMN_DAYS (120), MAB_COLUMN_DT (seconds, 60: the explicit vertical diffusion needs κ Δt / Δz² < 0.5, i.e. 0.15 for 2 m cells and κ = 1e-2),
# MAB_COLUMN_SUPPLY_DAYS (20), MAB_COLUMN_PO4 (0.5, PISCES: the half-saturation of the phosphate limitation is 0.8 so 0.5 limits strongly, see HANDOFF_BGC.md), MAB_COLUMN_BOTTOM_FLUX (0, mmol N/m²/day of NH₄ released into the bottom cell, in PISCES also PO₄/16; it is the only supply that depends on depth),
# MAB_COLUMN_LATITUDE (39.5), MAB_COLUMN_TRANSMISSIVITY (0.55), MAB_COLUMN_EVERY (10), BGC_OUT.
using OceanBioME, Oceananigans, Printf, Statistics
using Oceananigans.Units
using Oceananigans.Fields: FunctionField, Field
using Oceananigans.Forcings: Relaxation

const MODEL    = lowercase(get(ENV, "MAB_COLUMN_MODEL", "lobster"))
const DEPTHS   = parse.(Float64, split(get(ENV, "MAB_COLUMN_DEPTHS", "10,20,40,60,100,200"), ","))
const DAYS     = parse(Float64, get(ENV, "MAB_COLUMN_DAYS", "120"))
const Δt       = parse(Float64, get(ENV, "MAB_COLUMN_DT", "60"))
const SUPPLY   = parse(Float64, get(ENV, "MAB_COLUMN_SUPPLY_DAYS", "20")) * 86400
const LAT      = parse(Float64, get(ENV, "MAB_COLUMN_LATITUDE", "39.5"))
const TRANSM   = parse(Float64, get(ENV, "MAB_COLUMN_TRANSMISSIVITY", "0.55"))
const EVERY    = parse(Float64, get(ENV, "MAB_COLUMN_EVERY", "10"))
const PO4      = parse(Float64, get(ENV, "MAB_COLUMN_PO4", "0.5"))       # PISCES phosphate, initial and supply target (mmol P/m³)
const FLUX     = parse(Float64, get(ENV, "MAB_COLUMN_BOTTOM_FLUX", "0"))     # benthic NH₄ flux, mmol N/m²/day (0 = none)
const OUT      = get(ENV, "BGC_OUT", ".")
MODEL in ("lobster", "pisces") || error("MAB_COLUMN_MODEL must be lobster or pisces")
mkpath(OUT)
io = open(joinpath(OUT, "column_$(MODEL).txt"), "w")
function say(args...)
    println(args...); flush(stdout)
    println(io, args...); flush(io)
end

# ---------------- forcing ----------------
const DAY0 = 32.0                                         # day of the year of 1 February
function daily_mean_insolation(day_of_year, φ)             # W/m², top of the atmosphere
    δ = deg2rad(23.44) * sin(2π * (day_of_year - 81) / 365)
    h₀ = acos(clamp(-tan(deg2rad(φ)) * tan(δ), -1, 1))
    d⁻² = 1 + 0.033 * cos(2π * day_of_year / 365)
    return 1361 / π * d⁻² * (h₀ * sin(deg2rad(φ)) * sin(δ) + cos(deg2rad(φ)) * cos(δ) * sin(h₀))
end
@inline PAR⁰(x, y, t) = 0.43 * TRANSM * daily_mean_insolation(DAY0 + t / 86400, LAT)
@inline MLD(t) = 10 + 30 * exp(-t / (45 * 86400))
@inline κ(x, y, z, t) = ifelse(-z < MLD(t), 1e-2, 1e-5)
@inline T_target(x, y, z, t) = 5 + 7 * (1 - exp(-t / (70 * 86400))) * exp(z / 25)

# ---------------- one column ----------------
function run_column(H)
    Nz = clamp(round(Int, H / 2), 10, 60)
    grid = RectilinearGrid(size = (1, 1, Nz), extent = (20, 20, H))
    clock = Clock(; time = 0.0)
    κ_field = FunctionField{Center, Center, Center}(κ, grid; clock)
    supply(c) = SUPPLY > 0 ? Relaxation(rate = 1 / SUPPLY, target = c) : (x, y, z, t) -> 0.0
    # benthic regeneration: FLUX (mmol N/m²/day) enters the bottom cell as NH₄ (and PO₄/16 in PISCES), so the supply per unit volume goes as 1/H
    Δz = H / Nz
    bottom_source(x, y, z, t, p) = ifelse(z < p.zb, p.s, 0.0)
    bottom(factor) = Forcing(bottom_source; parameters = (zb = -H + Δz, s = factor * FLUX / Δz / 86400))
    if MODEL == "lobster"
        biogeochemistry = LOBSTER(grid; surface_PAR = PAR⁰, scale_negatives = true, open_bottom = true)
        forcing = (; NO₃ = supply(4.0), (FLUX > 0 ? (; NH₄ = bottom(1.0)) : NamedTuple())...)
        auxiliary_fields = NamedTuple()                  # LOBSTER without carbonate chemistry needs neither T nor S
    else
        mld = Field{Center, Center, Nothing}(grid)
        biogeochemistry = PISCES(; grid, latitude = PrescribedLatitude(LAT), surface_photosynthetically_active_radiation = PAR⁰,
                                   mixed_layer_depth = mld, scale_negatives = true, invalid_fill_value = 0.0)
        forcing = (; NO₃ = supply(4.0), PO₄ = FLUX > 0 ? (supply(PO4), bottom(1 / 16)) : supply(PO4), Si = supply(6.0), Fe = supply(0.3),
                     DIC = supply(2100.0), Alk = supply(2300.0),
                     T = Relaxation(rate = 1 / 3600, target = T_target), S = Relaxation(rate = 1 / 3600, target = 33.0),
                     (FLUX > 0 ? (; NH₄ = bottom(1.0)) : NamedTuple())...)
        auxiliary_fields = NamedTuple()
    end
    model = NonhydrostaticModel(grid; clock, closure = ScalarDiffusivity(ν = 1e-5, κ = κ_field), biogeochemistry, forcing, auxiliary_fields)
    if MODEL == "lobster"
        set!(model, P = 0.05, Z = 0.05, NO₃ = 4.0, NH₄ = 0.1)
    else
        set!(model, P = 0.3, PChl = 0.072, PFe = 0.003, D = 0.05, DChl = 0.012, DFe = 0.0005, DSi = 0.005, Z = 0.1, M = 0.3, DOC = 2.0,
                    POC = 2.0, SFe = 0.05, GOC = 1.0, BFe = 0.03, PSi = 0.05, CaCO₃ = 1e-10, NO₃ = 4.0, NH₄ = 0.2, PO₄ = PO4, Fe = 0.3,
                    Si = 6.0, DIC = 2100.0, Alk = 2300.0, O₂ = 330.0, T = 5.0, S = 33.0)
    end
    simulation = Simulation(model; Δt, stop_time = DAYS * 86400)
    if MODEL == "pisces"
        mld_field = biogeochemistry.underlying_biogeochemistry.mixed_layer_depth
        simulation.callbacks[:mld] = Callback(sim -> set!(mld_field, -min(MLD(sim.model.clock.time), H)), IterationInterval(1))
    end
    if haskey(ENV, "MAB_COLUMN_DEBUG")
        debug_names = MODEL == "pisces" ? (:T, :S, :DIC, :Alk, :Si, :CaCO₃, :NO₃, :P) : (:P, :NO₃)
        show_fields(sim) = say(@sprintf("  iteration %d (%.2f h): ", sim.model.clock.iteration, sim.model.clock.time / 3600),
                               join(("$n = $(round.(Array(interior(sim.model.tracers[n]))[1, 1, [1, end]], sigdigits = 5))" for n in debug_names), "; "))
        simulation.callbacks[:debug] = Callback(show_fields, IterationInterval(6))
    end
    days = Float64[]; rows = NamedTuple[]
    function record(sim)
        tr = sim.model.tracers
        chl = MODEL == "lobster" ? 1.31 .* Array(interior(tr.P))[1, 1, :] : Array(interior(tr.PChl))[1, 1, :] .+ Array(interior(tr.DChl))[1, 1, :]
        no3 = Array(interior(tr.NO₃))[1, 1, :]
        extra = MODEL == "pisces" ? (diatoms = Array(interior(tr.D))[1, 1, end] / (Array(interior(tr.P))[1, 1, end] + Array(interior(tr.D))[1, 1, end]),
                                     chl_C = chl[end] / (12 * (Array(interior(tr.P))[1, 1, end] + Array(interior(tr.D))[1, 1, end]))) :
                                    (diatoms = NaN, chl_C = 1.31 / (12 * 6.56))
        push!(days, sim.model.clock.time / 86400)
        push!(rows, (; surface = chl[end], mean = mean(chl), no3 = no3[end], extra...))
    end
    simulation.callbacks[:record] = Callback(record, TimeInterval(1day))
    run!(simulation)
    return days, rows
end

say("model $MODEL, column depths $DEPTHS m, $DAYS days, Δt $(Δt) s, supply time scale $(SUPPLY / 86400) days, latitude $(LAT)°, transmissivity $TRANSM")
say("daily-mean PAR (W/m²) on days 0, 30, 60, 90: ", join((@sprintf("%.1f", PAR⁰(0, 0, d * 86400.0)) for d in (0, 30, 60, 90)), ", "))
results = Dict{Float64, Any}()
for H in DEPTHS
    say("\ncolumn of $H m")
    days, rows = run_column(H)
    results[H] = (days, rows)
    for n in 1:round(Int, EVERY):length(days)
        r = rows[n]
        say(@sprintf("  day %5.1f  surface chlorophyll %7.3f  column mean %7.3f  surface nitrate %6.2f  diatom share %.2f  Chl:C %.4f", days[n], r.surface, r.mean, r.no3, r.diatoms, r.chl_C))
    end
end
say("\nmeans over days 30-90 (mg Chl/m³)")
say("  depth (m)  surface chlorophyll  column-mean chlorophyll  surface nitrate  diatom share  Chl:C (g/g)")
for H in DEPTHS
    days, rows = results[H]
    sel = findall(d -> 30 ≤ d ≤ 90, days)
    m(f) = isempty(sel) ? NaN : mean(getfield.(rows[sel], f))
    say(@sprintf("  %8.0f  %19.3f  %22.3f  %15.2f  %12.2f  %.4f", H, m(:surface), m(:mean), m(:no3), m(:diatoms), m(:chl_C)))
end
