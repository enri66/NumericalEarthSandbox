# Why is the PISCES column of bgc/column_compare.jl nearly empty? One column, the same forcing as column_compare.jl, and every MAB_PD_EVERY days
# the terms that decide the phytoplankton budget in the top cell and the bottom cell:
#   - the concentrations (nanophytoplankton P and diatoms D in mmol C/m³), the chlorophyll:carbon ratio (g/g) and the iron quota (μmol Fe/mmol C),
#   - the nutrient limitations of the growth: LN (nitrogen), LPO₄, LFe (iron quota), LSi (diatoms) and their minimum L,
#   - the gross growth rate μ (1/day, light and nutrients), the total loss rate (growth minus the tendency, i.e. mortality + grazing + exudation),
#     and the net tendency, all per unit of carbon,
#   - the tendencies of the zooplankton Z and M, and the iron concentration.
#   julia --project=<env with OceanBioME> bgc/pisces_diagnose.jl
# Environment: MAB_PD_DEPTH (20), MAB_PD_DAYS (45), MAB_PD_EVERY (3), MAB_PD_DT (60), MAB_PD_SUPPLY_DAYS (20), MAB_PD_FE (iron supply target, μmol/m³, 0.3),
# MAB_PD_PFE0 (initial iron quota of the phytoplankton, μmol Fe/mmol C, 0.01), MAB_PD_PO4 (0.5), MAB_COLUMN_LATITUDE (39.5), MAB_COLUMN_TRANSMISSIVITY (0.55), BGC_OUT.
using OceanBioME, Oceananigans, Printf, Statistics
using Oceananigans.Units
using Oceananigans.Fields: FunctionField, Field
using Oceananigans.Forcings: Relaxation
using Oceananigans.Biogeochemistry: biogeochemical_auxiliary_fields

const H       = parse(Float64, get(ENV, "MAB_PD_DEPTH", "20"))
const DAYS    = parse(Float64, get(ENV, "MAB_PD_DAYS", "45"))
const EVERY   = parse(Float64, get(ENV, "MAB_PD_EVERY", "3"))
const Δt      = parse(Float64, get(ENV, "MAB_PD_DT", "60"))
const SUPPLY  = parse(Float64, get(ENV, "MAB_PD_SUPPLY_DAYS", "20")) * 86400
const FE      = parse(Float64, get(ENV, "MAB_PD_FE", "0.3"))
const PFE0    = parse(Float64, get(ENV, "MAB_PD_PFE0", "0.01"))
const PO4     = parse(Float64, get(ENV, "MAB_PD_PO4", "0.5"))
const LAT     = parse(Float64, get(ENV, "MAB_COLUMN_LATITUDE", "39.5"))
const TRANSM  = parse(Float64, get(ENV, "MAB_COLUMN_TRANSMISSIVITY", "0.55"))
const OUT     = get(ENV, "BGC_OUT", ".")
mkpath(OUT)
io = open(joinpath(OUT, "pisces_diagnose.txt"), "w")
function say(args...)
    println(args...); flush(stdout)
    println(io, args...); flush(io)
end

const DAY0 = 32.0
function daily_mean_insolation(day_of_year, φ)
    δ = deg2rad(23.44) * sin(2π * (day_of_year - 81) / 365)
    h₀ = acos(clamp(-tan(deg2rad(φ)) * tan(δ), -1, 1))
    d⁻² = 1 + 0.033 * cos(2π * day_of_year / 365)
    return 1361 / π * d⁻² * (h₀ * sin(deg2rad(φ)) * sin(δ) + cos(deg2rad(φ)) * cos(δ) * sin(h₀))
end
@inline PAR⁰(x, y, t) = 0.43 * TRANSM * daily_mean_insolation(DAY0 + t / 86400, LAT)
@inline MLD(t) = 10 + 30 * exp(-t / (45 * 86400))
@inline κ(x, y, z, t) = ifelse(-z < MLD(t), 1e-2, 1e-5)
@inline T_target(x, y, z, t) = 5 + 7 * (1 - exp(-t / (70 * 86400))) * exp(z / 25)

Nz = clamp(round(Int, H / 2), 10, 60)
grid = RectilinearGrid(size = (1, 1, Nz), extent = (20, 20, H))
clock = Clock(; time = 0.0)
κ_field = FunctionField{Center, Center, Center}(κ, grid; clock)
supply(c) = SUPPLY > 0 ? Relaxation(rate = 1 / SUPPLY, target = c) : (x, y, z, t) -> 0.0
mld = Field{Center, Center, Nothing}(grid)
biogeochemistry = PISCES(; grid, latitude = PrescribedLatitude(LAT), surface_photosynthetically_active_radiation = PAR⁰,
                           mixed_layer_depth = mld, scale_negatives = true, invalid_fill_value = 0.0)
forcing = (; NO₃ = supply(4.0), PO₄ = supply(PO4), Si = supply(6.0), Fe = supply(FE), DIC = supply(2100.0), Alk = supply(2300.0),
             T = Relaxation(rate = 1 / 3600, target = T_target), S = Relaxation(rate = 1 / 3600, target = 33.0))
model = NonhydrostaticModel(grid; clock, closure = ScalarDiffusivity(ν = 1e-5, κ = κ_field), biogeochemistry, forcing)
set!(model, P = 0.3, PChl = 0.072, PFe = 0.3 * PFE0, D = 0.05, DChl = 0.012, DFe = 0.05 * PFE0, DSi = 0.005, Z = 0.1, M = 0.3, DOC = 2.0,
            POC = 2.0, SFe = 0.05, GOC = 1.0, BFe = 0.03, PSi = 0.05, CaCO₃ = 1e-10, NO₃ = 4.0, NH₄ = 0.2, PO₄ = PO4, Fe = FE,
            Si = 6.0, DIC = 2100.0, Alk = 2300.0, O₂ = 330.0, T = 5.0, S = 33.0)

simulation = Simulation(model; Δt, stop_time = DAYS * 86400)
mld_field = biogeochemistry.underlying_biogeochemistry.mixed_layer_depth
simulation.callbacks[:mld] = Callback(sim -> set!(mld_field, -min(MLD(sim.model.clock.time), H)), IterationInterval(1))

const U = biogeochemistry.underlying_biogeochemistry
const NANO, DIAT = U.phytoplankton.nano, U.phytoplankton.diatoms

function diagnose(sim)
    m = sim.model
    tr = m.tracers
    aux = biogeochemical_auxiliary_fields(m.biogeochemistry)
    t = m.clock.time
    day = 86400
    say(@sprintf("\nday %5.1f  PAR %.1f W/m²  mixed layer %.1f m  euphotic depth %.1f m", t / day, PAR⁰(0, 0, t), min(MLD(t), H), aux.zₑᵤ[1, 1, Nz]))
    for (label, k) in (("top cell", Nz), ("bottom cell", 1))
        say("  $label: NO₃ $(@sprintf("%.2f", tr.NO₃[1, 1, k])), NH₄ $(@sprintf("%.2f", tr.NH₄[1, 1, k])), PO₄ $(@sprintf("%.2f", tr.PO₄[1, 1, k])), Si $(@sprintf("%.2f", tr.Si[1, 1, k])), Fe $(@sprintf("%.3f", tr.Fe[1, 1, k])), T $(@sprintf("%.1f", tr.T[1, 1, k])), Z $(@sprintf("%.3f", tr.Z[1, 1, k])), M $(@sprintf("%.3f", tr.M[1, 1, k]))")
        for (name, ph, val, I, IChl, IFe) in (("nano", NANO, Val(:P), tr.P, tr.PChl, tr.PFe), ("diatoms", DIAT, Val(:D), tr.D, tr.DChl, tr.DFe))
            C, Chl, Fe = I[1, 1, k], IChl[1, 1, k], IFe[1, 1, k]
            L, LFe, LPO₄, LN = ph.nutrient_limitation(val, 1, 1, k, grid, U, ph, m.clock, tr, aux)
            μ = ph.growth_rate(val, 1, 1, k, grid, U, ph, m.clock, tr, aux, L) * day
            tend = m.biogeochemistry(1, 1, k, grid, val, m.clock, merge(tr, aux)) * day
            μC = (1 - ph.exudated_fraction) * μ                        # carbon growth rate, the first term of the tendency
            say(@sprintf("    %-8s C %6.3f  Chl:C %.4f  Fe quota %.4f  |  LN %.3f  LPO₄ %.3f  LFe %.3f  L %.3f  |  μ %.3f/d  loss %.3f/d  net %+.3f/d",
                         name, C, Chl / (12 * C + eps()), Fe / (C + eps()), LN, LPO₄, LFe, L, μ, μC - tend / (C + eps()), tend / (C + eps())))
        end
    end
    return nothing
end
simulation.callbacks[:diagnose] = Callback(diagnose, TimeInterval(EVERY * day))

say("PISCES column of $H m, $DAYS days, Δt $Δt s, supply time $(SUPPLY / 86400) days, Fe target $FE, PO₄ target $PO4, initial phytoplankton iron quota $PFE0")
diagnose(simulation)
run!(simulation)
