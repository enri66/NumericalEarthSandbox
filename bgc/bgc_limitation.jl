# Growth-limitation diagnostics for LOBSTER's PhytoZoo (OceanBioME 0.19 and 0.20, same formulas), included by mab_bgc.jl with MAB_BGC_LIMITATION=true.
# The phytoplankton growth rate of the model is μ = μ₀ · L_light · L_nutrient · L_temperature · P, with
#   L_light       = PAR / (k_PAR + PAR)                                   (MondoLightLimitation, the PhytoZoo default; PAR in W/m²)
#   L_NO₃         = NO₃ exp(-ψ NH₄) / (NO₃ + k_NO₃)
#   L_NH₄         = max(0, NH₄ / (k_NH₄ + NH₄))
#   L_nutrient    = (L_NO₃ + L_NH₄) / 2                                   (so the maximum rate is μ₀ / 2)
#   L_temperature = Q₁₀^(T / 10), or 1 when the model has no Q₁₀ (the default of PhytoZoo, and of LOBSTER)
# `limitation_outputs` returns, as live operations on the model's own fields and parameters, PAR and each term, and the specific
# growth rate μ / P in 1/day. The formulas are copied from OceanBioME's phyto_zoo.jl (phytoplankton_growth, nitrogen_limitation); they
# must be re-checked against it when OceanBioME changes. Also saved: NO₃ and NH₄ (mmol/m³) and `daylight`, the fraction of the samples
# with PAR above 1 W/m², so that the light term of the daytime is light / daylight (the daily mean of the light term includes the night). Written by a daily average of hourly samples (see mab_bgc.jl).
using Oceananigans.AbstractOperations: KernelFunctionOperation
using Oceananigans.Biogeochemistry: biogeochemical_auxiliary_fields

@inline light_term(i, j, k, grid, PAR, kPAR) = @inbounds PAR[i, j, k] / (kPAR + PAR[i, j, k])
@inline nitrate_term(i, j, k, grid, NO₃, NH₄, kNO₃, ψ) = @inbounds NO₃[i, j, k] * exp(-ψ * NH₄[i, j, k]) / (NO₃[i, j, k] + kNO₃)
@inline ammonium_term(i, j, k, grid, NH₄, kNH₄) = @inbounds max(0, NH₄[i, j, k] / (kNH₄ + NH₄[i, j, k]))
@inline daylight_term(i, j, k, grid, PAR) = @inbounds ifelse(PAR[i, j, k] > 1, one(eltype(grid)), zero(eltype(grid)))
@inline temperature_term(i, j, k, grid, T, ::Nothing) = one(eltype(grid))
@inline temperature_term(i, j, k, grid, T, Q₁₀) = @inbounds Q₁₀ ^ (T[i, j, k] / 10)
@inline function growth_term(i, j, k, grid, PAR, NO₃, NH₄, T, μ₀, kPAR, kNO₃, kNH₄, ψ, Q₁₀)
    Lₗ = light_term(i, j, k, grid, PAR, kPAR)
    Lₙ = (nitrate_term(i, j, k, grid, NO₃, NH₄, kNO₃, ψ) + ammonium_term(i, j, k, grid, NH₄, kNH₄)) / 2
    return 86400 * μ₀ * Lₗ * Lₙ * temperature_term(i, j, k, grid, T, Q₁₀)
end

function limitation_outputs(ocean_model, biogeochemistry)
    underlying = biogeochemistry
    while !hasproperty(underlying, :plankton) && hasproperty(underlying, :underlying_biogeochemistry)
        underlying = underlying.underlying_biogeochemistry
    end
    plankton = underlying.plankton
    nameof(typeof(plankton.light_limitation)) == :MondoLightLimitation ||
        error("bgc_limitation.jl assumes MondoLightLimitation, the model has $(typeof(plankton.light_limitation))")
    kPAR, μ₀, Q₁₀, ψ = plankton.light_half_saturation, plankton.phytoplankton_maximum_growth_rate, plankton.temperature_coefficient, plankton.nitrate_ammonia_inhibition
    kNO₃, kNH₄ = plankton.nutrient_half_saturations.nitrate, plankton.nutrient_half_saturations.ammonia
    auxiliary = biogeochemical_auxiliary_fields(biogeochemistry)
    haskey(auxiliary, :PAR) || error("the biogeochemistry has no PAR auxiliary field")
    grid, tr, PAR = ocean_model.grid, ocean_model.tracers, auxiliary.PAR
    op(f, args...) = KernelFunctionOperation{Center, Center, Center}(f, grid, args...)
    return (; PAR, NO₃ = tr.NO₃, NH₄ = tr.NH₄,
              daylight    = op(daylight_term, PAR),
              light       = op(light_term, PAR, kPAR),
              nitrate     = op(nitrate_term, tr.NO₃, tr.NH₄, kNO₃, ψ),
              ammonium    = op(ammonium_term, tr.NH₄, kNH₄),
              temperature = op(temperature_term, tr.T, Q₁₀),
              growth      = op(growth_term, PAR, tr.NO₃, tr.NH₄, tr.T, μ₀, kPAR, kNO₃, kNH₄, ψ, Q₁₀))
end
