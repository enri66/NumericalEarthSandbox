# CATKE's tracer diffusivity split into its two mixing-length branches, to see which process mixes where. CATKE's
# tracer mixing length is ℓc = max(ℓ★, ℓʰ), capped at the column depth, and κc = ℓc w★ with w★ = √e (at the face):
#   ℓ★ = σ(Ri) min(Cˢ d, Cᵇ h, w★/N), the shear branch (σ the Richardson-number stability function);
#   ℓʰ, the convective branch: Cᶜ w★³/Jᵇ where the water is unstably stratified, the penetrative length Cᵉ Jᵇ/(w★ N²)
#       in the stable cell right below, and zero elsewhere or when the surface is not losing buoyancy (Jᵇ the
#       time-averaged surface buoyancy flux CATKE keeps).
# catke_diagnostics(model) returns κc, the shear and convective parts κc_shear = ℓ★ w★ and κc_convective = ℓʰ w★
# (κc is the larger of the two), N², the squared vertical shear S², the TKE e and Jᵇ.
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity, stable_length_scaleᶜᶜᶠ,
    stability_functionᶜᶜᶠ, convective_length_scaleᶜᶜᶠ, turbulent_velocityᶜᶜᶜ, shearᶜᶜᶠ
using Oceananigans.BuoyancyFormulations: ∂z_b
using Oceananigans.Operators: ℑzᵃᵃᶠ
using Oceananigans.Grids: static_column_depthᶜᶜᵃ

@inline function κc_shearᶜᶜᶠ(i, j, k, grid, closure, velocities, tracers, buoyancy)
    ml = closure.mixing_length
    σ = stability_functionᶜᶜᶠ(i, j, k, grid, closure, ml.Cᵘⁿc, ml.Cˡᵒc, ml.Cʰⁱc, velocities, tracers, buoyancy)
    ℓ = σ * stable_length_scaleᶜᶜᶠ(i, j, k, grid, closure, tracers.e, velocities, tracers, buoyancy)
    ℓ = min(static_column_depthᶜᶜᵃ(i, j, grid), ifelse(isnan(ℓ), zero(ℓ), ℓ))
    return ℓ * ℑzᵃᵃᶠ(i, j, k, grid, turbulent_velocityᶜᶜᶜ, closure, tracers.e)
end

@inline function κc_convectiveᶜᶜᶠ(i, j, k, grid, closure, velocities, tracers, buoyancy, Jᵇ)
    ml = closure.mixing_length
    ℓ = convective_length_scaleᶜᶜᶠ(i, j, k, grid, closure, ml.Cᶜc, ml.Cᵉc, ml.Cˢᵖ, velocities, tracers, buoyancy, Jᵇ)
    ℓ = min(static_column_depthᶜᶜᵃ(i, j, grid), ifelse(isnan(ℓ), zero(ℓ), ℓ))
    return ℓ * ℑzᵃᵃᶠ(i, j, k, grid, turbulent_velocityᶜᶜᶜ, closure, tracers.e)
end

function catke_diagnostics(model)
    closures = model.closure isa Tuple ? model.closure : (model.closure,)
    fields   = model.closure isa Tuple ? model.closure_fields : (model.closure_fields,)
    n = only(findall(c -> c isa CATKEVerticalDiffusivity, closures))
    closure, cf = closures[n], fields[n]
    grid, velocities, tracers, buoyancy = model.grid, model.velocities, model.tracers, model.buoyancy
    ccf(f, args...) = Field(KernelFunctionOperation{Center, Center, Face}(f, grid, args...))
    return (κc            = cf.κc,
            κc_shear      = ccf(κc_shearᶜᶜᶠ, closure, velocities, tracers, buoyancy),
            κc_convective = ccf(κc_convectiveᶜᶜᶠ, closure, velocities, tracers, buoyancy, cf.Jᵇ),
            N²            = ccf(∂z_b, buoyancy, tracers),
            S²            = ccf(shearᶜᶜᶠ, velocities.u, velocities.v),
            e             = tracers.e,
            Jᵇ            = cf.Jᵇ)
end
