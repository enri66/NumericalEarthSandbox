# ==================================================================
# Depth-dependent quadratic bottom drag.
#
# WHY A SINGLE Cᴰ IS WRONG IN PRINCIPLE. Quadratic drag τ = −Cᴰ |u| u is a
# parameterisation of the turbulent bottom boundary layer, and Cᴰ is only
# meaningful relative to the height at which the velocity is measured. Assuming a
# logarithmic layer,
#
#     Cᴰ(z_ref) = ( κ / ln(z_ref / z₀) )²
#
# with κ = 0.41 and z₀ the roughness length. In a model the velocity is the lowest
# WET CELL's, so z_ref is half that cell's thickness — which varies enormously over
# a shelf-to-abyss domain. On the MAB grid the lowest wet cell runs from ~1 m on
# the inner shelf to ~680 m in the deep ocean, so a single Cᴰ is simultaneously too
# small inshore and too large offshore. This gives each column the Cᴰ its own
# resolution implies.
#
# THE STABILITY PROBLEM THIS CREATES, AND THE FIX. Explicit quadratic drag is
# stable only while Cᴰ |u| Δt / Δz_bot < 1, and the log law asks for MORE drag
# exactly where Δz_bot is SMALLEST — the two requirements point opposite ways.
# Measured on this grid: at Cᴰ = 0.012 the worst inner-shelf cell reaches a ratio
# of 2.4 and the run goes NaN. So the explicit scheme cannot deliver the drag the
# log law wants inshore, which is precisely where it matters.
#
# The fix is to treat the drag SEMI-IMPLICITLY, which Oceananigans supports
# directly through `IMEXFluxBoundaryCondition(Fₑ, λ)`, whose flux is affine in the
# boundary-cell value:
#
#     J(φᵦ) = Fₑ + λ φᵦ
#
# Writing the drag as Fₑ = 0 and λ = −Cᴰ |u|ⁿ gives J = −Cᴰ |u|ⁿ uⁿ⁺¹: the speed
# lags but the damped velocity is implicit, which is the standard semi-implicit
# treatment and is unconditionally stable. The λ term is folded into the vertical
# tridiagonal solve, so the model MUST carry a vertically implicit closure — see
# `implicit_drag_requires_vertical_closure` below.
#
# CAVEAT WORTH KEEPING. The log law is a bottom-boundary-layer result and stops
# being physical once z_ref leaves that layer. With a 680 m bottom cell in the
# abyss, z_ref = 340 m is far outside it and the formula is being extrapolated.
# It returns a small Cᴰ there, which is harmless because abyssal tidal velocities
# are ~1 cm/s and the drag is negligible either way — but the number should not be
# read as physically derived. `Cᴰ_bounds` clamps the range to keep it honest.
# ==================================================================

using Oceananigans
using Oceananigans.Grids: static_column_depthᶜᶜᵃ
using Oceananigans.Operators: Δzᶜᶜᶜ, ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
using Oceananigans.ImmersedBoundaries: immersed_cell
using Oceananigans.BoundaryConditions: IMEXFluxBoundaryCondition
using Printf, Statistics

"""
    bottom_cell_thickness(grid)

Thickness [m] of the lowest WET cell in each column, and the index of that cell.
`NaN` / `0` where the column is entirely dry.
"""
function bottom_cell_thickness(grid)
    Nx, Ny, Nz = size(grid)
    Δ = fill(NaN, Nx, Ny)
    k_bot = zeros(Int, Nx, Ny)
    # NOTE the nesting. `for i in …, j in …, k in …` is a SINGLE fused loop in
    # Julia, so a `break` inside it exits ALL of them — which silently left every
    # column but the first as NaN. The k loop must be separate for `break` to mean
    # "found this column's lowest wet cell".
    for i in 1:Nx, j in 1:Ny
        for k in 1:Nz
            if !immersed_cell(i, j, k, grid)
                Δ[i, j] = Δzᶜᶜᶜ(i, j, k, grid)
                k_bot[i, j] = k
                break
            end
        end
    end
    return Δ, k_bot
end

"""
    log_law_drag_coefficient(grid; z₀ = 0.003, κ = 0.41, Cᴰ_bounds = (1e-3, 1e-2))

A field of quadratic drag coefficients, `Cᴰ = (κ / ln(z_ref/z₀))²`, with `z_ref`
half the thickness of each column's lowest wet cell. Returns a
`Field{Center, Center, Nothing}(grid)`, clamped to `Cᴰ_bounds`.

`z₀ = 3 mm` is a common shelf value; 1 mm (smooth sand) to 10 mm (rough/rippled)
is the usual range and moves `Cᴰ` by roughly ±40%.

A `Field`, not a plain array: it carries its own architecture, so it is placed
correctly on GPU and multi-GPU without the boundary-condition machinery needing
to know anything about it.
"""
function log_law_drag_coefficient(grid; z₀ = 0.003, κ = 0.41, Cᴰ_bounds = (1e-3, 1e-2))
    Δ, _ = bottom_cell_thickness(grid)
    Nx, Ny = size(Δ)
    values = fill(Cᴰ_bounds[1], Nx, Ny)
    for i in 1:Nx, j in 1:Ny
        isfinite(Δ[i, j]) && Δ[i, j] > 0 || continue
        z_ref = Δ[i, j] / 2
        z_ref > z₀ || continue                       # below the roughness scale: use the floor
        c = (κ / log(z_ref / z₀))^2
        values[i, j] = clamp(c, Cᴰ_bounds[1], Cᴰ_bounds[2])
    end
    Cᴰ = Field{Center, Center, Nothing}(grid)
    interior(Cᴰ, :, :, 1) .= values
    return Cᴰ
end

"""
    report_drag_field(Cᴰ, grid; Δt, u_scale = 1.0)

Print the distribution of a variable drag field by depth class, alongside the
EXPLICIT stability ratio `Cᴰ u Δt / Δz_bot` it would imply — the number that
decides whether the drag has to be implicit.
"""
function report_drag_field(Cᴰ, grid; Δt, u_scale = 1.0)
    Cᴰ_values = Array(interior(Cᴰ, :, :, 1))
    Δ, _ = bottom_cell_thickness(grid)
    Nx, Ny = size(Δ)
    H = [static_column_depthᶜᶜᵃ(i, j, grid) for i in 1:Nx, j in 1:Ny]
    ratio = @. Cᴰ_values * u_scale * Δt / Δ

    @printf("  drag field: Cᴰ ∈ [%.2e, %.2e], explicit stability ratio at |u| = %.1f m/s, Δt = %.0f s\n",
            minimum(skipmissing(filter(isfinite, Cᴰ_values))), maximum(filter(isfinite, Cᴰ_values)), u_scale, Δt)
    println("  class            n    median Δz_bot   median Cᴰ    max Cᴰ   median ratio    MAX ratio")
    for (nm, lo, hi) in (("shelf <50 m", 0.0, 50.0), ("shelf 50–200", 50.0, 200.0),
                         ("slope 200–1000", 200.0, 1000.0), ("deep >1000", 1000.0, 1e9))
        sel = @. isfinite(Δ) & (Δ > 0) & (H > lo) & (H <= hi)
        count(sel) < 10 && continue
        @printf("  %-15s %5d %14.2f %12.2e %9.2e %14.2f %12.2f\n",
                nm, count(sel), median(Δ[sel]), median(Cᴰ_values[sel]), maximum(Cᴰ_values[sel]),
                median(ratio[sel]), maximum(ratio[sel]))
    end
    unstable = count(@. isfinite(ratio) & (ratio > 1))
    @printf("  cells where EXPLICIT drag would be unstable (ratio > 1): %d\n", unstable)
    return nothing
end

#####
##### Drag kernels taking a Field coefficient
#####
# Speed at the u and v points, matching NumericalEarth.Oceans' own definitions.
@inline _ϕ²(i, j, k, grid, ϕ) = @inbounds ϕ[i, j, k]^2
@inline _spᶠᶜᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, _ϕ², Φ.v))
@inline _spᶜᶠᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, _ϕ², Φ.u))

# Cᴰ lives at (Center, Center) while u is at (Face, Center) and v at (Center, Face),
# so indexing it directly is a half-cell approximation. Cᴰ varies on the scale of
# the bathymetry, far longer than a cell, so this is well below the uncertainty in
# z₀ — and interpolating it would need halo-filled Cᴰ for no measurable gain.

# --- EXPLICIT form (stability-limited; kept for comparison) ---
@inline u_variable_bottom_drag(i, j, grid, c, Φ, μ) =
    @inbounds -μ[i, j, 1] * Φ.u[i, j, 1] * _spᶠᶜᶜ(i, j, 1, grid, Φ)
@inline v_variable_bottom_drag(i, j, grid, c, Φ, μ) =
    @inbounds -μ[i, j, 1] * Φ.v[i, j, 1] * _spᶜᶠᶜ(i, j, 1, grid, Φ)
@inline u_variable_immersed_drag(i, j, k, grid, c, Φ, μ) =
    @inbounds -μ[i, j, 1] * Φ.u[i, j, k] * _spᶠᶜᶜ(i, j, k, grid, Φ)
@inline v_variable_immersed_drag(i, j, k, grid, c, Φ, μ) =
    @inbounds -μ[i, j, 1] * Φ.v[i, j, k] * _spᶜᶠᶜ(i, j, k, grid, Φ)

# --- SEMI-IMPLICIT form: λ = −Cᴰ|u|ⁿ, so the flux is −Cᴰ|u|ⁿ uⁿ⁺¹ ---
@inline u_drag_coefficient(i, j, grid, c, Φ, μ) =
    @inbounds -μ[i, j, 1] * _spᶠᶜᶜ(i, j, 1, grid, Φ)
@inline v_drag_coefficient(i, j, grid, c, Φ, μ) =
    @inbounds -μ[i, j, 1] * _spᶜᶠᶜ(i, j, 1, grid, Φ)
@inline u_immersed_drag_coefficient(i, j, k, grid, c, Φ, μ) =
    @inbounds -μ[i, j, 1] * _spᶠᶜᶜ(i, j, k, grid, Φ)
@inline v_immersed_drag_coefficient(i, j, k, grid, c, Φ, μ) =
    @inbounds -μ[i, j, 1] * _spᶜᶠᶜ(i, j, k, grid, Φ)

"""
    variable_drag_boundary_conditions(Cᴰ; implicit = true)

Return `(u_bottom, u_immersed, v_bottom, v_immersed)` boundary conditions applying
quadratic drag with the coefficient field `Cᴰ`, a `Field{Center, Center, Nothing}(grid)`.

`implicit = true` (the default, and the reason this file exists) uses
`IMEXFluxBoundaryCondition`, so the drag is unconditionally stable and the log
law's large inshore `Cᴰ` is actually attainable. It REQUIRES the model to carry a
vertically implicit closure, because the implicit term is solved by the vertical
tridiagonal solver; without one it is silently inert. Pass
`closure = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν = …)`
alongside whatever else you use.

`implicit = false` gives the explicit form, for comparison — subject to
`Cᴰ |u| Δt / Δz_bot < 1`.
"""
function variable_drag_boundary_conditions(Cᴰ; implicit::Bool = true)
    if implicit
        u_bot = IMEXFluxBoundaryCondition(0, u_drag_coefficient;
                                          discrete_form = true, parameters = Cᴰ)
        v_bot = IMEXFluxBoundaryCondition(0, v_drag_coefficient;
                                          discrete_form = true, parameters = Cᴰ)
        u_imm = ImmersedBoundaryCondition(bottom =
            IMEXFluxBoundaryCondition(0, u_immersed_drag_coefficient;
                                      discrete_form = true, parameters = Cᴰ))
        v_imm = ImmersedBoundaryCondition(bottom =
            IMEXFluxBoundaryCondition(0, v_immersed_drag_coefficient;
                                      discrete_form = true, parameters = Cᴰ))
    else
        u_bot = FluxBoundaryCondition(u_variable_bottom_drag; discrete_form = true, parameters = Cᴰ)
        v_bot = FluxBoundaryCondition(v_variable_bottom_drag; discrete_form = true, parameters = Cᴰ)
        u_imm = ImmersedBoundaryCondition(bottom =
            FluxBoundaryCondition(u_variable_immersed_drag; discrete_form = true, parameters = Cᴰ))
        v_imm = ImmersedBoundaryCondition(bottom =
            FluxBoundaryCondition(v_variable_immersed_drag; discrete_form = true, parameters = Cᴰ))
    end
    return (u_bottom = u_bot, u_immersed = u_imm, v_bottom = v_bot, v_immersed = v_imm)
end
