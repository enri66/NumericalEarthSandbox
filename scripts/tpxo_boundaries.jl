# ==================================================================
# Assemble TPXO tidal boundary forcing for a LatitudeLongitudeGrid.
#
# This is the piece that joins the two halves: it samples the atlas at
# the model's own open-boundary nodes, stores the harmonic constants
# there, and returns the closures that `GravityWaveRadiationBoundary-
# Condition` (Flather) wants.
#
# WHAT FLATHER NEEDS, AND WHY BOTH PIECES
#
#     Uᵇ = Uᵉˣᵗ ± √(gH) (ηᵇ − ηᵉˣᵗ)
#
# so each boundary needs a PAIR per boundary point: the external
# barotropic transport and the external sea level. Chapman on η radiates
# but carries no data, so Flather is the only entry point for the
# barotropic tide — which is the whole tide, for a tidal run. Omitting
# the U/V conditions and driving only the 3-D velocity gets the signal
# wiped by the split-explicit corrector, which overwrites the depth mean
# of the boundary velocity with whatever the barotropic solver says.
# (Script 9a measured this: precisely zero response.)
#
# WHERE EACH QUANTITY IS SAMPLED. This is the part that is easy to get
# half-a-cell wrong, and a half-cell error at 1/12° is ~4 km of tidal
# phase. The C-grid locations are:
#
#   west/east boundary   U at (Face, Center) → λ on the grid's λ FACE,
#                        φ at the Center of row j
#                        η at (Center, Center) → the boundary-adjacent
#                        Centre, which is the ηᵇ the condition compares to
#   south/north boundary V at (Center, Face) → φ on the grid's φ FACE,
#                        λ at the Center of column i
#                        η at the boundary-adjacent Centre
#
# TPXO's own staggering is handled inside `tpxo.jl`: each node type
# carries its own coordinate vectors, so asking for a u-transport at an
# arbitrary (λ, φ) interpolates from TPXO's u nodes, not its z nodes.
#
# TRANSPORT, NOT VELOCITY. TPXO gives depth-integrated transport
# [m² s⁻¹] and Flather acts on barotropic transport, so the atlas value
# goes in unconverted. That preserves the tidal mass flux, which is what
# sets shelf amplitude. See the discussion in tpxo.jl.
#
# CONSISTENCY. The time series is built by `reconstruct` from
# tidal_harmonics.jl using the SAME (ω, f, Θ) as the body force, so the
# boundary tide and the astronomical tide cannot drift out of phase.
# ==================================================================

using Oceananigans
using Printf

#####
##### Harmonic constants sampled on one boundary
#####

"""
    BoundaryTidalConstants

Harmonic constants for one open boundary, sized `(n_points, n_constituents)`:

- `A_n`, `G_n`: amplitude [m² s⁻¹] and Greenwich phase lag [rad] of the
  BOUNDARY-NORMAL transport
- `A_η`, `G_η`: amplitude [m] and phase lag [rad] of sea level
- `substituted`: per point, true where the atlas stencil was entirely
  land and a nearby wet value had to be substituted

`sign` is +1 or −1 and multiplies the normal transport so that the stored
value is the grid-normal component the Flather condition expects (TPXO's
`u` is west→east and `v` is south→north, which is already the grid's
positive direction on every side, so this is +1 in practice; it exists so
a caller can flip it without editing the sampling).
"""
struct BoundaryTidalConstants
    A_n :: Matrix{Float64}
    G_n :: Matrix{Float64}
    A_η :: Matrix{Float64}
    G_η :: Matrix{Float64}
    substituted :: Vector{Bool}
    sign :: Float64
end

"""
    tpxo_boundary_constants(grid, window, constituents; side, sign = 1)

Sample the atlas at the open-boundary nodes of `side`
(`:west`, `:east`, `:south`, `:north`) and return [`BoundaryTidalConstants`].

`window` is a [`TPXOWindow`] that must cover the boundary.
"""
function tpxo_boundary_constants(grid, window, constituents; side::Symbol, sign = 1.0)

    names = Symbol[Symbol(c) for c in constituents]
    nc = length(names)

    λc = λnodes(grid, Center()) ; λf = λnodes(grid, Face())
    φc = φnodes(grid, Center()) ; φf = φnodes(grid, Face())
    Nλ, Nφ = length(λc), length(φc)

    # (λ, φ) of the normal-transport node and of the adjacent η node, per point
    if side === :west
        n = Nφ
        λ_n, φ_n = fill(λf[1], n), φc
        λ_η, φ_η = fill(λc[1], n), φc
        normal = :u
    elseif side === :east
        n = Nφ
        λ_n, φ_n = fill(λf[Nλ+1], n), φc
        λ_η, φ_η = fill(λc[Nλ], n), φc
        normal = :u
    elseif side === :south
        n = Nλ
        λ_n, φ_n = λc, fill(φf[1], n)
        λ_η, φ_η = λc, fill(φc[1], n)
        normal = :v
    elseif side === :north
        n = Nλ
        λ_n, φ_n = λc, fill(φf[Nφ+1], n)
        λ_η, φ_η = λc, fill(φc[Nφ], n)
        normal = :v
    else
        throw(ArgumentError("side must be :west, :east, :south or :north, got $side"))
    end

    A_n = zeros(n, nc) ; G_n = zeros(n, nc)
    A_η = zeros(n, nc) ; G_η = zeros(n, nc)
    substituted = falses(n)

    transport = normal === :u ? tpxo_u_transport : tpxo_v_transport

    for (m, name) in enumerate(names), p in 1:n
        z_n, k_n = transport(window, name, λ_n[p], φ_n[p])
        z_η, k_η = tpxo_elevation(window, name, λ_η[p], φ_η[p])

        # A dry query is reported rather than silently zeroed: a zero transport
        # at an open boundary is a physical statement (a wall), not a missing
        # value, and conflating the two would hide a genuine sampling failure.
        if !isfinite(real(z_n)) || !isfinite(real(z_η))
            A_n[p, m] = 0.0 ; G_n[p, m] = 0.0
            A_η[p, m] = 0.0 ; G_η[p, m] = 0.0
            substituted[p] = true
        else
            a, g = amplitude_and_phase(z_n) ; A_n[p, m] = a ; G_n[p, m] = g
            a, g = amplitude_and_phase(z_η) ; A_η[p, m] = a ; G_η[p, m] = g
            (k_n == -1 || k_η == -1) && (substituted[p] = true)
        end
    end

    return BoundaryTidalConstants(A_n, G_n, A_η, G_η, substituted, Float64(sign))
end

#####
##### The Flather closures
#####

"""
    flather_pair(constants, harmonics)

Return a function with the `discrete_form` signature that
`GravityWaveRadiationBoundaryCondition` expects, yielding the
`(Uᵉˣᵗ, ηᵉˣᵗ)` pair at boundary index `p` and the clock's time.

For a west/east boundary the signature is `(j, k, grid, clock, fields)`;
for south/north it is `(i, k, grid, clock, fields)`. Either way the FIRST
argument is the along-boundary index, so one closure serves both.
"""
function flather_pair(c::BoundaryTidalConstants, harmonics::TidalHarmonics)
    p_rec = reconstruction_parameters(harmonics)
    n = size(c.A_n, 1)
    A_n, G_n, A_η, G_η, s = c.A_n, c.G_n, c.A_η, c.G_η, c.sign

    return function (p, k, grid, clock, fields)
        q = clamp(p, 1, n)
        t = clock.time
        U = zero(eltype(A_n)) ; η = zero(eltype(A_η))
        @inbounds for m in axes(A_n, 2)
            f  = p_rec.f[m] ; ω = p_rec.ω[m] ; Θ = p_rec.Θ[m]
            U += f * A_n[q, m] * cos(ω * t + Θ - G_n[q, m])
            η += f * A_η[q, m] * cos(ω * t + Θ - G_η[q, m])
        end
        return (s * U, η)
    end
end

"""
    tidal_boundary_conditions(grid, window, harmonics; sides, ramp_time = 0)

Build the `(U, V, η)` boundary conditions for a tidal run: Flather on the
barotropic transports, fed from the atlas, and Chapman on the free
surface.

`sides` is a collection of `:west`, `:east`, `:south`, `:north`; omitted
sides are left closed.

`ramp_time > 0` eases the boundary tide in as `tanh(t / ramp_time)`,
matching what `TidalBodyForce` does and for the same reason — a tide
switched on as a step rings the domain's gravity modes. The SAME ramp
should be used for both, or the two forcings are briefly inconsistent.

Returns `(; U, V, η, constants)` where `constants` is a `Dict` of the
per-side [`BoundaryTidalConstants`], kept so diagnostics can inspect what
was actually imposed.
"""
function tidal_boundary_conditions(grid, window, harmonics::TidalHarmonics;
                                   sides = (:west, :east, :south, :north),
                                   ramp_time = 0,
                                   verbose::Bool = true)

    names = harmonics.names
    constants = Dict{Symbol, BoundaryTidalConstants}()
    pairs_ = Dict{Symbol, Any}()

    for side in sides
        c = tpxo_boundary_constants(grid, window, names; side)
        constants[side] = c
        base = flather_pair(c, harmonics)
        pairs_[side] = ramp_time > 0 ?
            (p, k, g, clock, f) -> begin
                r = tanh(clock.time / ramp_time)
                U, η = base(p, k, g, clock, f)
                (r * U, r * η)
            end : base

        if verbose
            nsub = count(c.substituted)
            m1 = findfirst(==(first(names)), collect(names))
            @printf("  %-6s %4d points, %d substituted; %s |U| %.2f–%.2f m²/s, |η| %.3f–%.3f m\n",
                    side, length(c.substituted), nsub, first(names),
                    minimum(c.A_n[:, m1]), maximum(c.A_n[:, m1]),
                    minimum(c.A_η[:, m1]), maximum(c.A_η[:, m1]))
        end
    end

    U_kw = Dict{Symbol, Any}()
    V_kw = Dict{Symbol, Any}()
    η_kw = Dict{Symbol, Any}()

    for side in sides
        bc = GravityWaveRadiationBoundaryCondition(pairs_[side]; discrete_form = true)
        if side === :west || side === :east
            U_kw[side] = bc
        else
            V_kw[side] = bc
        end
        η_kw[side] = SurfaceWaveRadiationBoundaryCondition()
    end

    U = FieldBoundaryConditions(grid, (Face(), Center(), nothing); U_kw...)
    V = FieldBoundaryConditions(grid, (Center(), Face(), nothing); V_kw...)
    η = FieldBoundaryConditions(grid, (Center(), Center(), Face()); η_kw...)

    return (; U, V, η, constants)
end
