# ==================================================================
# Tidal harmonics: the astronomy shared by the BODY FORCE and the
# BOUNDARY FORCING.
#
# A regional tidal simulation is driven two ways at once:
#
#   1. a BODY FORCE  — the astronomical (equilibrium) tide, entering the
#      momentum equations as ∂u/∂t = … − g∇(η − η_eq), i.e. +g∇η_eq;
#   2. a BOUNDARY FORCE — sea level and barotropic transport prescribed
#      at the open boundaries from a tidal atlas (TPXO), fed to the
#      Flather/Chapman conditions.
#
# The two are only *consistent* if they share three things: ONE
# constituent list, ONE reference epoch, and ONE set of nodal factors
# f(t) and equilibrium arguments V(t). That shared core is this file.
# Everything physical lives in `astronomical_arguments`, which returns
# (ω, f, Θ) per constituent; the equilibrium tide and the TPXO
# reconstruction are then two different ways of weighting the SAME
# (ω, f, Θ). Get the sharing wrong and the body tide and the boundary
# tide drift out of phase with each other — the classic symptom being a
# standing wave that refuses to match the cotidal chart no matter how
# the amplitudes are tuned.
#
# This is what ROMS does with POT_TIDES + SSH_TIDES/UV_TIDES, and MOM6
# with MOM_tidal_forcing + OBC_TIDE_* (in particular their
# OBC_TIDE_ADD_EQ_PHASE / OBC_TIDE_ADD_NODAL switches, which exist
# precisely to enforce this consistency).
#
# CONSTITUENTS (10) — the 8 primary astronomical tides plus the two
# dominant long-period ones:
#   semidiurnal   M2  S2  N2  K2
#   diurnal       K1  O1  P1  Q1
#   long-period   Mf  Mm
# The compound/overtides a TPXO atlas also carries (M4, MS4, MN4, 2N2,
# S1) are deliberately EXCLUDED: on a shelf they are generated
# internally by nonlinear advection, bottom drag and bathymetry, so
# forcing them at the boundary as well would double-count them.
#
# SOURCES (all published; nothing transliterated)
#   - Schureman, P. (1958) "Manual of Harmonic Analysis and Prediction
#     of Tides", USC&GS Special Publication 98.
#   - Kowalik, Z. & J. Luick (2019) "Modern Theory and Practice of Tide
#     Analysis and Tidal Power": Eq. I.71 (astronomical longitudes),
#     Table I.4 (equilibrium phase arguments), Table I.6 (nodal f, u).
#   - Cartwright, D. E. & R. J. Tayler (1971) / Cartwright & Edden
#     (1973) — the equilibrium amplitudes.
#   - Love-number factors (1 + k − h): Wahr (1981) elastic Earth.
#   Cross-checked against MOM6 `MOM_tidal_forcing.F90`, which uses the
#   same tables for the same 10 constituents.
# ==================================================================

using Dates
using Printf
using LinearAlgebra

#####
##### The constituent table
#####

"""
    TidalConstituentData(ω, A, love, s)

Static properties of one tidal constituent.

- `ω`: angular frequency [rad s⁻¹]
- `A`: Cartwright–Tayler equilibrium amplitude [m]
- `love`: Love-number factor (1 + k − h) for the elastic Earth [-]
- `s`: longitude multiplier, which also identifies the species —
  `2` semidiurnal, `1` diurnal, `0` long-period.

`s` is stored rather than a species label because it is what the
formulae actually use: the equilibrium tide of a constituent varies as
`cos(ωt + Θ + s·λ)`, and its latitude structure follows from `s`.
"""
struct TidalConstituentData
    ω    :: Float64
    A    :: Float64
    love :: Float64
    s    :: Int
end

const TIDAL_CONSTITUENTS = (
    # semidiurnal  (s = 2)
    M2 = TidalConstituentData(1.4051890e-4, 0.242334, 0.693, 2),
    S2 = TidalConstituentData(1.4544410e-4, 0.112743, 0.693, 2),
    N2 = TidalConstituentData(1.3787970e-4, 0.046397, 0.693, 2),
    K2 = TidalConstituentData(1.4584234e-4, 0.030684, 0.693, 2),
    # diurnal      (s = 1)
    K1 = TidalConstituentData(0.7292117e-4, 0.141565, 0.736, 1),
    O1 = TidalConstituentData(0.6759774e-4, 0.100661, 0.695, 1),
    P1 = TidalConstituentData(0.7252295e-4, 0.046848, 0.706, 1),
    Q1 = TidalConstituentData(0.6495854e-4, 0.019273, 0.695, 1),
    # long-period  (s = 0)
    Mf = TidalConstituentData(0.053234e-4,  0.042041, 0.693, 0),
    Mm = TidalConstituentData(0.026392e-4,  0.022191, 0.693, 0),
)

"All ten constituents, in the conventional order."
const ALL_CONSTITUENTS = (:M2, :S2, :N2, :K2, :K1, :O1, :P1, :Q1, :Mf, :Mm)

"Period of a constituent [s]."
constituent_period(name::Symbol) = 2π / TIDAL_CONSTITUENTS[name].ω

#####
##### Astronomical longitudes  (Kowalik & Luick Eq. I.71, after Schureman 1958)
#####

"""
    astronomical_longitudes(t::DateTime)

The four mean longitudes the tidal equilibrium arguments are built from,
in radians:

- `s`: mean longitude of the Moon
- `h`: mean longitude of the Sun
- `p`: mean longitude of lunar perigee
- `N`: longitude of the Moon's ascending node

Time is measured in Julian centuries from 1900-01-01 00:00 UTC. Valid
for the Gregorian calendar; the polynomial fits are good to well under a
degree over the 20th–21st centuries, which is far finer than the tidal
phases need.
"""
function astronomical_longitudes(t::DateTime)
    D = Dates.value(t - DateTime(1900, 1, 1)) / 86_400_000   # days since 1900-01-01
    T = D / 36_525                                           # Julian centuries

    s = mod(277.0248 + 481267.8906 * T + 0.0011    * T^2, 360)
    h = mod(280.1895 +  36000.7689 * T + 3.0310e-4 * T^2, 360)
    p = mod(334.3853 +   4069.0340 * T - 0.0103    * T^2, 360)
    N = mod(259.1568 -   1934.1420 * T + 0.0021    * T^2, 360)

    return (s = deg2rad(s), h = deg2rad(h), p = deg2rad(p), N = deg2rad(N))
end

"""
    equilibrium_phase(name, L)

The equilibrium phase argument `V₀` [rad] of constituent `name`, given
the astronomical longitudes `L` from [`astronomical_longitudes`].
Kowalik & Luick Table I.4.
"""
function equilibrium_phase(name::Symbol, L)
    s, h, p = L.s, L.h, L.p

    name === :M2 && return 2 * (h - s)
    name === :S2 && return zero(h)
    name === :N2 && return -3s + 2h + p
    name === :K2 && return 2h
    name === :K1 && return  h + π/2
    name === :O1 && return -2s + h - π/2
    name === :P1 && return -h - π/2
    name === :Q1 && return -3s + h + p - π/2
    name === :Mf && return 2s
    name === :Mm && return s - p

    throw(ArgumentError("unknown tidal constituent $name"))
end

"""
    nodal_factors(name, N)

Amplitude factor `f` [-] and phase correction `u` [rad] describing the
modulation of constituent `name` by the 18.6-year nodal cycle, given the
longitude `N` of the ascending node. Kowalik & Luick Table I.6.

`f` departs from 1 by up to 41% (Mf) and 29% (K2), so it is not
optional: a run that ignores it can be badly wrong in amplitude for the
nodal-sensitive constituents even when every phase is right.
"""
function nodal_factors(name::Symbol, N)
    # purely solar constituents have no lunar-node modulation
    (name === :S2 || name === :P1) && return (one(N), zero(N))

    cosN, sinN = cos(N), sin(N)

    name === :M2 && return (1 - 0.037 * cosN, deg2rad( -2.1) * sinN)
    name === :N2 && return (1 - 0.037 * cosN, deg2rad( -2.1) * sinN)
    name === :K2 && return (1.024 + 0.286 * cosN, deg2rad(-17.7) * sinN)
    name === :K1 && return (1.006 + 0.115 * cosN, deg2rad( -8.9) * sinN)
    name === :O1 && return (1.009 + 0.187 * cosN, deg2rad( 10.8) * sinN)
    name === :Q1 && return (1.009 + 0.187 * cosN, deg2rad( 10.8) * sinN)
    name === :Mf && return (1.043 + 0.414 * cosN, deg2rad(-23.7) * sinN)
    name === :Mm && return (1 - 0.130 * cosN, zero(N))

    throw(ArgumentError("unknown tidal constituent $name"))
end

#####
##### The shared harmonic state
#####

"""
    TidalHarmonics

The resolved astronomy for a chosen constituent set at a chosen epoch —
the single object both the body force and the boundary forcing read, so
that they cannot disagree.

Fields, all `NTuple`s over the `N` constituents:

- `ω`: angular frequency [rad s⁻¹]
- `f`: nodal amplitude factor [-]
- `Θ`: equilibrium argument plus nodal phase, `V₀ + u` [rad], evaluated
  at `reference_time`
- `Aα`: Cartwright amplitude × Love factor [m] — used by the
  equilibrium tide only; the TPXO reconstruction supplies its own
  amplitudes
- `s`: longitude multiplier / species

Model time `t` in every routine below is **seconds since
`reference_time`**, matching an Oceananigans clock started at that
date. The epoch is absorbed into `Θ`, so a constituent's contribution
is always `f · A · cos(ω t + Θ + …)`.
"""
struct TidalHarmonics{N, FT}
    names          :: NTuple{N, Symbol}
    ω              :: NTuple{N, FT}
    f              :: NTuple{N, FT}
    Θ              :: NTuple{N, FT}
    Aα             :: NTuple{N, FT}
    s              :: NTuple{N, Int}
    reference_time :: DateTime
end

"""
    TidalHarmonics(reference_time; constituents = ALL_CONSTITUENTS,
                   nodal_time = reference_time, FT = Float64)

Resolve the astronomy for `constituents` at `reference_time`.

`nodal_time` sets where the 18.6-year nodal factors are evaluated. It
defaults to `reference_time`, which is right for runs up to ~a year; for
a multi-year run pass the midpoint so the error is split, or rebuild the
harmonics periodically.
"""
function TidalHarmonics(reference_time::DateTime;
                        constituents = ALL_CONSTITUENTS,
                        nodal_time::DateTime = reference_time,
                        FT = Float64)

    names = Tuple(Symbol(c) for c in constituents)

    for name in names
        haskey(TIDAL_CONSTITUENTS, name) ||
            throw(ArgumentError("unknown tidal constituent $name; known: $(ALL_CONSTITUENTS)"))
    end

    L  = astronomical_longitudes(reference_time)
    Ln = astronomical_longitudes(nodal_time)

    ω  = Tuple(FT(TIDAL_CONSTITUENTS[n].ω) for n in names)
    s  = Tuple(TIDAL_CONSTITUENTS[n].s     for n in names)
    Aα = Tuple(FT(TIDAL_CONSTITUENTS[n].A * TIDAL_CONSTITUENTS[n].love) for n in names)

    fu = Tuple(nodal_factors(n, Ln.N) for n in names)
    f  = Tuple(FT(x[1]) for x in fu)
    Θ  = Tuple(FT(equilibrium_phase(n, L) + x[2]) for (n, x) in zip(names, fu))

    return TidalHarmonics(names, ω, f, Θ, Aα, s, reference_time)
end

Base.length(::TidalHarmonics{N}) where N = N

function Base.show(io::IO, h::TidalHarmonics{N, FT}) where {N, FT}
    print(io, "TidalHarmonics{$N, $FT} at $(h.reference_time)\n")
    print(io, "┌──────┬────────────┬─────────┬──────────┬──────────┐\n")
    print(io, "│ name │ period (h) │  f      │  Θ (°)   │  A·α (m) │\n")
    print(io, "├──────┼────────────┼─────────┼──────────┼──────────┤\n")
    for n in 1:N
        @printf(io, "│ %-4s │ %10.4f │ %7.4f │ %8.3f │ %8.5f │\n",
                h.names[n], 2π / h.ω[n] / 3600, h.f[n],
                mod(rad2deg(h.Θ[n]), 360), h.Aα[n])
    end
    print(io, "└──────┴────────────┴─────────┴──────────┴──────────┘")
end

#####
##### GPU-safe parameter bundles
#####
# `Symbol` and `DateTime` fields are not wanted inside a kernel, so the
# forcing functions take these stripped, isbits NamedTuples instead.

"Parameters for [`equilibrium_tide`] and its gradient."
equilibrium_tide_parameters(h::TidalHarmonics) =
    (ω = h.ω, f = h.f, Θ = h.Θ, Aα = h.Aα, s = h.s,
     t₀ = Dates.value(h.reference_time))

"Parameters for [`reconstruct`] — no `Aα`, amplitudes come from the atlas."
reconstruction_parameters(h::TidalHarmonics) =
    (ω = h.ω, f = h.f, Θ = h.Θ, t₀ = Dates.value(h.reference_time))

#####
##### Model time
#####
# Every routine below wants "seconds since the harmonics' reference time".
#
# `ocean_simulation` defaults to `Clock(grid)`, a NUMERIC clock starting at
# time = 0, so a model started on the reference date reports exactly that and
# the `Real` method is a no-op. But `ocean_simulation` also accepts a
# `DateTime`-based clock, which is what a coupled run may use — and then `t`
# arrives as a `DateTime`. Dispatching here means the tidal forcing is correct
# under either clock instead of throwing (or, worse, being reinterpreted) the
# day someone switches. `t₀` is the reference time in milliseconds.

@inline elapsed_seconds(t::Real, p) = t
@inline elapsed_seconds(t::Dates.DateTime, p) = (Dates.value(t) - p.t₀) / 1000

#####
##### The equilibrium tide
#####
# Species geometry. For longitude multiplier s the equilibrium tide of a
# constituent is
#
#     η = A·α·f · G_s(φ) · cos(ω t + Θ + s·λ)
#
# with
#     s = 2 (semidiurnal)   G = cos²φ
#     s = 1 (diurnal)       G = sin 2φ
#     s = 0 (long-period)   G = ½ − (3/2) sin²φ
#
# This is the compact form of the cos/sin structure pair MOM6 tabulates:
# cosΘ·cos sλ − sinΘ·sin sλ = cos(Θ + sλ). Writing it as one cosine
# makes the analytic λ- and φ-derivatives below trivial, which is what
# lets the body force be an ordinary continuous-form `Forcing` with no
# auxiliary field and no callback.

@inline function species_structure(s::Int, φ)
    s == 2 && return cos(φ)^2
    s == 1 && return sin(2φ)
    return 1//2 - 3//2 * sin(φ)^2
end

"d/dφ of [`species_structure`]."
@inline function species_structure_derivative(s::Int, φ)
    s == 2 && return -sin(2φ)
    s == 1 && return 2cos(2φ)
    return -3//2 * sin(2φ)
end

"""
    equilibrium_tide(λ, φ, t, p)

Equilibrium tide elevation `η_eq` [m] at longitude `λ` and latitude `φ`
(both **degrees**, as Oceananigans reports `LatitudeLongitudeGrid`
nodes) and model time `t` [s since the harmonics' reference time).

`p` is [`equilibrium_tide_parameters`]. The associated tide-generating
potential is `Φ = g · η_eq`.
"""
@inline function equilibrium_tide(λ, φ, t, p)
    λr, φr = deg2rad(λ), deg2rad(φ)
    τ = elapsed_seconds(t, p)
    η = zero(λr)
    @inbounds for n in eachindex(p.ω)
        η += p.Aα[n] * p.f[n] * species_structure(p.s[n], φr) *
             cos(p.ω[n] * τ + p.Θ[n] + p.s[n] * λr)
    end
    return η
end

"""
    equilibrium_tide_gradient(λ, φ, t, p)

`(∂η_eq/∂λ, ∂η_eq/∂φ)` [m rad⁻¹], differentiated analytically. Angles
in degrees as above; the derivatives are with respect to *radians*,
which is what the metric terms below expect.
"""
@inline function equilibrium_tide_gradient(λ, φ, t, p)
    λr, φr = deg2rad(λ), deg2rad(φ)
    τ = elapsed_seconds(t, p)
    ∂λ = zero(λr)
    ∂φ = zero(λr)
    @inbounds for n in eachindex(p.ω)
        Afn = p.Aα[n] * p.f[n]
        Θn  = p.ω[n] * τ + p.Θ[n] + p.s[n] * λr
        ∂λ -= Afn * species_structure(p.s[n], φr) * p.s[n] * sin(Θn)
        ∂φ += Afn * species_structure_derivative(p.s[n], φr) * cos(Θn)
    end
    return (∂λ, ∂φ)
end

#####
##### Harmonic reconstruction — the TPXO boundary side
#####

"""
    reconstruct(amplitude, phase_lag, t, p)

Predict a tidal quantity from atlas harmonic constants:

    x(t) = Σₖ fₖ Aₖ cos(ωₖ t + Θₖ − Gₖ)

where `amplitude[k]` = `Aₖ` and `phase_lag[k]` = `Gₖ` [rad] are the
Greenwich amplitude and phase lag from the atlas (TPXO), and `(ω, f, Θ)`
come from `p` = [`reconstruction_parameters`] — *the same* `(ω, f, Θ)`
the equilibrium tide uses. That shared triple is the whole consistency
guarantee.

Works for elevation (m), transport (m² s⁻¹), or velocity (m s⁻¹) alike;
it is the harmonic sum that is generic, not the variable.
"""
@inline function reconstruct(amplitude, phase_lag, t, p)
    τ = elapsed_seconds(t, p)
    x = zero(eltype(amplitude))
    @inbounds for n in eachindex(p.ω)
        x += p.f[n] * amplitude[n] * cos(p.ω[n] * τ + p.Θ[n] - phase_lag[n])
    end
    return x
end

#####
##### Harmonic analysis — the inverse, for verification and de-tiding
#####

"""
    harmonic_analysis(t, x, h::TidalHarmonics; mean = true)

Least-squares fit of the constituents in `h` to the time series `x`
sampled at times `t` [s since `h.reference_time`]. Returns
`(; amplitude, phase_lag, mean_value, residual)` with `amplitude` in the
units of `x`, `phase_lag` the Greenwich lag in **degrees** on `[0, 360)`,
and `residual` the RMS misfit.

The returned constants are directly comparable with a tidal atlas or
with published tide-gauge harmonic constants, because the fit undoes the
same `f` and `Θ` that [`reconstruct`] applies:

    x(t) = x̄ + Σₖ [aₖ cos ωₖt + bₖ sin ωₖt]
    Aₖ = √(aₖ² + bₖ²) / fₖ ,    Gₖ = Θₖ + atan(bₖ, aₖ)

`t` must span enough time to separate the requested constituents — the
Rayleigh criterion, `Δt_record > 2π/|ωᵢ − ωⱼ|`. The closest pair here is
K2–S2, needing ~182 days; S2–M2 needs ~14.8 days; M2 alone needs only a
couple of days. Asking for constituents the record cannot resolve
produces a badly conditioned fit rather than an error, so
`condition_number` is returned for inspection.
"""
function harmonic_analysis(t::AbstractVector, x::AbstractVector, h::TidalHarmonics{N};
                           mean::Bool = true) where N

    length(t) == length(x) || throw(DimensionMismatch("t and x must have equal length"))

    ncol = 2N + (mean ? 1 : 0)
    length(t) >= ncol ||
        throw(ArgumentError("need at least $ncol samples to fit $N constituents"))

    M = Matrix{Float64}(undef, length(t), ncol)
    for n in 1:N
        @. M[:, 2n-1] = cos(h.ω[n] * t)
        @. M[:, 2n  ] = sin(h.ω[n] * t)
    end
    mean && (M[:, end] .= 1)

    coefficients = M \ collect(Float64, x)
    residual = sqrt(sum(abs2, M * coefficients .- x) / length(x))

    amplitude = zeros(N)
    phase_lag = zeros(N)
    for n in 1:N
        a, b = coefficients[2n-1], coefficients[2n]
        amplitude[n] = sqrt(a^2 + b^2) / h.f[n]
        phase_lag[n] = mod(rad2deg(h.Θ[n] + atan(b, a)), 360)
    end

    return (; amplitude, phase_lag,
              mean_value = mean ? coefficients[end] : 0.0,
              residual,
              condition_number = cond(M),
              names = h.names)
end

#####
##### Body force on a LatitudeLongitudeGrid
#####
# THE SIGN. Write the momentum equation with the tidal term folded into
# the pressure gradient:
#
#     ∂u/∂t = … − g ∇(η − η_eq)
#
# which makes the equilibrium (static) solution η = η_eq manifest — that
# is the defining property of η_eq, so this form is the definition rather
# than a convention. Hence the tidal body force is
#
#     F = + g ∇η_eq
#
# and on a sphere, with λ, φ in radians,
#
#     Fᵘ = + (g / (R cos φ)) ∂η_eq/∂λ ,    Fᵛ = + (g / R) ∂η_eq/∂φ
#
# both [m s⁻²]. Equivalently the tide-generating POTENTIAL is
# Ω = −g·η_eq (water piles up in the potential well), and the force is
# −∇Ω. It is easy to get this backwards by calling g·η_eq "the
# potential" and then writing −∇ of it; the closed-basin test in
# scripts/verify_tidal_body_force.jl exists to catch exactly that, and
# did: it returned a least-squares gain of −1.015 against the expected
# +1. MOM6 agrees — `MOM_PressureForce_FV.F90` applies
# `za ← za − g·(e_tidal_eq + e_tidal_sal)` to the geopotential whose
# gradient enters the momentum equation as −∇za, which is +g∇η_eq.
#
# η_eq is z-independent, so this is a barotropic body force applied
# uniformly over the column — which is what the astronomical tide is.
#
# SELF-ATTRACTION AND LOADING is included only as the scalar
# approximation η_sal ≈ β·η_eq, folded into `scalar_sal` as a factor
# (1 + β) on the whole potential. The honest version is
# η_sal ∝ the *total* ocean surface load, which is a global integral and
# cannot be computed in a regional domain at all. For a MAB run this
# matters little: the tide here is overwhelmingly set by what arrives
# through the open boundaries, and the TPXO boundary data already
# contains the real SAL implicitly. `scalar_sal = 0` (the default) is the
# defensible choice for a regional run; the knob exists so that a
# large-domain configuration can switch it on.

"""
    TidalBodyForce(harmonics; gravitational_acceleration = 9.80665,
                   radius = 6371e3, scalar_sal = 0, ramp_time = 0,
                   FT = Float64)

Build `(u = Fᵘ, v = Fᵛ)` continuous-form `Forcing`s implementing the
equilibrium tidal body force on a `LatitudeLongitudeGrid`, ready to pass
to a model's `forcing` keyword (or to merge with other forcings).

No auxiliary field and no callback: the harmonic sum and its analytic
gradient are evaluated inline per cell, ~10 cosines for the full
constituent set.

`ramp_time > 0` multiplies the force by `tanh(t / ramp_time)`, easing it
in from rest instead of switching it on as a step. Starting a basin
impulsively rings its gravest gravity mode, and with no closure there is
nothing to damp that seiche — it then contaminates any harmonic analysis
of the response. One tidal period is a good ramp. (Script 9a ramps its
wind stress over a day for the same reason.) Default `0` is no ramp.
"""
function TidalBodyForce(harmonics::TidalHarmonics;
                        gravitational_acceleration = 9.80665,
                        radius = 6371e3,
                        scalar_sal = 0,
                        ramp_time = 0,
                        FT = Float64)

    p = equilibrium_tide_parameters(harmonics)

    parameters = (; p.ω, p.f, p.Θ, p.s, p.t₀,
                  Aα = map(A -> FT(A * (1 + scalar_sal)), p.Aα),
                  g  = FT(gravitational_acceleration),
                  R  = FT(radius),
                  ramp = FT(ramp_time))

    u_forcing = Forcing(tidal_body_force_u; parameters)
    v_forcing = Forcing(tidal_body_force_v; parameters)

    return (u = u_forcing, v = v_forcing)
end

# tanh ramp, guarding the 0/0 at t = 0 when no ramp is requested.
@inline function tidal_ramp(τ, ramp)
    ramp <= 0 && return one(τ)
    return tanh(τ / ramp)
end

@inline function tidal_body_force_u(λ, φ, z, t, p)
    ∂λ, _ = equilibrium_tide_gradient(λ, φ, t, p)
    τ = elapsed_seconds(t, p)
    return p.g * ∂λ / (p.R * cos(deg2rad(φ))) * tidal_ramp(τ, p.ramp)
end

@inline function tidal_body_force_v(λ, φ, z, t, p)
    _, ∂φ = equilibrium_tide_gradient(λ, φ, t, p)
    τ = elapsed_seconds(t, p)
    return p.g * ∂φ / p.R * tidal_ramp(τ, p.ramp)
end
