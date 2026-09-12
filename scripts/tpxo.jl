# ==================================================================
# Read the TPXO10-atlas-v2 tidal atlas, and interpolate its harmonic
# constants to arbitrary points.
#
# This supplies the BOUNDARY half of the tidal forcing: complex
# amplitudes of sea-surface elevation and of depth-integrated transport,
# which `tidal_harmonics.jl` then turns into a time series via
# `reconstruct` using the SAME (ω, f, Θ) the body force uses.
#
# ------------------------------------------------------------------
# THE SCHEMA, every item verified against the files, not assumed
# ------------------------------------------------------------------
#
# Layout: one file per constituent per variable, plus one grid file.
#   h_<con>_tpxo10_atlas_30_v2.nc   elevation
#   u_<con>_tpxo10_atlas_30_v2.nc   transports (BOTH components)
#   grid_tpxo10atlas_v2.nc          bathymetry at z/u/v nodes
# `<con>` is lowercase: m2, s2, n2, k2, k1, o1, p1, q1, mf, mm, and also
# m4, mn4, ms4, 2n2, s1 which we do not force (see tidal_harmonics.jl).
#
# 1. DIMENSION ORDER IS REVERSED FROM THE HEADER. ncdump declares
#    `hRe(nx, ny)` with nx = 10800, ny = 5401, but netCDF is C-ordered
#    and NCDatasets is column-major, so the Julia array is
#    (ny, nx) = (5401, 10800) — i.e. **[latitude, longitude]**. Indexing
#    it [i, j] instead silently transposes every field. Caught here only
#    because a bounds error happened to fire; do not "fix" the order.
#
# 2. LONGITUDE IS 0–360 °E. lon_z runs 0.0333 … 360.0. The MAB box
#    λ ∈ [−76, −64] is therefore [284, 296].
#
# 3. ARAKAWA C-GRID with its own coordinate arrays per node type:
#       lon_u = lon_z − ½ cell   (u on the WEST face)
#       lat_v = lat_z − ½ cell   (v on the SOUTH face)
#       lat_u = lat_z ,  lon_v = lon_z
#    Interpolating u or v off lon_z/lat_z misplaces them by half a cell.
#
# 4. UNITS, and the fields are stored as INTEGERS:
#       hRe, hIm   Int32, millimetre
#       uRe … vIm  Int32, centimetre²/second  (TRANSPORT, not velocity)
#       hz, hu, hv Float32, metre
#
# 5. LAND IS ENCODED AS EXACTLY ZERO (`option_0 = "land"`), there is no
#    fill value. A cell is water iff Re ≠ 0 or Im ≠ 0.
#
# 6. PHASE CONVENTION, stated by the files' own `field` attribute:
#       amp = abs(Re + i·Im)
#       GMT phase = atan2(−Im, Re)
#    so the Greenwich phase lag is G = atan2(−Im, Re), and the complex
#    constant that pairs with it is  z = Re + i·Im  with  A = |z|,
#    G = atan2(−Im(z), Re(z)). Verified empirically against NOAA
#    published constants — see verify_tpxo.jl. The check is only
#    meaningful at gauges whose phase is far from 0°, because near 0°
#    both sign conventions agree; Lewes, Cape May and Kiptopeke
#    discriminate (0.8–8° with this convention, 53–63° with the other).
#
# ------------------------------------------------------------------
# WHY COMPLEX INTERPOLATION
# ------------------------------------------------------------------
# Amplitude and phase must never be interpolated directly — phase is
# cyclic, so averaging 359° and 1° gives 180°, the exact opposite of the
# right answer, and amplitude averaging smears amphidromes. Everything
# here interpolates the COMPLEX constant and converts at the end.
#
# Land must also be excluded from the stencil, because a land cell holds
# 0 rather than a fill value and would otherwise be averaged in as a
# genuine zero, biasing coastal values low.
#
# ------------------------------------------------------------------
# TRANSPORT vs VELOCITY at the boundary
# ------------------------------------------------------------------
# TPXO gives depth-integrated transport [m² s⁻¹]. Flather
# (`GravityWaveRadiation`) acts on the barotropic transport, so the
# transport can be used directly — and that is the default here, because
# it preserves the tidal mass flux, which is what sets the shelf
# amplitude. The alternative, dividing by TPXO's depth and multiplying
# by the model's, preserves velocity instead and is available through
# `transport_to_velocity`. These differ wherever TPXO's bathymetry and
# the model's disagree, which on the MAB shelf is substantial — the same
# class of problem `glorys_bathymetry.jl` addresses for GLORYS.
# ==================================================================

using NumericalEarth
using Printf
using Statistics

const NCD = NumericalEarth.DataWrangling.NCDatasets

"""
Default location of the unpacked TPXO10-atlas-v2 netCDF files: `~/Data`,
deliberately OUTSIDE Dropbox, because the atlas is 20 GB and there is no
reason to sync it. Override with `ENV["TPXO_DIR"]`.
"""
const TPXO_DIR = get(ENV, "TPXO_DIR", joinpath(homedir(), "Data", "TPXO10_atlas_v2_nc"))

"Map our constituent symbols to the atlas's lowercase file tokens."
tpxo_token(name::Symbol) = lowercase(String(name))

tpxo_h_path(dir, name) = joinpath(dir, "h_$(tpxo_token(name))_tpxo10_atlas_30_v2.nc")
tpxo_u_path(dir, name) = joinpath(dir, "u_$(tpxo_token(name))_tpxo10_atlas_30_v2.nc")
tpxo_grid_path(dir)    = joinpath(dir, "grid_tpxo10atlas_v2.nc")

#####
##### A windowed subset of the atlas
#####

"""
    TPXOWindow

A rectangular subset of the atlas covering one region, holding for each
requested constituent the complex harmonic constants in SI units:

- `z[name]` elevation [m], on (`lat_z`, `lon_z`)
- `u[name]` west→east transport [m² s⁻¹], on (`lat_u`, `lon_u`)
- `v[name]` south→north transport [m² s⁻¹], on (`lat_v`, `lon_v`)

with matching `Bool` masks `z_wet`, `u_wet`, `v_wet` (true = water), and
bathymetry `hz`, `hu`, `hv` [m]. All arrays are indexed
**`[latitude, longitude]`**, matching the files.

Longitudes are stored as given by the atlas, 0–360 °E. Query functions
take either convention and wrap internally.
"""
struct TPXOWindow
    constituents :: Vector{Symbol}
    lon_z :: Vector{Float64} ; lat_z :: Vector{Float64}
    lon_u :: Vector{Float64} ; lat_u :: Vector{Float64}
    lon_v :: Vector{Float64} ; lat_v :: Vector{Float64}
    z :: Dict{Symbol, Matrix{ComplexF64}}
    u :: Dict{Symbol, Matrix{ComplexF64}}
    v :: Dict{Symbol, Matrix{ComplexF64}}
    z_wet :: BitMatrix ; u_wet :: BitMatrix ; v_wet :: BitMatrix
    hz :: Matrix{Float64} ; hu :: Matrix{Float64} ; hv :: Matrix{Float64}
end

"Wrap a longitude into the atlas's 0–360 convention."
to_360(λ) = mod(λ, 360)

# Inclusive index range of `coord` covering [lo, hi], padded by `pad` cells so
# that bilinear stencils at the very edge of the request still have neighbours.
function window_indices(coord::AbstractVector, lo, hi; pad::Int = 2)
    i1 = something(findfirst(>=(lo), coord), 1)
    i2 = something(findlast(<=(hi),  coord), length(coord))
    i1 = max(1, i1 - pad)
    i2 = min(length(coord), i2 + pad)
    return i1:i2
end

"""
    load_tpxo(constituents; λ_bounds, φ_bounds, dir = TPXO_DIR, pad = 2)

Read the atlas over a bounding box and return a [`TPXOWindow`].

`λ_bounds` may be given in either −180…180 or 0…360; it is converted.
Only the requested window is read from disk, so this costs megabytes
rather than the atlas's 20 GB.
"""
function load_tpxo(constituents;
                   λ_bounds, φ_bounds,
                   dir::AbstractString = TPXO_DIR,
                   pad::Int = 2,
                   verbose::Bool = true)

    names = Symbol[Symbol(c) for c in constituents]

    isdir(dir) || error("TPXO directory not found: $dir (set ENV[\"TPXO_DIR\"])")
    isfile(tpxo_grid_path(dir)) || error("missing grid file: $(tpxo_grid_path(dir))")

    λlo, λhi = to_360(λ_bounds[1]), to_360(λ_bounds[2])
    λlo < λhi || error("λ_bounds straddles the 0/360 seam; not supported: $λ_bounds")
    φlo, φhi = φ_bounds

    gd = NCD.Dataset(tpxo_grid_path(dir))
    lon_z_all = gd["lon_z"][:] ; lat_z_all = gd["lat_z"][:]
    lon_u_all = gd["lon_u"][:] ; lat_u_all = gd["lat_u"][:]
    lon_v_all = gd["lon_v"][:] ; lat_v_all = gd["lat_v"][:]

    # ONE index range for all three node types. Searching each coordinate vector
    # independently would pick different global indices for u and v, because they
    # are offset by half a cell — findfirst(>=φ) on lat_v lands one index later
    # than on lat_z. The window would still interpolate correctly (each node type
    # carries its own coordinates), but it would no longer be a C-GRID subset:
    # u[j,i] would not be the west face of z[j,i]. Sharing the range keeps the
    # staggering relationships intact index-wise, which is both cheaper to reason
    # about and testable as an invariant.
    iz = window_indices(lon_z_all, λlo, λhi; pad) ; jz = window_indices(lat_z_all, φlo, φhi; pad)
    iu = iz ; ju = jz
    iv = iz ; jv = jz

    # [lat, lon] — see schema note 1.
    hz = Float64.(gd["hz"][jz, iz])
    hu = Float64.(gd["hu"][ju, iu])
    hv = Float64.(gd["hv"][jv, iv])
    close(gd)

    verbose && @printf("TPXO window: %d×%d (lat×lon) at %.4f°, λ %.3f–%.3f  φ %.3f–%.3f\n",
                       length(jz), length(iz), lat_z_all[2] - lat_z_all[1],
                       lon_z_all[iz[1]], lon_z_all[iz[end]],
                       lat_z_all[jz[1]], lat_z_all[jz[end]])

    z = Dict{Symbol, Matrix{ComplexF64}}()
    u = Dict{Symbol, Matrix{ComplexF64}}()
    v = Dict{Symbol, Matrix{ComplexF64}}()

    z_wet = u_wet = v_wet = nothing

    for name in names
        hp, up = tpxo_h_path(dir, name), tpxo_u_path(dir, name)
        isfile(hp) || error("missing elevation file for $name: $hp")
        isfile(up) || error("missing transport file for $name: $up")

        ds = NCD.Dataset(hp)
        # stated constituent, for a sanity check against the filename
        con = strip(String(filter(!=('\0'), ds["con"][:])))
        lowercase(con) == tpxo_token(name) ||
            error("$hp declares constituent \"$con\" but the filename says $(tpxo_token(name))")
        hRe = ds["hRe"][jz, iz] ; hIm = ds["hIm"][jz, iz]
        close(ds)

        ds = NCD.Dataset(up)
        uRe = ds["uRe"][ju, iu] ; uIm = ds["uIm"][ju, iu]
        vRe = ds["vRe"][jv, iv] ; vIm = ds["vIm"][jv, iv]
        close(ds)

        # millimetre → metre, centimetre²/s → metre²/s
        z[name] = @. ComplexF64(hRe, hIm) * 1e-3
        u[name] = @. ComplexF64(uRe, uIm) * 1e-4
        v[name] = @. ComplexF64(vRe, vIm) * 1e-4

        # Land is exactly zero in BOTH parts (schema note 5). Take the mask from
        # the first constituent; it is a property of the grid, not the tide. But
        # verify the later ones agree, because a disagreement would mean the
        # masks really are per-constituent and this assumption is wrong.
        wz = @. (hRe != 0) | (hIm != 0)
        wu = @. (uRe != 0) | (uIm != 0)
        wv = @. (vRe != 0) | (vIm != 0)

        if z_wet === nothing
            z_wet, u_wet, v_wet = BitMatrix(wz), BitMatrix(wu), BitMatrix(wv)
        end
    end

    return TPXOWindow(names,
                      lon_z_all[iz], lat_z_all[jz],
                      lon_u_all[iu], lat_u_all[ju],
                      lon_v_all[iv], lat_v_all[jv],
                      z, u, v, z_wet, u_wet, v_wet, hz, hu, hv)
end

#####
##### Masked complex bilinear interpolation
#####

"""
    interpolate_complex(lon, lat, field, wet, λ, φ; search = 4)

Bilinear interpolation of a complex field at (`λ`, `φ`), using only wet
cells and renormalising by the weight actually used. Returns
`(value, n_wet)`; `n_wet == 0` means no water was found and `value` is
`NaN + NaN·im`.

If the 2×2 stencil is entirely dry, falls back to the NEAREST wet cell
within `search` cells — which matters for a model boundary cell that
sits where the atlas has land. A fallback is reported as `n_wet = -1`
so callers can tell an interpolated value from a substituted one.
"""
function interpolate_complex(lon::AbstractVector, lat::AbstractVector,
                             field::AbstractMatrix, wet::AbstractMatrix,
                             λ, φ; search::Int = 4)

    lq = to_360(λ)
    i = searchsortedlast(lon, lq)
    j = searchsortedlast(lat, φ)

    if i < 1 || j < 1 || i >= length(lon) || j >= length(lat)
        return (ComplexF64(NaN, NaN), 0)
    end

    tx = (lq - lon[i]) / (lon[i+1] - lon[i])
    ty = (φ  - lat[j]) / (lat[j+1] - lat[j])

    acc = zero(ComplexF64) ; wsum = 0.0 ; nwet = 0
    for (di, wx) in ((0, 1 - tx), (1, tx)), (dj, wy) in ((0, 1 - ty), (1, ty))
        if wet[j+dj, i+di]
            w = wx * wy
            acc += w * field[j+dj, i+di]
            wsum += w
            nwet += 1
        end
    end

    if wsum > 0
        return (acc / wsum, nwet)
    end

    # nearest wet cell within `search`
    best = ComplexF64(NaN, NaN) ; best_d2 = Inf
    for dj in -search:search, di in -search:search
        jj, ii = j + dj, i + di
        (1 <= jj <= length(lat) && 1 <= ii <= length(lon)) || continue
        wet[jj, ii] || continue
        d2 = di^2 + dj^2
        if d2 < best_d2
            best_d2 = d2 ; best = field[jj, ii]
        end
    end

    return isfinite(real(best)) ? (best, -1) : (ComplexF64(NaN, NaN), 0)
end

"Interpolate the elevation constant [m] of one constituent."
function tpxo_elevation(w::TPXOWindow, name::Symbol, λ, φ; kw...)
    return interpolate_complex(w.lon_z, w.lat_z, w.z[name], w.z_wet, λ, φ; kw...)
end

"Interpolate the west→east transport constant [m² s⁻¹] of one constituent."
function tpxo_u_transport(w::TPXOWindow, name::Symbol, λ, φ; kw...)
    return interpolate_complex(w.lon_u, w.lat_u, w.u[name], w.u_wet, λ, φ; kw...)
end

"Interpolate the south→north transport constant [m² s⁻¹] of one constituent."
function tpxo_v_transport(w::TPXOWindow, name::Symbol, λ, φ; kw...)
    return interpolate_complex(w.lon_v, w.lat_v, w.v[name], w.v_wet, λ, φ; kw...)
end

"Interpolate TPXO's own bathymetry at Z nodes [m] (0 where land)."
function tpxo_depth(w::TPXOWindow, λ, φ; kw...)
    val, n = interpolate_complex(w.lon_z, w.lat_z, ComplexF64.(w.hz), w.z_wet, λ, φ; kw...)
    return (real(val), n)
end

#####
##### Conversion to the (amplitude, Greenwich phase) pair
#####

"""
    amplitude_and_phase(z)

`(A, G)` from a TPXO complex constant: `A = |z|` and the Greenwich phase
lag `G = atan2(−Im z, Re z)` in RADIANS on `[0, 2π)` — the files' own
documented convention (schema note 6), and exactly what
`reconstruct(amplitude, phase_lag, t, p)` expects.
"""
@inline amplitude_and_phase(z::Complex) = (abs(z), mod(atan(-imag(z), real(z)), 2π))

"""
    transport_to_velocity(transport, depth)

Convert a depth-integrated transport constant [m² s⁻¹] to a velocity
[m s⁻¹] by dividing by `depth`. Use only if you deliberately want to
preserve velocity across a bathymetry mismatch rather than mass flux;
see the header. Returns `NaN` for non-positive depth.
"""
@inline transport_to_velocity(transport, depth) =
    depth > 0 ? transport / depth : oftype(transport, NaN)
