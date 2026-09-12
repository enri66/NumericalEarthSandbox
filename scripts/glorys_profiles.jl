# ==================================================================
# Horizontally uniform stratification profiles from GLORYS.
#
# Used by the stratified tide experiments. The profile must be HORIZONTALLY
# UNIFORM: any horizontal density gradient spins up a baroclinic circulation from
# rest, and that flow would contaminate the tidal harmonic analysis the
# experiments rely on. So the question is only WHICH profile — and in April the
# MAB has no representative one. Measured from GLORYS on 2019-04-01:
#
#                 surface T   surface S   N² at 10 m   set by
#   shelf <200 m     6.9 °C      31.1      6.0e-4      salinity (T nearly uniform)
#   deep >1000 m    18.6 °C      36.1      1.7e-5      weak upper ocean; thermocline
#   whole box       15.8 °C      34.9      1.5e-4      an averaging artefact
#
# Shelf and deep water differ by 12 °C and 5 psu at the surface — the shelf-break
# front — and the whole-box mean describes neither. Hence two separate profiles,
# testing two different mechanisms:
#   :deep  — realistic stratification at the shelf break and slope, where
#            barotropic-to-internal tide conversion happens (the energy sink)
#   :shelf — realistic stratification on the shelf, where the model's error is
#            (direct effect on shelf wave propagation)
#
# CAVEATS
#   - GLORYS supplies POTENTIAL temperature (`thetao`) and PRACTICAL salinity
#     (`so`); TEOS-10 formally wants conservative temperature and absolute
#     salinity. The difference is a few hundredths of a degree and ~0.16 g/kg,
#     nearly uniform, and immaterial for a stratification experiment.
#   - Beyond the data's depth range the profile is held CONSTANT. For the shelf
#     class that range ends near 190 m, so a shelf-profile run is "shelf
#     stratification above ~190 m, unstratified below".
# ==================================================================

using NumericalEarth
using Dates, Statistics

const NCD_profiles = NumericalEarth.DataWrangling.NCDatasets

"""
    glorys_mean_profile(dir, date; region = :deep,
                        box = "-76.0_-64.0_34.0_42.0")

Area-weighted (cos φ) mean GLORYS potential temperature and practical salinity
over the columns of one depth class, on the GLORYS depth levels (positive down).
`region` is `:deep` (column depth > 1000 m), `:shelf` (< 200 m) or `:all`.

Column depth is taken as the deepest GLORYS level holding data. Returns
`(; depth, T, S, ncolumns)`, truncated where fewer than 20 columns contribute.
"""
function glorys_mean_profile(dir, date::DateTime; region::Symbol = :deep,
                             box = "-76.0_-64.0_34.0_42.0")

    stamp = Dates.format(date, "yyyy-mm-ddTHH-MM-SS")
    tag   = "GLORYSDaily_$(stamp)_$(stamp)_$(box).nc"
    pT, pS = joinpath(dir, "thetao_" * tag), joinpath(dir, "so_" * tag)
    isfile(pT) || error("missing GLORYS temperature file: $pT")
    isfile(pS) || error("missing GLORYS salinity file: $pS")

    dsT = NCD_profiles.Dataset(pT) ; dsS = NCD_profiles.Dataset(pS)
    T   = dsT["thetao"][:, :, :, 1] ; S = dsS["so"][:, :, :, 1]
    dep = Float64.(dsT["depth"][:]) ; lat = Float64.(dsT["latitude"][:])
    close(dsT) ; close(dsS)

    # netCDF is C-ordered and NCDatasets column-major: the Julia array is
    # (lon, lat, depth). Checked, not assumed — the TPXO reader was bitten by this.
    size(T, 3) == length(dep) || error("unexpected GLORYS dimension order: $(size(T))")

    valid = .!ismissing.(T)
    H = [(k = findlast(valid[i, j, :]); k === nothing ? 0.0 : dep[k])
         for i in axes(T, 1), j in axes(T, 2)]
    w = [cosd(lat[j]) for i in axes(T, 1), j in axes(T, 2)]

    in_class = region === :deep  ? (h -> h > 1000) :
               region === :shelf ? (h -> 0 < h < 200) :
               region === :all   ? (h -> h > 0) :
               throw(ArgumentError("region must be :deep, :shelf or :all, got $region"))

    depth = Float64[] ; Tp = Float64[] ; Sp = Float64[]
    ncolumns = count(in_class, H)
    for k in eachindex(dep)
        m = [valid[i, j, k] && in_class(H[i, j]) for i in axes(T, 1), j in axes(T, 2)]
        count(m) < 20 && break
        push!(depth, dep[k])
        push!(Tp, sum(w[m] .* Float64.(T[:, :, k][m])) / sum(w[m]))
        push!(Sp, sum(w[m] .* Float64.(S[:, :, k][m])) / sum(w[m]))
    end

    return (; depth, T = Tp, S = Sp, ncolumns)
end

"""
    profile_function(depth, values)

A function of Oceananigans' `z` (negative downward) that interpolates `values`
linearly in depth and holds the end values beyond the data range.
"""
function profile_function(depth, values)
    d = collect(depth) ; v = collect(values)
    return function (z)
        dz = -z
        dz <= d[1]   && return v[1]
        dz >= d[end] && return v[end]
        k = searchsortedlast(d, dz)
        θ = (dz - d[k]) / (d[k+1] - d[k])
        return (1 - θ) * v[k] + θ * v[k+1]
    end
end
