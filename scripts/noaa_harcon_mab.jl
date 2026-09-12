# ==================================================================
# NOAA CO-OPS published harmonic constants for Mid-Atlantic Bight tide
# gauges -- the ground truth the model tide is validated against.
#
# These are OBSERVED constants from multi-year records (the `comments`
# field of the API records the analysis period; most are a vector
# analysis of 5 years). They serve two purposes:
#
#   1. VALIDATION. Model sea level, harmonically analysed at a gauge
#      location, should reproduce these amplitudes and Greenwich phase
#      lags.
#   2. PINNING THE TPXO PHASE CONVENTION. TPXO stores complex Re/Im
#      parts, and whether the Greenwich lag is atan2(-Im, Re) or
#      atan2(Im, Re) decides which way the tide propagates. Comparing
#      TPXO interpolated to a gauge against that gauge's published
#      constants settles the sign empirically instead of by assumption.
#
# `amplitude` is in METRES and `phase` is the GREENWICH phase lag in
# DEGREES -- i.e. exactly the (A, G) convention of
# `reconstruct(amplitude, phase_lag, t, p)` in tidal_harmonics.jl, once
# the phase is converted to radians.
#
# CAVEAT on using these: a coastal gauge at 1/12 deg often lands in a
# land or immersed cell, and gauges inside estuaries feel resonance the
# model cannot resolve. `kind` records this. Only :open stations are
# fair comparisons for a shelf model; Atlantic City and Duck are the
# two best MAB targets. Mf and Mm are reported as 0.0 at most of these
# stations -- NOAA does not resolve the long-period constituents here,
# so they cannot be validated against gauges.
#
# Regenerated with:
#   curl -s "https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/\
#            stations/<ID>/harcon.json?units=metric"
#   plus stations/<ID>.json for lat/lon.
# Fetched 2026-09-10.
# ==================================================================

"""
NOAA CO-OPS observed harmonic constants, keyed by station ID.
Each value has `name`, `lat`, `lon`, `kind`, `note`, and `con`, a
`Dict{Symbol, Tuple{Float64, Float64}}` of `constituent => (amplitude_m,
greenwich_phase_deg)`.
"""
const NOAA_HARCON = Dict(
    "8447930" => (name = "Woods Hole, MA", lat = 41.5236, lon = -70.6711, kind = :complex,
        note = "Vineyard Sound -- strongly unresolved, outside the shelf regime",
        con = Dict(:M2 => (0.2290, 35.3), :S2 => (0.0550, 36.9), :N2 => (0.0760, 20.9), :K2 => (0.0140, 31.1), :K1 => (0.0680, 191.2), :O1 => (0.0620, 201.5), :P1 => (0.0260, 201.1), :Q1 => (0.0140, 191.4), :Mf => (0.0000, 0.0), :Mm => (0.0000, 0.0))),
    "8461490" => (name = "New London, CT", lat = 41.3717, lon = -72.0956, kind = :complex,
        note = "Long Island Sound -- unresolved",
        con = Dict(:M2 => (0.3590, 58.7), :S2 => (0.0640, 70.4), :N2 => (0.0850, 34.0), :K2 => (0.0160, 72.1), :K1 => (0.0720, 180.1), :O1 => (0.0500, 209.1), :P1 => (0.0240, 193.7), :Q1 => (0.0130, 196.7), :Mf => (0.0000, 0.0), :Mm => (0.0000, 0.0))),
    "8510560" => (name = "Montauk, NY", lat = 41.0483, lon = -71.9594, kind = :complex,
        note = "mouth of Long Island Sound -- reclassified from :open. Geographically it \
                looks like open coast, but it sits at the entrance to the LIS/Block Island/\
                Nantucket Sound system, which is strongly resonant and unresolved at 1/12 deg. \
                Evidence: its TPXO stencil is half land, TPXO itself misses the gauge phase by \
                42 deg, and the M2 amphidrome sits ~40 km away, so the model cannot place the \
                node correctly either. Not a fair target for a shelf model.",
        con = Dict(:M2 => (0.2830, 48.2), :S2 => (0.0590, 58.6), :N2 => (0.0760, 25.6), :K2 => (0.0150, 59.8), :K1 => (0.0700, 178.9), :O1 => (0.0480, 207.6), :P1 => (0.0240, 191.6), :Q1 => (0.0130, 201.6), :Mf => (0.0000, 0.0), :Mm => (0.0000, 0.0))),
    "8531680" => (name = "Sandy Hook, NJ", lat = 40.4669, lon = -74.0094, kind = :estuary,
        note = "inside Sandy Hook / Raritan Bay -- NY Bight apex amplification",
        con = Dict(:M2 => (0.6790, 5.6), :S2 => (0.1300, 32.4), :N2 => (0.1560, 350.0), :K2 => (0.0360, 32.5), :K1 => (0.1030, 173.9), :O1 => (0.0510, 172.1), :P1 => (0.0330, 174.1), :Q1 => (0.0120, 189.6), :Mf => (0.0000, 0.0), :Mm => (0.0000, 0.0))),
    "8534720" => (name = "Atlantic City, NJ", lat = 39.3567, lon = -74.4180, kind = :open,
        note = "open Atlantic coast -- best MAB target",
        con = Dict(:M2 => (0.5810, 355.5), :S2 => (0.1150, 20.1), :N2 => (0.1370, 337.3), :K2 => (0.0310, 18.3), :K1 => (0.1090, 181.3), :O1 => (0.0740, 167.4), :P1 => (0.0340, 179.1), :Q1 => (0.0140, 165.8), :Mf => (0.0000, 0.0), :Mm => (0.0000, 0.0))),
    "8536110" => (name = "Cape May, NJ", lat = 38.9683, lon = -74.9600, kind = :estuary,
        note = "Delaware Bay mouth",
        con = Dict(:M2 => (0.6970, 28.6), :S2 => (0.1220, 56.0), :N2 => (0.1510, 8.8), :K2 => (0.0340, 55.1), :K1 => (0.1040, 200.2), :O1 => (0.0820, 186.8), :P1 => (0.0350, 197.3), :Q1 => (0.0150, 179.4), :Mf => (0.0000, 0.0), :Mm => (0.0000, 0.0))),
    "8557380" => (name = "Lewes, DE", lat = 38.7828, lon = -75.1193, kind = :estuary,
        note = "inside Delaware Bay",
        con = Dict(:M2 => (0.5880, 31.2), :S2 => (0.1050, 57.2), :N2 => (0.1320, 9.3), :K2 => (0.0280, 57.9), :K1 => (0.1020, 202.2), :O1 => (0.0810, 189.4), :P1 => (0.0330, 198.1), :Q1 => (0.0150, 178.9), :Mf => (0.0000, 0.0), :Mm => (0.0000, 0.0))),
    "8570283" => (name = "Ocean City Inlet, MD", lat = 38.3283, lon = -75.0917, kind = :inlet,
        note = "behind a barrier inlet -- unresolved at 1/12 deg",
        con = Dict(:M2 => (0.3040, 8.6), :S2 => (0.0550, 31.5), :N2 => (0.0700, 348.5), :K2 => (0.0150, 30.5), :K1 => (0.0550, 210.5), :O1 => (0.0520, 200.0), :P1 => (0.0190, 198.9), :Q1 => (0.0080, 191.3), :Mf => (0.0000, 0.0), :Mm => (0.0270, 69.5))),
    "8632200" => (name = "Kiptopeke, VA", lat = 37.1652, lon = -75.9884, kind = :estuary,
        note = "lower Chesapeake Bay",
        con = Dict(:M2 => (0.3830, 32.7), :S2 => (0.0660, 57.5), :N2 => (0.0860, 14.3), :K2 => (0.0190, 54.6), :K1 => (0.0570, 194.3), :O1 => (0.0470, 215.5), :P1 => (0.0180, 196.6), :Q1 => (0.0120, 195.7), :Mf => (0.0000, 0.0), :Mm => (0.0000, 0.0))),
    "8638610" => (name = "Sewells Point, VA", lat = 36.9428, lon = -76.3286, kind = :estuary,
        note = "Hampton Roads; also west of the model box",
        con = Dict(:M2 => (0.3500, 47.1), :S2 => (0.0640, 76.2), :N2 => (0.0800, 27.5), :K2 => (0.0170, 71.8), :K1 => (0.0510, 199.8), :O1 => (0.0400, 223.8), :P1 => (0.0150, 198.0), :Q1 => (0.0110, 209.5), :Mf => (0.0000, 0.0), :Mm => (0.0000, 0.0))),
    "8651370" => (name = "Duck, NC", lat = 36.1833, lon = -75.7467, kind = :open,
        note = "open Atlantic coast, Outer Banks -- good MAB target",
        con = Dict(:M2 => (0.4740, 358.0), :S2 => (0.0860, 22.1), :N2 => (0.1150, 338.5), :K2 => (0.0220, 21.7), :K1 => (0.0850, 172.9), :O1 => (0.0560, 190.9), :P1 => (0.0280, 171.8), :Q1 => (0.0150, 187.2), :Mf => (0.0000, 0.0), :Mm => (0.0000, 0.0))),
)

"Station IDs that are fair comparisons for a shelf-resolving model."
open_coast_stations() = [k for (k, v) in NOAA_HARCON if v.kind === :open]

"""
    noaa_constants(station_id, constituents)

`(amplitude, phase_lag)` vectors for `constituents` at `station_id`,
amplitudes in metres and phase lags in RADIANS -- ready to hand to
[`reconstruct`].
"""
function noaa_constants(station_id::AbstractString, constituents)
    con = NOAA_HARCON[station_id].con
    # Tolerate case differences: NOAA returns MF/MM, our tables use Mf/Mm.
    lookup = Dict(lowercase(String(k)) => v for (k, v) in con)
    get_con(c) = get(lookup, lowercase(String(c))) do
        throw(KeyError("constituent $c not available at station $station_id"))
    end
    amplitude = [get_con(c)[1] for c in constituents]
    phase_lag = [deg2rad(get_con(c)[2]) for c in constituents]
    return amplitude, phase_lag
end

"NOAA published angular speeds [deg/hour], an INDEPENDENT check on the frequency table."
const NOAA_SPEEDS_DEG_PER_HOUR = Dict(
    :M2 => 28.984104, :S2 => 30.000000, :N2 => 28.439730, :K2 => 30.082138,
    :K1 => 15.041069, :O1 => 13.943035, :P1 => 14.958931, :Q1 => 13.398661,
    :Mf =>  1.098033, :Mm =>  0.544375)
