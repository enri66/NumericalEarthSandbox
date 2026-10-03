# Near-inertial currents and vertical shear of a script-05 run (MAB_MOORINGS=pioneer) against the OOI Coastal Pioneer
# New England Shelf current profilers. Does the model's storm response put too much near-inertial shear across the base
# of the mixed layer? For each mooring: the ADCP's hourly profiles (ERDDAP CSV, see OOI_DIR) and the model's hourly
# column at the nearest wet cell, both interpolated to a common 4 m (shelf) or 8 m (slope) grid between ZTOP (below the
# ADCP's surface side-lobe contamination) and a few bins above the bottom, then
#   near-inertial velocity: an FFT band-pass of each depth's hourly series to periods 16-22 h (f at 40°N: 18.6 h;
#       the M2/S2 and K1/O1 tides lie outside the band), gaps filled linearly first;
#   near-inertial kinetic energy, depth-mean over the upper layer, as time series;
#   squared vertical shear of the hourly velocity and of its near-inertial part, over each bin pair, as time-mean
#       profiles for the whole record and for the storm windows.
# Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/catke_fall/moor OOI_DIR=/t0/workdir/enrique/Data/OOI/pioneer \
#       julia --project=. scripts/moorings_vs_ooi.jl
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))
const FFTW = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "FFTW")]

const TAG     = get(ENV, "MAB_TAG", "moor")
const PREFIX  = run_prefix(TAG)
const OOI_DIR = get(ENV, "OOI_DIR", joinpath(homedir(), "Data", "OOI", "pioneer"))
const BAND    = (16.0, 22.0)                                   # near-inertial periods (hours)
const STORMS  = (("Dorian, 6-12 Sep", DateTime(2019, 9, 6), DateTime(2019, 9, 13)),
                 ("17-23 Oct", DateTime(2019, 10, 17), DateTime(2019, 10, 24)))

# name in the run => (ERDDAP dataset, vertical grid spacing, top of the usable ADCP range)
const MOORINGS = Dict("OSSM" => ("ooi-cp04ossm-mfd35-01-adcpsj000", 8.0, 40.0),
                      "PMUO" => ("ooi-cp02pmuo-rii01-02-adcpsl010", 8.0, 40.0),
                      "PMCO" => ("ooi-cp02pmco-rii01-02-adcptg010", 4.0, 12.0),
                      "CNSM" => ("ooi-cp01cnsm-mfd35-01-adcptf000", 4.0, 12.0))

# ---------------- ADCP ----------------
# Hourly (u, v) on the ADCP's bins; with QC columns only flags 1 (pass) and 2 (not evaluated) are kept
function read_adcp(dataset)
    lines = readlines(joinpath(OOI_DIR, dataset * ".csv"))
    header = split(lines[1], ","); col = Dict(h => i for (i, h) in enumerate(header))
    qc = haskey(col, "eastward_sea_water_velocity_qc_agg")
    data = Dict{Tuple{DateTime, Float64}, NTuple{2, Float64}}()
    for line in lines[3:end]
        v = split(line, ",")
        qc && !(v[col["eastward_sea_water_velocity_qc_agg"]] in ("1", "2")) && continue
        t = floor(DateTime(v[col["time"]][1:19]), Hour)
        u = tryparse(Float64, v[col["eastward_sea_water_velocity"]]); w = tryparse(Float64, v[col["northward_sea_water_velocity"]])
        (isnothing(u) || isnothing(w) || abs(u) > 3 || abs(w) > 3) && continue
        d = -parse(Float64, v[col["z"]])
        d > 0 && (data[(t, d)] = (u, w))                       # bins above the sea surface are not data
    end
    times = sort(unique(first.(keys(data)))); depths = sort(unique(last.(keys(data))))
    U = fill(NaN, length(depths), length(times)); V = fill(NaN, length(depths), length(times))
    ti = Dict(t => n for (n, t) in enumerate(times)); di = Dict(d => n for (n, d) in enumerate(depths))
    for ((t, d), (u, w)) in data
        U[di[d], ti[t]] = u; V[di[d], ti[t]] = w
    end
    return times, depths, U, V
end

# ---------------- model ----------------
function read_model(name)
    for f in filter(f -> occursin(r"_moorings_rank\d+\.jld2$", f) && startswith(basename(f), basename(PREFIX) * "_moorings"),
                    readdir(dirname(PREFIX); join = true))
        out = JLD2.jldopen(f) do file
            haskey(file, name) || return nothing
            t0 = DateTime(file["start_date"])
            times = [round(t0 + Second(round(Int, s)), Hour) for s in file["$name/time"]]   # the first record is one time step after a pickup
            (; times, depths = reverse(-file["z_center"]), U = reverse(file["$name/u"]; dims = 1),
               V = reverse(file["$name/v"]; dims = 1), cell = file["$name/cell"])
        end
        isnothing(out) || return out
    end
    error("no model column for $name")
end

# ---------------- processing ----------------
# Profiles onto the common grid (NaN outside each profile's valid range), then onto a common hourly time axis
function regrid(depths, A, grid)
    out = fill(NaN, length(grid), size(A, 2))
    for n in axes(A, 2)
        ok = isfinite.(A[:, n]); d = depths[ok]; a = A[ok, n]
        length(d) < 2 && continue
        for (m, z) in enumerate(grid)
            d[1] <= z <= d[end] || continue
            k = clamp(searchsortedlast(d, z), 1, length(d) - 1)
            out[m, n] = a[k] + (a[k+1] - a[k]) * (z - d[k]) / (d[k+1] - d[k])
        end
    end
    return out
end
on_axis(times, A, axis) = (ti = Dict(t => n for (n, t) in enumerate(times));
                           [haskey(ti, t) ? A[m, ti[t]] : NaN for m in axes(A, 1), t in axis])

# Linear gap filling (series with less than half the record valid stay NaN), then an FFT band-pass
function fill_gaps(x)
    ok = findall(isfinite, x)
    length(ok) < length(x) ÷ 2 && return fill(NaN, length(x))
    y = copy(x)
    for n in eachindex(y)
        isfinite(y[n]) && continue
        a = findlast(<(n), ok); b = findfirst(>(n), ok)
        y[n] = isnothing(a) ? x[ok[b]] : isnothing(b) ? x[ok[a]] : x[ok[a]] + (x[ok[b]] - x[ok[a]]) * (n - ok[a]) / (ok[b] - ok[a])
    end
    return y
end
function bandpass(x; dt = 1.0)
    y = fill_gaps(x); any(isnan, y) && return y
    Y = FFTW.fft(y .- sum(y) / length(y))
    freq = FFTW.fftfreq(length(y), 1 / dt)                     # cycles per hour
    Y[.!(1 / BAND[2] .<= abs.(freq) .<= 1 / BAND[1])] .= 0
    return real(FFTW.ifft(Y))
end
bandpass_rows(A) = reduce(vcat, [bandpass(A[m, :])' for m in axes(A, 1)])

shear²(U, V, Δz) = ((U[2:end, :] .- U[1:end-1, :]) .^ 2 .+ (V[2:end, :] .- V[1:end-1, :]) .^ 2) ./ Δz^2
nanmean(x) = (f = filter(isfinite, x); isempty(f) ? NaN : sum(f) / length(f))
rowmean(A, cols = axes(A, 2)) = [nanmean(A[m, cols]) for m in axes(A, 1)]

# ---------------- comparison ----------------
results = Dict{String, Any}()
for name in sort(collect(keys(MOORINGS)))
    dataset, Δz, ztop = MOORINGS[name]
    at, ad, AU, AV = read_adcp(dataset)
    m = read_model(name)
    zbottom = min(maximum(ad), m.cell[5]) - 3Δz
    grid = collect(ztop:Δz:zbottom)
    axis = collect(max(first(at), first(m.times)):Hour(1):min(last(at), last(m.times)))
    keep = [true; diff(m.times) .> Hour(0)]                    # drop a duplicate hour from the rounding
    m = merge(m, (times = m.times[keep], U = m.U[:, keep], V = m.V[:, keep]))
    obs = (U = on_axis(at, regrid(ad, AU, grid), axis), V = on_axis(at, regrid(ad, AV, grid), axis))
    mdl = (U = on_axis(m.times, regrid(m.depths, m.U, grid), axis), V = on_axis(m.times, regrid(m.depths, m.V, grid), axis))
    ni(x) = (U = bandpass_rows(x.U), V = bandpass_rows(x.V))
    obs_ni, mdl_ni = ni(obs), ni(mdl)
    results[name] = (; grid, axis, Δz, obs, mdl, obs_ni, mdl_ni, cell = m.cell)

    @printf("\n== %s: model cell at (%.3f, %.3f), %.0f m deep; ADCP %s; common grid %g-%g m every %g m, %d hours\n",
            name, m.cell[3], m.cell[4], m.cell[5], dataset, grid[1], grid[end], Δz, length(axis))
    upper = findall(z -> z <= min(grid[1] + 60, grid[end]), grid)
    KE(x) = [nanmean(0.5 .* (x.U[upper, n] .^ 2 .+ x.V[upper, n] .^ 2)) for n in eachindex(axis)]
    ko, km = KE(obs_ni), KE(mdl_ni)
    @printf("   near-inertial KE, mean over %g-%g m and the record: ADCP %.2e, model %.2e J/kg (model/ADCP %.2f)\n",
            grid[upper[1]], grid[upper[end]], nanmean(ko), nanmean(km), nanmean(km) / nanmean(ko))
    for (label, a, b) in STORMS
        w = findall(t -> a <= t < b, axis)
        @printf("   %-17s near-inertial KE: ADCP %.2e, model %.2e (ratio %.2f)\n", label, nanmean(ko[w]), nanmean(km[w]),
                nanmean(km[w]) / nanmean(ko[w]))
    end
    println("   squared shear (10⁻⁵ s⁻²) by depth: total hourly, ADCP / model;  near-inertial, ADCP / model;  same in storm windows")
    So, Sm = shear²(obs.U, obs.V, Δz), shear²(mdl.U, mdl.V, Δz)
    No, Nm = shear²(obs_ni.U, obs_ni.V, Δz), shear²(mdl_ni.U, mdl_ni.V, Δz)
    storm = findall(t -> any(a <= t < b for (_, a, b) in STORMS), axis)
    for k in axes(So, 1)
        k > 12 && break
        @printf("   %5.0f-%-5.0f m %7.2f / %-7.2f %7.2f / %-7.2f   storms %7.2f / %-7.2f %7.2f / %-7.2f\n", grid[k], grid[k+1],
                1e5 * nanmean(So[k, :]), 1e5 * nanmean(Sm[k, :]), 1e5 * nanmean(No[k, :]), 1e5 * nanmean(Nm[k, :]),
                1e5 * nanmean(So[k, storm]), 1e5 * nanmean(Sm[k, storm]), 1e5 * nanmean(No[k, storm]), 1e5 * nanmean(Nm[k, storm]))
    end
    results[name] = merge(results[name], (; ko, km, So, Sm, No, Nm))
end

# ---------------- figure ----------------
mnames = sort(collect(keys(results)))
fig = Figure(size = (1700, 380 * length(mnames)), fontsize = 14)
Label(fig[0, 1:3], "$(basename(PREFIX)) vs OOI Pioneer current profilers: near-inertial ($(Int(BAND[1]))-$(Int(BAND[2])) h) currents and shear",
      fontsize = 18)
for (row, name) in enumerate(mnames)
    r = results[name]
    tdays = [Dates.value(t - r.axis[1]) / 86_400_000 for t in r.axis]
    ax = Axis(fig[row, 1], title = "$name: near-inertial KE, $(r.grid[1])-$(min(r.grid[1] + 60, r.grid[end])) m (J/kg)",
              xlabel = "days since $(Dates.format(r.axis[1], "yyyy-mm-dd"))")
    lines!(ax, tdays, r.ko; color = :black, label = "ADCP"); lines!(ax, tdays, r.km; color = :firebrick, label = "model")
    row == 1 && axislegend(ax; position = :lt)
    zmid = (r.grid[1:end-1] .+ r.grid[2:end]) ./ 2
    for (c, (A, B, t)) in enumerate(((r.So, r.Sm, "total hourly"), (r.No, r.Nm, "near-inertial")))
        ax = Axis(fig[row, 1 + c], title = "$name: mean squared shear, $t (s⁻²)", xscale = log10, yreversed = true, ylabel = "depth (m)")
        lines!(ax, max.(rowmean(A), 1e-9), zmid; color = :black, label = "ADCP")
        lines!(ax, max.(rowmean(B), 1e-9), zmid; color = :firebrick, label = "model")
    end
end
out = PREFIX * "_moorings_vs_ooi.png"
save(out, fig; px_per_unit = 1.2)
println("saved ", out)
