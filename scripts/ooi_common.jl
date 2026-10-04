# Shared pieces of the OOI Pioneer mooring analyses (ni_decay.jl, ni_budget.jl): the current profilers' hourly
# profiles (ERDDAP CSV in OOI_DIR) and a run's hourly mooring columns (MAB_MOORINGS=pioneer), regridding to a common
# depth grid, gap filling and the near-inertial band-pass.
include(joinpath(@__DIR__, "mab_analysis_common.jl"))
const FFTW = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "FFTW")]

const OOI_DIR = get(ENV, "OOI_DIR", joinpath(homedir(), "Data", "OOI", "pioneer"))
const BAND    = (16.0, 22.0)                                  # near-inertial periods (hours); f at 40°N: 18.6 h
const T_INERTIAL = 2π / (2 * 7.292e-5 * sind(40.0)) / 3600     # hours
# name in a run's mooring output => (ERDDAP dataset, common-grid spacing, top of the usable ADCP range), in metres
const MOORINGS = Dict("OSSM" => ("ooi-cp04ossm-mfd35-01-adcpsj000", 8.0, 40.0),
                      "PMUO" => ("ooi-cp02pmuo-rii01-02-adcpsl010", 8.0, 40.0),
                      "PMCO" => ("ooi-cp02pmco-rii01-02-adcptg010", 4.0, 12.0),
                      "CNSM" => ("ooi-cp01cnsm-mfd35-01-adcptf000", 4.0, 12.0))

function read_adcp(dataset)
    lines = readlines(joinpath(OOI_DIR, dataset * ".csv"))
    col = Dict(h => i for (i, h) in enumerate(split(lines[1], ",")))
    qc = haskey(col, "eastward_sea_water_velocity_qc_agg")
    data = Dict{Tuple{DateTime, Float64}, NTuple{2, Float64}}()
    for line in lines[3:end]
        v = split(line, ",")
        qc && !(v[col["eastward_sea_water_velocity_qc_agg"]] in ("1", "2")) && continue
        u = tryparse(Float64, v[col["eastward_sea_water_velocity"]]); w = tryparse(Float64, v[col["northward_sea_water_velocity"]])
        (isnothing(u) || isnothing(w) || abs(u) > 3 || abs(w) > 3) && continue
        d = -parse(Float64, v[col["z"]]); d > 0 || continue
        data[(floor(DateTime(v[col["time"]][1:19]), Hour), d)] = (u, w)
    end
    times = sort(unique(first.(keys(data)))); depths = sort(unique(last.(keys(data))))
    U = fill(NaN, length(depths), length(times)); V = fill(NaN, length(depths), length(times))
    ti = Dict(t => n for (n, t) in enumerate(times)); di = Dict(d => n for (n, d) in enumerate(depths))
    for ((t, d), (u, w)) in data
        U[di[d], ti[t]] = u; V[di[d], ti[t]] = w
    end
    return (; times, depths, U, V)
end

function read_model(prefix, name)
    for f in filter(f -> occursin(r"_moorings_rank\d+\.jld2$", f) && startswith(basename(f), basename(prefix) * "_moorings"),
                    readdir(dirname(prefix); join = true))
        out = JLD2.jldopen(f) do file
            haskey(file, name) || return nothing
            t0 = DateTime(file["start_date"])
            times = [round(t0 + Second(round(Int, s)), Hour) for s in file["$name/time"]]
            keep = [true; diff(times) .> Hour(0)]
            (; times = times[keep], depths = reverse(-file["z_center"]), U = reverse(file["$name/u"][:, keep]; dims = 1),
               V = reverse(file["$name/v"][:, keep]; dims = 1))
        end
        isnothing(out) || return out
    end
    return nothing
end

function regrid(depths, A, grid)
    out = fill(NaN, length(grid), size(A, 2))
    for n in axes(A, 2)
        ok = isfinite.(A[:, n]) .& (A[:, n] .!= 0); d = depths[ok]; a = A[ok, n]
        length(d) < 2 && continue
        for (m, z) in enumerate(grid)
            d[1] <= z <= d[end] || continue
            k = clamp(searchsortedlast(d, z), 1, length(d) - 1)
            out[m, n] = a[k] + (a[k+1] - a[k]) * (z - d[k]) / (d[k+1] - d[k])
        end
    end
    return out
end
function fill_gaps(x)
    ok = findall(isfinite, x); length(ok) < length(x) ÷ 2 && return fill(NaN, length(x))
    y = copy(x)
    for n in eachindex(y)
        isfinite(y[n]) && continue
        a = findlast(<(n), ok); b = findfirst(>(n), ok)
        y[n] = isnothing(a) ? x[ok[b]] : isnothing(b) ? x[ok[a]] : x[ok[a]] + (x[ok[b]] - x[ok[a]]) * (n - ok[a]) / (ok[b] - ok[a])
    end
    return y
end
function bandpass(x)
    y = fill_gaps(x); any(isnan, y) && return y
    Y = FFTW.fft(y .- sum(y) / length(y)); fr = FFTW.fftfreq(length(y), 1.0)
    Y[.!(1 / BAND[2] .<= abs.(fr) .<= 1 / BAND[1])] .= 0
    return real(FFTW.ifft(Y))
end
nanmean(x) = (f = filter(isfinite, x); isempty(f) ? NaN : sum(f) / length(f))
