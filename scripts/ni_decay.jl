# How fast does near-inertial energy decay after a storm, in the model runs and at the OOI Pioneer current profilers?
# For each mooring and each source (ADCP, and every run in MAB_TAGS that saved MAB_MOORINGS=pioneer columns), the
# near-inertial kinetic energy is computed as in moorings_vs_ooi.jl (16-22 h band-pass on a common depth grid,
# averaged over the upper 60 m of the usable range), smoothed with a one-inertial-period running mean to remove the
# 2f ripple, and after each storm an exponential is fitted to it from its peak (searched over the storm window) over
# the next FIT_DAYS days: KE ∝ exp(-t / τ). Reported: peak KE and τ (days).
# Usage:
#   MAB_TAGS=/t0/.../catke_fall/moor,/t0/.../res_test/r12,/t0/.../res_test/r24 OOI_DIR=/t0/workdir/enrique/Data/OOI/pioneer \
#       julia --project=. scripts/ni_decay.jl
include(joinpath(@__DIR__, "mab_analysis_common.jl"))
const FFTW = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "FFTW")]

const TAGS     = split(get(ENV, "MAB_TAGS", ""), ",")
const OOI_DIR  = get(ENV, "OOI_DIR", joinpath(homedir(), "Data", "OOI", "pioneer"))
const FIT_DAYS = parse(Float64, get(ENV, "FIT_DAYS", "4"))
const BAND     = (16.0, 22.0)
const STORMS   = (("Dorian", DateTime(2019, 9, 6), DateTime(2019, 9, 10)),
                  ("mid-Sep", DateTime(2019, 9, 14), DateTime(2019, 9, 18)),
                  ("17 Oct", DateTime(2019, 10, 16), DateTime(2019, 10, 20)))
const MOORINGS = Dict("OSSM" => ("ooi-cp04ossm-mfd35-01-adcpsj000", 8.0, 40.0),
                      "PMUO" => ("ooi-cp02pmuo-rii01-02-adcpsl010", 8.0, 40.0),
                      "PMCO" => ("ooi-cp02pmco-rii01-02-adcptg010", 4.0, 12.0),
                      "CNSM" => ("ooi-cp01cnsm-mfd35-01-adcptf000", 4.0, 12.0))
const T_INERTIAL = 2π / (2 * 7.292e-5 * sind(40.0)) / 3600     # hours

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

# Near-inertial KE (smoothed) on an hourly axis
function ni_series(src, grid, axis)
    ti = Dict(t => n for (n, t) in enumerate(src.times))
    take(A) = [haskey(ti, t) ? A[m, ti[t]] : NaN for m in axes(A, 1), t in axis]
    U = take(regrid(src.depths, src.U, grid)); V = take(regrid(src.depths, src.V, grid))
    Ub = reduce(vcat, [bandpass(U[m, :])' for m in axes(U, 1)]); Vb = reduce(vcat, [bandpass(V[m, :])' for m in axes(V, 1)])
    ke = [nanmean(0.5 .* (Ub[:, n] .^ 2 .+ Vb[:, n] .^ 2)) for n in eachindex(axis)]
    w = round(Int, T_INERTIAL ÷ 2)
    return [nanmean(ke[max(1, n - w):min(end, n + w)]) for n in eachindex(ke)]
end

# Peak in the storm window, then a least-squares line through log KE over FIT_DAYS days
function decay(axis, ke, a, b)
    win = findall(t -> a <= t <= b, axis); isempty(win) && return (NaN, NaN)
    p = win[argmax(ke[win])]
    fit = p:min(length(ke), p + round(Int, 24FIT_DAYS))
    t = (fit .- p) ./ 24; y = log.(ke[fit])
    ok = isfinite.(y); count(ok) < 24 && return (ke[p], NaN)
    slope = cov(t[ok], y[ok]) / var(t[ok])
    return ke[p], -1 / slope
end

println("Near-inertial KE decay after storms: peak KE (10⁻³ J/kg) and e-folding time τ (days) over $(FIT_DAYS) days after the peak")
for name in sort(collect(keys(MOORINGS)))
    dataset, Δz, ztop = MOORINGS[name]
    adcp = read_adcp(dataset)
    grid = collect(ztop:Δz:ztop + 60)
    axis = collect(DateTime(2019, 8, 30):Hour(1):DateTime(2019, 10, 27))
    sources = [("ADCP", adcp); [(basename(t), read_model(String(t), name)) for t in TAGS if !isempty(t)]]
    println("\n== $name ($(grid[1])-$(grid[end]) m)")
    @printf("   %-10s", "source"); for (s, _, _) in STORMS; @printf(" %22s", s); end; println()
    for (label, src) in sources
        isnothing(src) && continue
        ke = ni_series(src, grid, axis)
        @printf("   %-10s", label)
        for (_, a, b) in STORMS
            k, τ = decay(axis, ke, a, b)
            @printf("   %6.2f  τ = %5.1f d", 1e3k, τ)
        end
        println()
    end
end
