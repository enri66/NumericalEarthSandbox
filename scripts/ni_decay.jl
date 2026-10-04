# How fast does near-inertial energy decay after a storm, in the model runs and at the OOI Pioneer current profilers?
# For each mooring and each source (ADCP, and every run in MAB_TAGS that saved MAB_MOORINGS=pioneer columns), the
# near-inertial kinetic energy is computed as in moorings_vs_ooi.jl (16-22 h band-pass on a common depth grid,
# averaged over the upper 60 m of the usable range), smoothed with a one-inertial-period running mean to remove the
# 2f ripple, and after each storm an exponential is fitted to it from its peak (searched over the storm window) over
# the next FIT_DAYS days: KE ∝ exp(-t / τ). Reported: peak KE and τ (days).
# Usage:
#   MAB_TAGS=/t0/.../catke_fall/moor,/t0/.../res_test/r12,/t0/.../res_test/r24 OOI_DIR=/t0/workdir/enrique/Data/OOI/pioneer \
#       julia --project=. scripts/ni_decay.jl
include(joinpath(@__DIR__, "ooi_common.jl"))

const TAGS     = split(get(ENV, "MAB_TAGS", ""), ",")
const FIT_DAYS = parse(Float64, get(ENV, "FIT_DAYS", "4"))
const STORMS   = (("Dorian", DateTime(2019, 9, 6), DateTime(2019, 9, 10)),
                  ("mid-Sep", DateTime(2019, 9, 14), DateTime(2019, 9, 18)),
                  ("17 Oct", DateTime(2019, 10, 16), DateTime(2019, 10, 20)))

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
