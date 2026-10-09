# Sea-surface height of several runs against GLORYS, from the de-tided daily output: the area-weighted domain mean of η
# (the volume budget, set by the open boundaries) and, for the pattern, the rms difference and the correlation of the
# anomaly from each one's own mean, over the shelf (< 200 m) and the deep water, day by day. Usage:
#   RUNS=/t0/workdir/enrique/runs/bottom_test/b12_pz,/t0/workdir/enrique/runs/bottom_test/b12_pz2 MAB_START_DATE=2019-08-29 \
#       julia --project=/t0/workdir/enrique/mpi05_ib scripts/ssh_compare.jl
# MAB_DAYS: days to report (default: every 2nd day available in all runs)
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const RUNS = split(ENV["RUNS"], ",")
names = basename.(RUNS)

grid = global_grid(RUNS[1] * "_eta_daily.jld2"); ug = grid.underlying_grid
Nx, Ny, _ = size(ug)
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center()))
B = grid.immersed_boundary.bottom_height
bh = [B[i, j, 1] for i in 1:Nx, j in 1:Ny]
wet = bh .< 0
area = [cosd(φ[j]) for i in 1:Nx, j in 1:Ny]                     # cell area up to a constant (uniform λ, φ spacing)
shelf = wet .& (-bh .< 200); deep = wet .& (-bh .> 1000)

series = [surface_series(r * "_eta_daily.jld2", "η") for r in RUNS]
common_days = reduce(intersect, [round.(Int, t ./ 86400) for (_, t) in series])
days = haskey(ENV, "MAB_DAYS") ? parse.(Int, split(ENV["MAB_DAYS"], ",")) : common_days[1:2:end]

wmean(X, m) = sum(X[m] .* area[m]) / sum(area[m])
function glorys_eta(d)
    day0 = start_date + Day(d) - Day(1)
    read(date) = begin
        ds = NCD.Dataset(gfile("zos", date))
        lon = Float64.(ds["longitude"][:]); lat = Float64.(ds["latitude"][:])
        A = [ismissing(x) ? NaN : Float64(x) for x in ds["zos"][:, :, 1]]
        close(ds); (lon, lat, A)
    end
    lon, lat, A0 = read(day0); _, _, A1 = read(day0 + Day(1))
    A = (A0 .+ A1) ./ 2
    return [wet[i, j] ? bilinear(A, lon, lat, λ[i], φ[j]) : NaN for i in 1:Nx, j in 1:Ny]
end
stats(X, G, m) = begin
    good = m .& isfinite.(X) .& isfinite.(G)
    x = X[good] .- wmean(X, good); g = G[good] .- wmean(G, good)
    (sqrt(mean((x .- g) .^ 2)), cor(x, g))
end

@printf("%-5s %-14s %10s %10s %10s %8s %10s %8s\n", "day", "run", "mean η", "GLORYS", "mean diff", "", "", "")
@printf("%-5s %-14s %10s %10s %10s %8s %10s %8s\n", "", "", "(m)", "(m)", "(m)", "shelf rms", "shelf r", "deep rms")
for d in days
    G = glorys_eta(d)
    for (r, (η, t)) in enumerate(series)
        n = argmin(abs.(t ./ 86400 .- d))
        X = ifelse.(wet, η[:, :, n], NaN)
        m = wmean(X, wet .& isfinite.(G)); g = wmean(G, wet .& isfinite.(G))
        e_s, c_s = stats(X, G, shelf); e_d, _ = stats(X, G, deep)
        @printf("%-5d %-14s %+10.3f %+10.3f %+10.3f %8.3f %10.2f %8.3f\n", d, names[r], m, g, m - g, e_s, c_s, e_d)
    end
end
