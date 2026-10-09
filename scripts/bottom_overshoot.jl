# A run's saltiest bottom cell, day by day, from the de-tided daily output: salinity and temperature in that cell, in
# the cell above, in the bottom cells of its four neighbours, and GLORYS's bottom salinity and temperature there; the
# currents (u, v, w) at the cell; and the full column on the worst day against GLORYS. Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/bottom_test/b12_partial MAB_START_DATE=2019-08-29 \
#       julia --project=/t0/workdir/enrique/mpi05_ib scripts/bottom_overshoot.jl
# MAB_CELL: "i,j" to follow (default: the cell holding the largest salinity over the whole run)
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const PREFIX = run_prefix(get(ENV, "MAB_TAG", "mab_H90"))
volume = PREFIX * "_volume_daily.jld2"
S = open_series(volume, "S"; backend = OnDisk()); T = open_series(volume, "T"; backend = OnDisk())
U = open_series(volume, "u"; backend = OnDisk()); V = open_series(volume, "v"; backend = OnDisk())
W = open_series(volume, "w"; backend = OnDisk())
ug = S.grid.underlying_grid
Nx, Ny, Nz = size(ug)
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center())); zc = collect(znodes(ug, Center()))
days = S.times ./ 86400

# the bottom (deepest active) cell of a column: inactive cells hold zero
bottom(s, i, j) = findfirst(!iszero, view(s, i, j, :))

frames = [Array(interior(S[n])) for n in eachindex(days)]
i, j = if haskey(ENV, "MAB_CELL")
    parse.(Int, split(ENV["MAB_CELL"], ","))
else
    best = argmax([maximum(f) for f in frames])
    Tuple(argmax(frames[best]))[1:2]
end
kb = bottom(frames[end], i, j)
@printf("%s: cell (%d, %d) at %.3f°E %.3f°N, bottom cell k = %d (z = %.1f m), %d active levels\n",
        basename(PREFIX), i, j, λ[i], φ[j], kb, zc[kb], Nz - kb + 1)

println("day   S_bottom  S_above  S_neighbours(W,E,S,N)        GLORYS_S  T_bottom  GLORYS_T   u      v      w(×1e4)")
for n in eachindex(days)
    s = frames[n]; t = Array(interior(T[n]))
    u = Array(interior(U[n])); v = Array(interior(V[n])); w = Array(interior(W[n]))
    nb = [(a, b) for (a, b) in ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1))]
    sn = [1 ≤ a ≤ Nx && 1 ≤ b ≤ Ny && !isnothing(bottom(s, a, b)) ? s[a, b, bottom(s, a, b)] : NaN for (a, b) in nb]
    gs = days[n] ≥ 1 ? glorys_at("so", days[n]) : nothing
    gt = days[n] ≥ 1 ? glorys_at("thetao", days[n]) : nothing
    # GLORYS at the cell, at the model's bottom-cell depth (nearest level above GLORYS's own bottom there)
    function gval(g)
        isnothing(g) && return NaN
        lon, lat, dep, A = g
        prof = [bilinear(view(A, :, :, k), lon, lat, λ[i], φ[j]) for k in eachindex(dep)]
        good = findall(isfinite, prof)
        isempty(good) && return NaN
        k = good[argmin(abs.(dep[good] .+ zc[kb]))]
        return prof[k]
    end
    @printf("%4.1f  %7.3f  %7.3f  %6.2f %6.2f %6.2f %6.2f   %7.3f  %7.3f  %7.3f  %+6.3f %+6.3f %+6.2f\n",
            days[n], s[i, j, kb], s[i, j, kb + 1], sn..., gval(gs), t[i, j, kb], gval(gt),
            (u[i, j, kb] + u[i + 1, j, kb]) / 2, (v[i, j, kb] + v[i, j + 1, kb]) / 2, 1e4 * w[i, j, kb + 1])
end

n = argmax([f[i, j, kb] for f in frames])
println("\ncolumn on day $(days[n]) (k, z, S, T):")
s = frames[n]; t = Array(interior(T[n]))
for k in Nz:-1:kb
    @printf("  %3d %8.1f  %7.3f  %7.3f\n", k, zc[k], s[i, j, k], t[i, j, k])
end
