# Where does a run's salinity peak, and is the place a shallow cell or a one-cell bay? For each frame of the raw surface
# output, and for each de-tided daily full-depth frame, the maximum salinity, its cell, the depth there and the number of
# wet neighbours (of 4); then the cells above a threshold at the last frame, and the same cells in a reference run.
# Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/bottom_test/b12_pz MAB_REFERENCE=/t0/workdir/enrique/runs/bottom_test/b12_partial \
#   MAB_START_DATE=2019-08-29 julia --project=/t0/workdir/enrique/mpi05_ib scripts/salinity_peak.jl
# MAB_SMAX: the threshold (psu, default 37.2: above the GLORYS maximum in the box)
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const PREFIX = run_prefix(get(ENV, "MAB_TAG", "mab_H90"))
const REFERENCE = get(ENV, "MAB_REFERENCE", "")
const SMAX = parse(Float64, get(ENV, "MAB_SMAX", "37.2"))

path = PREFIX * ".jld2"
grid = global_grid(path); ug = grid.underlying_grid
Nx, Ny, _ = size(ug)
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center()))
B = grid.immersed_boundary.bottom_height
bh = [B[i, j, 1] for i in 1:Nx, j in 1:Ny]
wet = bh .< 0
neighbours(i, j) = count(((a, b),) -> 1 ≤ a ≤ Nx && 1 ≤ b ≤ Ny && wet[a, b], ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)))
describe(i, j) = @sprintf("(%d, %d) %.3f°E %.3f°N, depth %.1f m, %d wet neighbours", i, j, λ[i], φ[j], -bh[i, j], neighbours(i, j))

S, times = surface_series(path, "S")
S[.!repeat(wet, 1, 1, size(S, 3))] .= NaN
println("surface salinity maximum, every 8th frame of $(basename(path)):")
for n in 1:8:length(times)
    s = S[:, :, n]
    c = argmax(replace(s, NaN => -Inf))
    @printf("  %s  %.2f psu at %s\n", Dates.format(start_date + Second(round(Int, times[n])), "mm-dd HH:MM"), s[c], describe(c.I...))
end

last = S[:, :, end]
hot = [(i, j) for i in 1:Nx, j in 1:Ny if wet[i, j] && last[i, j] > SMAX]
@printf("\n%d surface cells above %.1f psu at the last frame (%s):\n", length(hot), SMAX,
        Dates.format(start_date + Second(round(Int, times[end])), "yyyy-mm-dd HH:MM"))
Sref = isempty(REFERENCE) ? nothing : surface_series(run_prefix(REFERENCE) * ".jld2", "S")[1][:, :, end]
for (i, j) in sort(hot; by = c -> -last[c...])[1:min(end, 15)]
    @printf("  %.2f psu at %s%s\n", last[i, j], describe(i, j),
            isnothing(Sref) ? "" : @sprintf("; %s: %.2f psu", basename(REFERENCE), Sref[i, j]))
end
