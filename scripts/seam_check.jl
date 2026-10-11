# Do a run's fields jump at the edges of its MPI ranks? For one frame of the raw surface output `<tag>.jld2`, the mean
# over wet points of the absolute second difference across each column (along x) and each row (along y), for u, v, T and
# the Rossby number ζ/f; the largest ones are listed, marked when they lie on a rank edge. Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/res_test/hr100_fx3 MAB_START_DATE=2019-08-29 \
#       julia --project=/t0/workdir/enrique/mpi05_ib scripts/seam_check.jl
# MAB_FRAME: the frame (default: the middle one)
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const PREFIX = run_prefix(get(ENV, "MAB_TAG", "mab_H90"))
path = PREFIX * ".jld2"
blocks = rank_blocks(path)
grid = global_grid(path); ug = grid.underlying_grid
Nx, Ny = ug.Nx, ug.Ny
B = grid.immersed_boundary.bottom_height
wet = [B[i, j, 1] < 0 for i in 1:Nx, j in 1:Ny]
px, py = size(blocks)
xedges = Set(k * Nx ÷ px for k in 1:px-1); yedges = Set(k * Ny ÷ py for k in 1:py-1)
println("$(px) × $(py) ranks; rank edges after columns $(sort(collect(xedges))) and rows $(sort(collect(yedges)))")

T, times = surface_series(path, "T"); u, _ = surface_series(path, "u"); v, _ = surface_series(path, "v")
n = haskey(ENV, "MAB_FRAME") ? parse(Int, ENV["MAB_FRAME"]) : length(times) ÷ 2
@printf("frame %d, day %.2f\n", n, times[n] / 86400)

φc = collect(φnodes(ug, Center())); φf = collect(φnodes(ug, Face()))
Δλ = deg2rad(λnodes(ug, Center())[2] - λnodes(ug, Center())[1]); Δφ = deg2rad(φc[2] - φc[1])
U = u[:, :, n]; V = v[:, :, n]
Ro = fill(NaN, Nx, Ny)
for i in 2:Nx, j in 2:Ny
    (wet[i-1, j-1] & wet[i, j-1] & wet[i-1, j] & wet[i, j]) || continue
    ζ = ((V[i, j] - V[i-1, j]) / Δλ - (U[i, j] * cosd(φc[j]) - U[i, j-1] * cosd(φc[j-1])) / Δφ) / (6.371e6 * cosd(φf[j]))
    Ro[i, j] = ζ / (2 * 7.292115e-5 * sind(φf[j]))
end

# fields on their own grids, with the wet mask of the cells they sit in
fields = (("u", U[1:Nx, :]), ("v", V[:, 1:Ny]), ("T", T[:, :, n]), ("ζ/f", Ro))
function score_columns(A)
    s = fill(NaN, Nx)
    for i in 2:Nx-1
        d = [abs(A[i+1, j] - 2A[i, j] + A[i-1, j]) for j in 1:Ny if wet[i-1, j] & wet[i, j] & wet[i+1, j] && all(isfinite, (A[i-1, j], A[i, j], A[i+1, j]))]
        isempty(d) || (s[i] = mean(d))
    end
    return s
end
function score_rows_proper(A)
    s = fill(NaN, Ny)
    for j in 2:Ny-1
        d = [abs(A[i, j+1] - 2A[i, j] + A[i, j-1]) for i in 1:Nx if wet[i, j-1] & wet[i, j] & wet[i, j+1] && all(isfinite, (A[i, j-1], A[i, j], A[i, j+1]))]
        isempty(d) || (s[j] = mean(d))
    end
    return s
end

for (name, A) in fields
    sc = score_columns(A); sr = score_rows_proper(A)
    med_c = median(filter(isfinite, sc)); med_r = median(filter(isfinite, sr))
    top_c = sort(filter(i -> isfinite(sc[i]), 2:Nx-1); by = i -> -sc[i])[1:8]
    top_r = sort(filter(j -> isfinite(sr[j]), 2:Ny-1); by = j -> -sr[j])[1:8]
    println("\n$name: median column score $(round(med_c; sigdigits = 3)), median row score $(round(med_r; sigdigits = 3))")
    println("  worst columns: ", join(["$i$(i in xedges || i - 1 in xedges ? "*" : "") ($(round(sc[i] / med_c; digits = 1))×)" for i in top_c], ", "))
    println("  worst rows:    ", join(["$j$(j in yedges || j - 1 in yedges ? "*" : "") ($(round(sr[j] / med_r; digits = 1))×)" for j in top_r], ", "))
    println("  at the rank edges (columns/rows next to an edge, × median): columns ",
            join(["$(e):$(round(sc[e] / med_c; digits = 1))/$(round(sc[e+1] / med_c; digits = 1))" for e in sort(collect(xedges))], " "),
            "; rows ", join(["$(e):$(round(sr[e] / med_r; digits = 1))/$(round(sr[e+1] / med_r; digits = 1))" for e in sort(collect(yedges))], " "))
end

# Short streaks: runs of at least MAB_STREAK_LEN points along a column (or a row) whose second difference across it is
# above MAB_STREAK_K times the field's median there, and what u, v and T do at the same points
const LEN = parse(Int, get(ENV, "MAB_STREAK_LEN", "8"))
const K = parse(Float64, get(ENV, "MAB_STREAK_K", "8"))
d2x(A, i, j) = abs(A[i+1, j] - 2A[i, j] + A[i-1, j])
d2y(A, i, j) = abs(A[i, j+1] - 2A[i, j] + A[i, j-1])
Uc = U[1:Nx, :]; Vc = V[:, 1:Ny]; Tn = T[:, :, n]
med = median(filter(isfinite, [d2x(Ro, i, j) for i in 2:Nx-1, j in 2:Ny]))
streaks = []
for i in 3:Nx-1
    run = Int[]
    for j in 2:Ny
        hot = all(isfinite, (Ro[i-1, j], Ro[i, j], Ro[i+1, j])) && d2x(Ro, i, j) > K * med
        if hot
            push!(run, j)
        else
            length(run) ≥ LEN && push!(streaks, (i, first(run), last(run)))
            run = Int[]
        end
    end
    length(run) ≥ LEN && push!(streaks, (i, first(run), last(run)))
end
λc = collect(λnodes(ug, Center()))
println("\n$(length(streaks)) column streaks in ζ/f (≥ $LEN points, second difference across x > $K × median)")
for (i, j1, j2) in sort(streaks; by = s -> -(s[3] - s[2]))[1:min(end, 15)]
    js = j1:j2
    ru = mean(d2x(Uc, i, j) for j in js) / median(filter(isfinite, [d2x(Uc, a, b) for a in 2:Nx-1, b in js]))
    rv = mean(d2x(Vc, i, j) for j in js) / median(filter(isfinite, [d2x(Vc, a, b) for a in 2:Nx-1, b in js]))
    rt = mean(d2x(Tn, i, j) for j in js if wet[i-1, j] & wet[i, j] & wet[i+1, j]; init = NaN) / median(filter(isfinite, [d2x(Tn, a, b) for a in 2:Nx-1, b in js]))
    @printf("  column %3d (%.3f°E, rank column %d, %s) rows %3d-%3d (%.2f-%.2f°N): across-x jump vs median: u %.1f×, v %.1f×, T %.1f×; depth %.0f m\n",
            i, λc[i], (i - 1) ÷ (Nx ÷ px) + 1, (i in xedges || i - 1 in xedges) ? "RANK EDGE" : "interior", j1, j2, φc[j1], φc[j2], ru, rv, rt, -B[i, (j1 + j2) ÷ 2, 1])
end
