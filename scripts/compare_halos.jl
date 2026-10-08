# Compare every cell, halos included, of the state written with MAB_STAGE=steps between a single-rank run and a split
# run: each split rank's cell is matched to the single-rank cell at the same whole-domain index (local index + offset).
# Reports, per field and step, the largest difference among the rank's own cells, its halo columns (x), its halo rows (y)
# and the corner cells, with where it is.
#   julia --project=. scripts/compare_halos.jl /path/s_1 /path/s_x  [steps, default 0,1]
using NumericalEarth, Printf
const jldopen = NumericalEarth.DataWrangling.jldopen

ref, other = ARGS[1], ARGS[2]
steps = length(ARGS) >= 3 ? parse.(Int, split(ARGS[3], ",")) : [0, 1]
rankfiles(tag, n) = sort(filter(f -> occursin(Regex("^" * basename(tag) * "_state" * string(n) * "_rank\\d+\\.jld2\$"), basename(f)), readdir(dirname(tag); join = true)))

for n in steps
    single = jldopen(only(rankfiles(ref, n)), "r") do io; Dict(k => io[k] for k in keys(io)); end
    for file in rankfiles(other, n)
        r = jldopen(file, "r") do io; Dict(k => io[k] for k in keys(io)); end
        io_, jo = r["I_OFF"], r["J_OFF"]
        println("== step $n, ", basename(file), "  (offsets ", io_, ", ", jo, ")")
        names = [chop(k; tail = 2) for k in keys(r) if endswith(k, "_p") && haskey(single, k)]
        for name in names
            A = single[name * "_p"]; HA = single[name * "_halo"]
            B = r[name * "_p"];      HB = r[name * "_halo"]
            ref = haskey(r, name) ? name : "eta"                          # interior extent of this rank (centre count for 2D fields)
            nx = size(r[ref], 1)
            ny = size(r[ref], 2)
            worst = Dict{String, Any}(c => (0.0, ()) for c in ("own", "xhalo", "yhalo", "corner"))
            for k in axes(B, 3), q in axes(B, 2), p in axes(B, 1)
                i = p - HB[1]; j = q - HB[2]                # local index (may be ≤ 0 or > n in the halo)
                gi = i + io_; gj = j + jo                   # whole-domain index
                pa = gi + HA[1]; qa = gj + HA[2]
                (1 <= pa <= size(A, 1) && 1 <= qa <= size(A, 2)) || continue
                d = abs(B[p, q, k] - A[pa, qa, k]); isfinite(d) || continue
                inx = 1 <= i <= nx; iny = 1 <= j <= ny
                cls = inx && iny ? "own" : (!inx && iny ? "xhalo" : (inx && !iny ? "yhalo" : "corner"))
                d > worst[cls][1] && (worst[cls] = (d, (i, j, k, gi, gj)))
            end
            @printf("  %-7s halo %s  ", name, string(Tuple(HB)))
            for cls in ("own", "xhalo", "yhalo", "corner")
                d, at = worst[cls]
                @printf("%s %.2e%s   ", cls, d, d > 0 ? " at local(i,j,k)=" * string(at[1:3]) * " global(i,j)=" * string(at[4:5]) : "")
            end
            println()
        end
    end
end
