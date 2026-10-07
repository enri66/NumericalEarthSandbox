# Compare the state written by script 05 with MAB_STAGE=steps for the same case run on one rank and split over ranks:
# for every dumped step and field, the largest difference, where it is, and how many cells differ at all. Each field is
# assembled from the ranks' blocks with the offsets stored in the files (I_OFF, J_OFF).
#   julia --project=. scripts/compare_layouts.jl /path/s_1 /path/s_x /path/s_y     (the run tags)
using NumericalEarth, Printf
const jldopen = NumericalEarth.DataWrangling.jldopen

function assemble(tag, n, name)
    files = filter(f -> occursin(Regex("^" * basename(tag) * "_state" * string(n) * "_rank\\d+\\.jld2\$"), basename(f)), readdir(dirname(tag); join = true))
    blocks = [jldopen(f, "r") do io; (io["I_OFF"], io["J_OFF"], io[name]); end for f in files]
    Nx = maximum(b[1] + size(b[3], 1) for b in blocks); Ny = maximum(b[2] + size(b[3], 2) for b in blocks); Nz = size(blocks[1][3], 3)
    A = fill(NaN, Nx, Ny, Nz)
    for (io, jo, a) in blocks; A[io+1:io+size(a, 1), jo+1:jo+size(a, 2), :] .= a; end
    return A
end

tags = ARGS
ref = tags[1]
steps = sort(unique(parse(Int, match(r"_state(\d+)_rank", f).captures[1]) for f in readdir(dirname(ref)) if startswith(f, basename(ref) * "_state")))
for other in tags[2:end]
    println("== ", basename(ref), " vs ", basename(other))
    for n in steps, name in ("eta", "u", "v", "T", "S", "e")
        A = assemble(ref, n, name); B = assemble(other, n, name)
        size(A) == size(B) || (println("  step $n $name: shapes differ ", size(A), " ", size(B)); continue)
        D = abs.(A .- B); D[.!isfinite.(D)] .= 0
        m, idx = findmax(D)
        @printf("  step %2d %-4s max|Δ| = %.3e at (i,j,k) = %s   cells differing: %d of %d\n", n, name, m, string(Tuple(idx)), count(>(0), D), length(D))
    end
end
