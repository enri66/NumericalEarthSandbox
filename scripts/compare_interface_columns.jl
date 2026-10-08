# For an x-split run against the single-rank run (MAB_STAGE=steps dumps with halos), the largest difference in each
# local column i near the rank interfaces, for every dumped field: rank r's last own columns and first halo columns,
# and its first own columns and last halo columns on the other side.
#   julia --project=. scripts/compare_interface_columns.jl /path/s_1 /path/s_x  [steps, default 0,1]  [width, default 3]
using NumericalEarth, Printf
const jldopen = NumericalEarth.DataWrangling.jldopen

ref, other = ARGS[1], ARGS[2]
steps = length(ARGS) >= 3 ? parse.(Int, split(ARGS[3], ",")) : [0, 1]
width = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 3
rankfiles(tag, n) = sort(filter(f -> occursin(Regex("^" * basename(tag) * "_state" * string(n) * "_rank\\d+\\.jld2\$"), basename(f)), readdir(dirname(tag); join = true)))
load(file) = jldopen(file, "r") do io; Dict(k => io[k] for k in keys(io)); end

for n in steps
    single = load(only(rankfiles(ref, n)))
    for file in rankfiles(other, n)
        r = load(file)
        io_, jo = r["I_OFF"], r["J_OFF"]
        nx = size(r["eta"], 1)
        cols = vcat(collect(1 - width:width), collect(nx - width + 1:nx + width))
        println("== step $n, ", basename(file), "  (offset ", io_, ")   local columns ", join(cols, " "))
        for name in sort([chop(k; tail = 2) for k in keys(r) if endswith(k, "_p") && haskey(single, k)])
            A = single[name * "_p"]; HA = single[name * "_halo"]
            B = r[name * "_p"];      HB = r[name * "_halo"]
            line = IOBuffer()
            for i in cols
                p = i + HB[1]; pa = i + io_ + HA[1]
                if !(1 <= p <= size(B, 1) && 1 <= pa <= size(A, 1))
                    print(line, "       -  ")
                    continue
                end
                d = 0.0
                for k in axes(B, 3), q in axes(B, 2)
                    j = q - HB[2]; qa = j + jo + HA[2]
                    1 <= qa <= size(A, 2) || continue
                    δ = abs(B[p, q, k] - A[pa, qa, k]); isfinite(δ) && (d = max(d, δ))
                end
                print(line, d == 0 ? "       0  " : @sprintf("%8.1e  ", d))
            end
            @printf("  %-8s %s\n", name, String(take!(line)))
        end
    end
end
