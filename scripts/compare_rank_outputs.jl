# Compare a script-05 run on one rank with the same run split over PX × 1 ranks: for every output file kind
# (<tag>_<kind>_rank0.jld2 …), every saved field and every saved iteration, the ranks' interiors are joined along x
# (halos stripped; they can differ in width between the builds) and compared with the one-rank field. Prints the largest difference and where it is per field (last iteration) and
# the first iteration with any difference.
#
#   julia --project=<env> compare_rank_outputs.jl <dir> <one-rank tag> <split tag> <PX> <Nx> <Ny>
using NumericalEarth.DataWrangling: jldopen   # JLD2 through NumericalEarth (not a direct dependency of the run environments)
using Printf

dir, tag1, tagN = ARGS[1:3]
PX, Nx, Ny = parse.(Int, ARGS[4:6])

# The interior of a saved array of a field with N interior points (+1 for a face on a Bounded side) along x and y
function interior_of(a, N, Ny)
    Hx, Hy = (size(a, 1) - N) ÷ 2, (size(a, 2) - Ny) ÷ 2
    return a[Hx+1:size(a, 1)-Hx, Hy+1:size(a, 2)-Hy, :]
end

kinds = sort(unique([replace(f, tag1 * "_" => "", r"_?rank0\.jld2$" => "") for f in readdir(dir)
                     if startswith(f, tag1 * "_") && endswith(f, "rank0.jld2") && !occursin("checkpoint", f)]))

for kind in kinds
    name(t, r) = joinpath(dir, isempty(kind) ? "$(t)_rank$(r).jld2" : "$(t)_$(kind)_rank$(r).jld2")
    isfile(name(tagN, 0)) || continue
    one = jldopen(name(tag1, 0))
    split = [jldopen(name(tagN, r)) for r in 0:PX-1]
    haskey(one, "timeseries") || (close(one); foreach(close, split); continue)
    fields = [k for k in keys(one["timeseries"]) if k != "t"]
    println("== $(kind == "" ? "(main)" : kind)")
    for f in fields
        iterations = sort(parse.(Int, [k for k in keys(one["timeseries/$f"]) if all(isdigit, k)]))
        first_diff, last_line = nothing, ""
        for it in iterations
            a₀ = one["timeseries/$f/$it"]
            ex, ey = (size(a₀, 1) - Nx) % 2, (size(a₀, 2) - Ny) % 2   # 1 for a face along x (y) on a Bounded side
            a = interior_of(a₀, Nx + ex, Ny + ey)
            b = cat([interior_of(s["timeseries/$f/$it"], Nx ÷ PX + (r == PX ? ex : 0), Ny + ey)
                     for (r, s) in enumerate(split)]...; dims = 1)
            if size(a) != size(b)
                last_line = "size $(size(a)) vs $(size(b))"
                break
            end
            d = abs.(a .- b)
            Δ, at = findmax(map(x -> isnan(x) ? Inf : x, d))
            Δ > 0 && isnothing(first_diff) && (first_diff = it)
            last_line = @sprintf("max|Δ| %.2e of %.2e at %s", Δ, maximum(abs, a), Tuple(at))
        end
        @printf("   %-22s %s  (last of %d iterations; first difference at %s)\n", f, last_line, length(iterations),
                isnothing(first_diff) ? "none" : string(first_diff))
    end
    close(one); foreach(close, split)
end
