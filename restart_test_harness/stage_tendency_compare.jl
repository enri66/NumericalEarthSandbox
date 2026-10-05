# julia stage_tendency_compare.jl before after  -> bitwise comparison of every stored array
using Serialization
A = deserialize("/tmp/stage_tendency_$(ARGS[1]).jls")
B = deserialize("/tmp/stage_tendency_$(ARGS[2]).jls")
for configuration in sort(collect(keys(A)))
    a, b = A[configuration], B[configuration]
    differing = [k for k in keys(a) if !isequal(a[k], b[k])]
    worst = isempty(differing) ? 0.0 : maximum(maximum(abs, a[k] .- b[k]) for k in differing)
    println(rpad(configuration, 26), isempty(differing) ? "bit-identical ($(length(a)) arrays)" :
                                                          "DIFFERS in $(join(sort(differing), ", ")), max |Δ| = $worst")
end
