using NumericalEarth.DataWrangling: jldopen
maxdiff = 0.0; nvals = 0
for r in 0:3
    a = jldopen("/t0/workdir/enrique/year_run/test_sponge/old_rank$(r).jld2"); b = jldopen("/t0/workdir/enrique/year_run/test_sponge/new_rank$(r).jld2")
    for name in keys(a["timeseries"])
        name == "t" && continue
        its = filter(k -> all(isdigit, k), keys(a["timeseries"][name])); itsb = filter(k -> all(isdigit, k), keys(b["timeseries"][name]))
        its == itsb || (println("rank $r $name: different iterations"); continue)
        for it in its
            x = a["timeseries"][name][it]; y = b["timeseries"][name][it]
            x isa AbstractArray || (println("skipping $name ($(typeof(x)))"); break)
            d = maximum(abs.(replace(x .- y, NaN => 0.0))); global maxdiff = max(maxdiff, d); global nvals += length(x)
            d > 0 && println("rank $r $name it $it: max |old - new| = $d")
        end
    end
    close(a); close(b)
end
println("compared $nvals values; max |old - new| = $maxdiff")
