# Prints the class table of river_chemistry.jl: julia --project=<env with NumericalEarth> bgc/river_chemistry_test.jl <RC4USCoast dir>
using NumericalEarth, Printf
include(joinpath(@__DIR__, "river_chemistry.jl"))
table, rivers = river_chemistry_table(ARGS[1]; box = (-75.83, -64.17, 34.17, 41.83))
for (class, list) in rivers; println(rpad(class, 9), "rivers: ", join(list, ", ")); end
println("\nmonth   class      ", join([rpad(t, 8) for t in keys(RC4_VARIABLES)]))
for m in (1, 2, 3, 4, 7), class in (:Delaware, :Hudson, :other)
    println(rpad(m, 8), rpad(class, 11), join([@sprintf("%-8.1f", table[t][class][m]) for t in keys(RC4_VARIABLES)]))
end
