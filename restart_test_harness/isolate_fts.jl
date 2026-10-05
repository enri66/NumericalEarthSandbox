# Isolate whether the bug is in the raw FieldTimeSeries/InMemory backend, or in NumericalEarth's atmosphere
# layer around it. Two ways of reaching t=300 on the same on-disk series with a window of 3:
#  "continuous": one FTS object, queried in sequence at 0,60,...,240,300 (like the uninterrupted run)
#  "cold":       a FRESH FTS object (like a restart rebuilding everything), queried directly at 240 then 300
using Oceananigans, Oceananigans.Units
using Oceananigans.OutputReaders: OnDisk, InMemory

grid = LatitudeLongitudeGrid(size = (12, 10, 8), halo = (7, 7, 7),
                              longitude = (-71, -70), latitude = (37, 38), z = (-100, 0))
path = joinpath(mktempdir(), "T_series.jld2")
times = 0:60:600.0
let f_tmp = Field{Center, Center, Nothing}(grid)
    f = FieldTimeSeries{Center, Center, Nothing}(grid, times; backend = OnDisk(), path, name = "T")
    for (n, t) in enumerate(times)
        set!(f_tmp, 15 + t / 60)
        set!(f, f_tmp, n)
    end
end

println("=== continuous access (one FTS object, queried in order) ===")
fts_cont = FieldTimeSeries(path, "T"; backend = InMemory(3), architecture = CPU())
v_cont = 0.0
for t in (0.0, 60.0, 120.0, 180.0, 240.0, 300.0)
    v_cont = interior(fts_cont[Time(t)])[1,1,1]
    println("t=$t  T=$v_cont  backend=", fts_cont.backend)
end

println("\n=== cold access (fresh FTS object, jump straight to 240 then 300) ===")
fts_cold = FieldTimeSeries(path, "T"; backend = InMemory(3), architecture = CPU())
v240 = interior(fts_cold[Time(240.0)])[1,1,1]
println("t=240.0  T=$v240  backend=", fts_cold.backend)
v300 = interior(fts_cold[Time(300.0)])[1,1,1]
println("t=300.0  T=$v300  backend=", fts_cold.backend)

println("\n################ raw FieldTimeSeries comparison at t=300")
println("continuous: T=$v_cont   cold: T=$v300   Δ=$(abs(v_cont - v300))")
