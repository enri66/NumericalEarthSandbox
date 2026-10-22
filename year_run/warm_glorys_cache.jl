# Build the inpainted GLORYS caches for the year-long run on a few independent processes (each date and variable is written by exactly
# one process). Script 05's ranks all read these files; if they are missing the 25 ranks inpaint and replace them at the same time,
# and a rank that has one open (mapped) when it is replaced dies with a bus error on NFS.
#   WARM_PART=i WARM_PARTS=n julia --project=... warm_glorys_cache.jl
using NumericalEarth, Oceananigans, Dates
using NumericalEarth.DataWrangling: Metadata, Metadatum, BoundingBox
using CopernicusMarine
dir = joinpath(homedir(), "Data", "NumericalEarth")
region = BoundingBox(longitude = (-76.0, -64.0), latitude = (34.0, 42.0))
part, parts = parse(Int, get(ENV, "WARM_PART", "1")), parse(Int, get(ENV, "WARM_PARTS", "1"))
dates = collect(DateTime(2019, 3, 31):Day(1):DateTime(2020, 4, 2))
jobs = [(name, date) for date in dates for name in (:u_velocity, :v_velocity, :temperature, :salinity, :free_surface)]
mine = jobs[part:parts:end]
t₀ = time()
for (n, (name, date)) in enumerate(mine)
    Field(Metadatum(name; dataset = GLORYSDaily(), date, dir, region))
    n % 50 == 0 && println("part $part: $n of $(length(mine)) done after $(round(Int, time() - t₀)) s")
end
println("part $part done: $(length(mine)) fields in $(round(Int, time() - t₀)) s")
