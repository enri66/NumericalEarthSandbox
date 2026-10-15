# Sea ice of a GSL run (script 05 with MAB_REGION=gsl, MAB_SEA_ICE=true) against GLORYS, from its checkpoints: ice area
# (Σ ℵ A) and volume (Σ ℵ h A) in the Gulf of St. Lawrence, on the Labrador side north of 51.7°N (the Labrador shelf,
# Lake Melville) and on the Newfoundland shelf and Grand Banks (the rest), leaving out the sponge band (8 cells from the
# box edges), where the
# ice is restored to GLORYS. Also where the surface salinity is lowest. Writes <tag>_ice_maps.jld2 (model and GLORYS ℵ, h
# and surface S on each checkpoint date) for plotting elsewhere.
#
#   julia --project=<env> gsl_ice_vs_glorys.jl <run dir> <tag> [data dir]
using NumericalEarth
using NumericalEarth.DataWrangling: jldopen, BoundingBox
using NumericalEarth.DataWrangling.GLORYS: GLORYSDaily
using Dates, Printf
const NCD = NumericalEarth.DataWrangling.NCDatasets

rundir, tag = ARGS[1], ARGS[2]
datadir = length(ARGS) ≥ 3 ? ARGS[3] : joinpath(homedir(), "Data", "NumericalEarth")
region = BoundingBox(longitude = (-71.5, -48.0), latitude = (45.0, 56.0))
start = DateTime(2013, 12, 1)

const Nλ, Nφ, H = 278, 128, 7
λc = [-71.5 + 2 / 12 + (i - 0.5) / 12 for i in 1:Nλ]
φc = [45.0 + 2 / 12 + (j - 0.5) / 12 for j in 1:Nφ]
area = [(6.371e6 * deg2rad(1 / 12))^2 * cosd(φc[j]) for i in 1:Nλ, j in 1:Nφ]
inner = [8 < i ≤ Nλ - 8 && 8 < j ≤ Nφ - 8 for i in 1:Nλ, j in 1:Nφ]
labrador = [φc[j] ≥ 51.7 for i in 1:Nλ, j in 1:Nφ]
gulf = [!labrador[i, j] && ((λc[i] < -59 && φc[j] > 46) || (-59 ≤ λc[i] < -56.5 && 47.5 < φc[j])) for i in 1:Nλ, j in 1:Nφ]
shelf = .!labrador .& .!gulf

function glorys(var, name, date)
    file = joinpath(datadir, NumericalEarth.DataWrangling.metadata_filename(GLORYSDaily(), name, date, region))
    NCD.Dataset(file) do ds
        data = Float64.(coalesce.(ds[var][:, :, 1], 0.0))
        λg, φg = Float64.(ds["longitude"][:]), Float64.(ds["latitude"][:])
        ig = [argmin(abs.(λg .- λ)) for λ in λc]; jg = [argmin(abs.(φg .- φ)) for φ in φc]
        [data[ig[i], jg[j]] for i in 1:Nλ, j in 1:Nφ]
    end
end

const Nz = 50
# the interior of a 2D field, or the surface level of a 3D one
interior2d(a) = a[H+1:H+Nλ, H+1:H+Nφ, size(a, 3) == 1 ? 1 : H + Nz]

files = sort(filter(f -> startswith(f, tag * "_checkpoint_rank0_iteration") && endswith(f, ".jld2"), readdir(rundir)),
             by = f -> parse(Int, match(r"iteration(\d+)", f)[1]))
maps = Dict{String, Any}()
@printf("%-10s %-8s %12s %12s %12s %12s\n", "date", "where", "area model", "area GLORYS", "vol model", "vol GLORYS")
for f in files
    ckpt = jldopen(joinpath(rundir, f))
    t = ckpt["simulation/model/clock"].time
    h = interior2d(ckpt["simulation/model/sea_ice/model/ice_thickness/data"])
    ℵ = interior2d(ckpt["simulation/model/sea_ice/model/ice_concentration/data"])
    S = interior2d(ckpt["simulation/model/ocean/model/tracers/S/data"])
    T = interior2d(ckpt["simulation/model/ocean/model/tracers/T/data"])
    close(ckpt)
    date = start + Second(round(Int, t))
    ℵg, hg = glorys("siconc", :sea_ice_concentration, date), glorys("sithick", :sea_ice_thickness, date)
    for (label, m) in (("gulf", inner .& gulf), ("labrador", inner .& labrador), ("shelf", inner .& shelf))
        @printf("%-10s %-8s %9.1f k² %9.1f k² %9.1f km³ %9.1f km³\n", Dates.format(date, "yyyy-mm-dd"), label,
                sum(ℵ .* area .* m) / 1e9, sum(ℵg .* area .* m) / 1e9,
                sum(ℵ .* h .* area .* m) / 1e9, sum(ℵg .* hg .* area .* m) / 1e9)
    end
    wet = S .> 0 .|| T .!= 0
    Smin, k = findmin(s -> s, [wet[i] ? S[i] : Inf for i in eachindex(S)])
    i, j = Tuple(CartesianIndices(S)[k])
    @printf("           lowest surface salinity %.2f at (i, j) = (%d, %d), %.2f°E %.2f°N; %d surface cells below 20\n",
            Smin, i, j, λc[i], φc[j], count(s -> 0 ≤ s < 20, S[wet]))
    maps[Dates.format(date, "yyyy-mm-dd")] = (; ℵ, h, S, T, ℵg, hg)
end
jldopen(joinpath(rundir, tag * "_ice_maps.jld2"), "w") do file
    file["λ"] = λc; file["φ"] = φc
    for (d, m) in maps, (k, v) in pairs(m)
        file["$d/$k"] = v
    end
end
