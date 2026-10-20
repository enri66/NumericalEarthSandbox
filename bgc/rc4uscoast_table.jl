# RC4USCoast (river chemistry, NCEI accession 0260455) monthly climatology 1990-2022: the rivers whose mouth is in or next to the MAB box.
# julia --project=<env with NumericalEarth> bgc/rc4uscoast_table.jl <dir with mclim_19902022_chem.nc and mclim_19902022_disc.nc>
using NumericalEarth, Printf, Statistics
const NCD = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "NCDatasets")]
dir = ARGS[1]
chem = NCD.NCDataset(joinpath(dir, "mclim_19902022_chem.nc")); disc = NCD.NCDataset(joinpath(dir, "mclim_19902022_disc.nc"))
println("discharge file variables: ", join(keys(disc), ", "))
name = String.(chem["river_name"][:]); mlat = Float64.(chem["mouth_lat"][:]); mlon = Float64.(chem["mouth_lon"][:])
nanmean(v) = (w = filter(isfinite, v); isempty(w) ? NaN : mean(w))
get_(v, i, m) = (x = chem[v][i, m]; ismissing(x) || !isfinite(x) ? NaN : Float64(x))
dvar = first(filter(k -> !(k in ("time", "time_bnds", "climatology_bnds", "RC4USCoast_ID", "original_ID", "lat", "lon", "mouth_lat", "mouth_lon", "region", "river_name", "source", "dataset")), keys(disc)))
Q(i, m) = (x = disc[dvar][i, m]; ismissing(x) || !isfinite(x) ? NaN : Float64(x))
sel = [i for i in eachindex(name) if -76.3 < mlon[i] < -64 && 34.0 < mlat[i] < 42.5]
@printf("%-18s %7s %7s | %6s | Feb-Mar mean: %6s %6s %6s %6s %6s %7s %7s %6s | annual NO3 %5s\n", "river", "mouthλ", "mouthφ", "Q m3/s", "NO3", "NOx", "NH4", "DON", "PO4", "Alk", "DIC", "DO", "")
for i in sel
    fm(v) = nanmean([get_(v, i, 2), get_(v, i, 3)])
    q = nanmean([Q(i, 2), Q(i, 3)])
    ann = nanmean([get_("no3", i, m) for m in 1:12])
    @printf("%-18s %7.2f %7.2f | %6.0f | %27.1f %6.1f %6.1f %6.1f %6.2f %7.0f %7.0f %6.0f | %5.1f\n", name[i], mlon[i], mlat[i], q, fm("no3"), fm("nox"), fm("nh4"), fm("don"), fm("po4"), fm("alk"), fm("dic"), fm("do"), ann)
end
println("\nmonthly NO3 (mmol/m3) Jan..Dec, Delaware / Connecticut / Raritan / Susquehanna:")
for r in ("Delaware", "Connecticut", "Raritan", "Susquehanna")
    i = findfirst(==(r), name); println(rpad(r, 12), join([@sprintf("%5.1f", get_("no3", i, m)) for m in 1:12], " "))
end
