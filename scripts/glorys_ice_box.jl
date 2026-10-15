# GLORYS sea-ice area and volume in a box, and the ice fluxes across its north and east edges, every 15 days through
# a winter. Reads the global daily GLORYS12 archive on antares (A3). Used to pick the sea-ice region and period
# (HANDOFF_SEA_ICE.md). Run on antares with an environment that has NCDatasets.
using NCDatasets, Dates
const R = "/Volumes/A3/enrique/GLORYS12"
λb, φb = (-70.5, -48.0), (45.0, 56.0)
fileof(d) = (dir = joinpath(R, Dates.format(d, "yyyy"), Dates.format(d, "mm")); tag = Dates.format(d, "yyyymmdd");
             f = filter(x -> occursin("_mean_$tag", x), readdir(dir)); joinpath(dir, f[1]))
ds0 = NCDataset(fileof(Date(2013,12,1))); lon = ds0["longitude"][:]; lat = ds0["latitude"][:]; close(ds0)
ii = findall(x -> λb[1] <= x <= λb[2], lon); jj = findall(x -> φb[1] <= x <= φb[2], lat)
Re = 6.371e6; d = deg2rad(1/12); A = [Re^2*d*d*cosd(lat[j]) for i in ii, j in jj]; dx = [Re*d*cosd(lat[j]) for j in jj]
println("date        area(10^3 km2)  vol(km3)  | north edge: ice area flux southward (10^3 km2/day), vol flux (km3/day) | east edge: westward area flux")
for dt in vcat(Date(2013,11,15), Date(2013,12,1):Day(15):Date(2014,5,31))
    NCDataset(fileof(dt)) do ds
        c = coalesce.(ds["siconc"][ii, jj, 1], 0.0); h = coalesce.(ds["sithick"][ii, jj, 1], 0.0)
        v = coalesce.(ds["vsi"][ii, jj, 1], 0.0); u = coalesce.(ds["usi"][ii, jj, 1], 0.0)
        jn = length(jj)
        aflux_n = -sum(c[:, jn] .* v[:, jn]) * dx[jn] * 86400 / 1e9
        vflux_n = -sum(c[:, jn] .* h[:, jn] .* v[:, jn]) * dx[jn] * 86400 / 1e9
        aflux_e = -sum(c[end, :] .* u[end, :] .* (Re*d)) * 86400 / 1e9
        println(dt, "   ", rpad(round(sum(c .* A)/1e9, digits=1), 8), "  ", rpad(round(sum(c .* h .* A)/1e9, digits=1), 8),
                "  |  ", rpad(round(aflux_n, digits=2), 6), "  ", round(vflux_n, digits=2), "  |  ", round(aflux_e, digits=2))
    end
end
