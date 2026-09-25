# Model vs GLORYS vertical sections of temperature and salinity, with the mixed-layer depth and the thermocline, and
# maps of both depths, from a run of 04_mab_glorys_tides_reservoirs.jl. Uses the de-tided daily full-depth output
# (`<tag>_volume_daily.jld2`), which carries its own grid, so it works for any vertical grid. Usage:
#   MAB_TAG=mab_H90 MAB_DAYS=30,60,87 julia --project=. scripts/transects_vs_glorys.jl
# MAB_DAYS:     days to analyse (default: the last daily frame)
# MAB_SECTIONS: "lat:<φ>" (zonal) or "lon:<λ>" (meridional) sections, comma separated
# MAB_ZMAX:     depth shown in the sections (m)
#
# Mixed-layer depth: where σ₀ first exceeds its 10 m value by 0.03 kg/m³ (de Boyer Montégut et al. 2004).
# Seasonal thermocline depth: the midpoint of the largest temperature decrease between two adjacent levels, between the
# mixed-layer base and 300 m, where that decrease is at least 0.02 °C/m (otherwise there is none).
using Oceananigans, NumericalEarth, Printf, Statistics, Dates, CairoMakie
using Oceananigans.Grids: znodes, λnodes, φnodes
const NCD = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "NCDatasets")]           # loaded by NumericalEarth
const SWP = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "SeawaterPolynomials")]  # loaded by Oceananigans

const OUT  = joinpath(homedir(), "Data", "mab_glorys_obc")
const GDIR = get(ENV, "MAB_DATA_DIR", joinpath(homedir(), "Data", "NumericalEarth"))
const TAG  = get(ENV, "MAB_TAG", "mab_H90")
const start_date = DateTime(2019, 4, 1)
const ZMAX = parse(Float64, get(ENV, "MAB_ZMAX", "500"))
const SECTIONS = [(Symbol(split(s, ":")[1]), parse(Float64, split(s, ":")[2]))
                  for s in split(get(ENV, "MAB_SECTIONS", "lat:35.0,lat:38.0,lat:40.5,lon:-70.5"), ",")]

const eos = SWP.TEOS10.TEOS10EquationOfState()
σ₀(T, S) = SWP.ρ(T, S, 0, eos) - 1000

# ---------------- model ----------------
vol = joinpath(OUT, TAG * "_volume_daily.jld2")
TV = FieldTimeSeries(vol, "T"; backend = OnDisk()); SV = FieldTimeSeries(vol, "S"; backend = OnDisk())
grid = TV.grid; ug = grid.underlying_grid
Nx, Ny, Nz = size(ug)
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center())); zc = collect(znodes(ug, Center()))
B = grid.immersed_boundary.bottom_height
bh = [B[i, j, 1] for i in 1:Nx, j in 1:Ny]
wet = bh .< 0
wet3 = [zc[k] > bh[i, j] for i in 1:Nx, j in 1:Ny, k in 1:Nz]
tdays = TV.times ./ 86400
days = haskey(ENV, "MAB_DAYS") ? parse.(Float64, split(ENV["MAB_DAYS"], ",")) : [tdays[end]]

# ---------------- GLORYS ----------------
gfile(var, date) = first(filter(x -> occursin("$(var)_GLORYSDaily_$(Dates.format(date, "yyyy-mm-dd"))T", x) &&
                                     endswith(x, "-76.0_-64.0_34.0_42.0.nc"), readdir(GDIR; join = true)))

function read_glorys(var, date)
    ds = NCD.Dataset(gfile(var, date))
    lon = Float64.(ds["longitude"][:]); lat = Float64.(ds["latitude"][:]); dep = Float64.(ds["depth"][:])
    A = [ismissing(x) ? NaN : Float64(x) for x in ds[var][:, :, :, 1]]
    close(ds)
    return lon, lat, dep, A
end

# A de-tided frame is centred on midnight, between two GLORYS daily means (centred on noon)
function glorys_at(var, d)
    day0 = start_date + Day(floor(Int, d)) - Day(1)
    lon, lat, dep, A = read_glorys(var, day0)
    _, _, _, A1 = read_glorys(var, day0 + Day(1))
    return lon, lat, dep, (A .+ A1) ./ 2
end

function bilinear(A, lon, lat, x, y)
    i = clamp(searchsortedlast(lon, x), 1, length(lon) - 1); j = clamp(searchsortedlast(lat, y), 1, length(lat) - 1)
    a = (x - lon[i]) / (lon[i+1] - lon[i]); b = (y - lat[j]) / (lat[j+1] - lat[j])
    w = ((1 - a) * (1 - b), a * (1 - b), (1 - a) * b, a * b); v = (A[i, j], A[i+1, j], A[i, j+1], A[i+1, j+1])
    num = sum(w[n] * v[n] for n in 1:4 if isfinite(v[n]); init = 0.0)
    den = sum(w[n] for n in 1:4 if isfinite(v[n]); init = 0.0)
    return den > 0 ? num / den : NaN
end

# GLORYS profile (native levels) at a model column; levels below the model's bottom are dropped
glorys_profile(A, lon, lat, dep, i, j) =
    [dep[k] < -bh[i, j] ? bilinear(view(A, :, :, k), lon, lat, λ[i], φ[j]) : NaN for k in eachindex(dep)]

# ---------------- mixed layer and thermocline from one profile ----------------
# depth: positive downward and increasing; T, S: NaN below the bottom
function column_depths(depth, T, S)
    ok = isfinite.(T) .& isfinite.(S)
    d = depth[ok]; t = T[ok]; s = S[ok]
    (length(d) < 3 || d[1] > 10 || d[end] < 10) && return (NaN, NaN)
    σ = σ₀.(t, s)
    m = searchsortedlast(d, 10.0)
    σref = m == length(d) ? σ[m] : σ[m] + (σ[m+1] - σ[m]) * (10 - d[m]) / (d[m+1] - d[m])
    mld = d[end]                                        # whole column mixed: the bottom
    for k in m+1:length(d)
        if σ[k] > σref + 0.03
            mld = d[k-1] + (d[k] - d[k-1]) * (σref + 0.03 - σ[k-1]) / (σ[k] - σ[k-1])
            mld = max(mld, 10.0)
            break
        end
    end
    best = 0.02; tcl = NaN
    for k in 1:length(d)-1
        (d[k] >= mld && d[k+1] <= 300) || continue
        g = (t[k] - t[k+1]) / (d[k+1] - d[k])
        g > best && (best = g; tcl = (d[k] + d[k+1]) / 2)
    end
    return mld, tcl
end

# ---------------- analysis ----------------
region(i, j) = -bh[i, j] < 200 ? 1 : -bh[i, j] < 1000 ? 2 : 3
region_names = ("shelf (<200 m)", "slope (200-1000 m)", "deep (>1000 m)")

for dwant in days
    n = argmin(abs.(tdays .- dwant)); d = tdays[n]
    @printf("\n==== %s, day %.1f (%s) ====\n", TAG, d, Dates.format(start_date + Second(round(Int, d * 86400)), "yyyy-mm-dd HH:MM"))
    T = Array(interior(TV[n])); S = Array(interior(SV[n]))
    T[.!wet3] .= NaN; S[.!wet3] .= NaN
    lon, lat, dep, GT = glorys_at("thetao", d)
    _, _, _, GS = glorys_at("so", d)
    depth_m = reverse(-zc)                              # model depths, top first

    # maps of mixed-layer and thermocline depth
    mld = fill(NaN, Nx, Ny, 2); tcl = fill(NaN, Nx, Ny, 2)
    for i in 1:Nx, j in 1:Ny
        wet[i, j] || continue
        mld[i, j, 1], tcl[i, j, 1] = column_depths(depth_m, reverse(T[i, j, :]), reverse(S[i, j, :]))
        mld[i, j, 2], tcl[i, j, 2] = column_depths(dep, glorys_profile(GT, lon, lat, dep, i, j), glorys_profile(GS, lon, lat, dep, i, j))
    end
    println("                         mixed-layer depth (m)                 seasonal thermocline (m)")
    println("region                   model  GLORYS  bias   rms   cells     model  GLORYS  bias   rms   cells")
    for r in 1:3
        cm = [(i, j) for i in 1:Nx, j in 1:Ny if wet[i, j] && region(i, j) == r && isfinite(mld[i, j, 1]) && isfinite(mld[i, j, 2])]
        ct = [(i, j) for i in 1:Nx, j in 1:Ny if wet[i, j] && region(i, j) == r && isfinite(tcl[i, j, 1]) && isfinite(tcl[i, j, 2])]
        a = [mld[i, j, 1] for (i, j) in cm]; b = [mld[i, j, 2] for (i, j) in cm]
        c = [tcl[i, j, 1] for (i, j) in ct]; e = [tcl[i, j, 2] for (i, j) in ct]
        stat(x, y) = isempty(x) ? (NaN, NaN, NaN, NaN) : (median(x), median(y), mean(x .- y), sqrt(mean((x .- y) .^ 2)))
        @printf("%-22s %6.1f %6.1f %+6.1f %6.1f %6d   %6.1f %6.1f %+6.1f %6.1f %6d\n", region_names[r],
                stat(a, b)..., length(cm), stat(c, e)..., length(ct))
    end

    fig = Figure(size = (1500, 900), fontsize = 16)
    Label(fig[0, 1:3], @sprintf("%s day %.1f: mixed-layer depth (top) and seasonal thermocline depth (bottom)", TAG, d), fontsize = 20)
    for (row, (A, name, cr)) in enumerate(((mld, "MLD", (0, 100)), (tcl, "thermocline", (0, 150))))
        for (col, (X, t)) in enumerate(((A[:, :, 1], "model"), (A[:, :, 2], "GLORYS"), (A[:, :, 1] .- A[:, :, 2], "model − GLORYS")))
            ax = Axis(fig[2row - 1, col], title = "$(name) (m): $(t)", aspect = DataAspect())
            hm = heatmap!(ax, λ, φ, X; colormap = col == 3 ? :balance : :deep, colorrange = col == 3 ? (-cr[2] / 2, cr[2] / 2) : cr,
                          nan_color = :gray85)
            for (kind, v) in SECTIONS
                kind === :lat ? hlines!(ax, v; color = :black, linewidth = 1) : vlines!(ax, v; color = :black, linewidth = 1)
            end
            (col == 1 || col == 3) && Colorbar(fig[2row, col == 1 ? (1:2) : 3], hm; vertical = false, flipaxis = false)
        end
    end
    out = joinpath(OUT, @sprintf("%s_mld_thermocline_day%02d.png", TAG, round(Int, d)))
    save(out, fig; px_per_unit = 1.2); println("saved ", out)

    # sections
    fig = Figure(size = (1900, 420 * length(SECTIONS)), fontsize = 15)
    Label(fig[0, 1:5], @sprintf("%s day %.1f: temperature and salinity sections, model vs GLORYS (white: mixed layer, dashed: thermocline)", TAG, d), fontsize = 19)
    println("\nsection                upper $(round(Int, ZMAX)) m RMS model − GLORYS:  T (°C)   S")
    for (row, (kind, v)) in enumerate(SECTIONS)
        if kind === :lat
            j = argmin(abs.(φ .- v)); cols = [(i, j) for i in 1:Nx]; x = λ; xl = "longitude"; lbl = @sprintf("%.2f°N", φ[j])
        else
            i = argmin(abs.(λ .- v)); cols = [(i, j) for j in 1:Ny]; x = φ; xl = "latitude"; lbl = @sprintf("%.2f°E", λ[i])
        end
        Tm = [T[i, j, k] for (i, j) in cols, k in 1:Nz]; Sm = [S[i, j, k] for (i, j) in cols, k in 1:Nz]
        kg = findall(<=(ZMAX + 500), dep)
        Tg = fill(NaN, length(cols), length(kg)); Sg = fill(NaN, length(cols), length(kg))
        for (c, (i, j)) in enumerate(cols)
            wet[i, j] || continue
            Tg[c, :] = glorys_profile(GT, lon, lat, dep, i, j)[kg]; Sg[c, :] = glorys_profile(GS, lon, lat, dep, i, j)[kg]
        end
        # RMS over the upper ZMAX m, GLORYS interpolated linearly in depth to the model levels
        dT = Float64[]; dS = Float64[]
        for (c, (i, j)) in enumerate(cols), k in 1:Nz
            (isfinite(Tm[c, k]) && -zc[k] <= ZMAX) || continue
            m = clamp(searchsortedlast(dep[kg], -zc[k]), 1, length(kg) - 1); f = clamp((-zc[k] - dep[m]) / (dep[m+1] - dep[m]), 0, 1)
            tg = (1 - f) * Tg[c, m] + f * Tg[c, m+1]; sg = (1 - f) * Sg[c, m] + f * Sg[c, m+1]
            isfinite(tg) && push!(dT, Tm[c, k] - tg); isfinite(sg) && push!(dS, Sm[c, k] - sg)
        end
        @printf("%-6s %-15s %30.2f %8.3f\n", String(kind), lbl, sqrt(mean(dT .^ 2)), sqrt(mean(dS .^ 2)))
        mline = [mld[i, j, 1] for (i, j) in cols]; gline = [mld[i, j, 2] for (i, j) in cols]
        mtl = [tcl[i, j, 1] for (i, j) in cols]; gtl = [tcl[i, j, 2] for (i, j) in cols]
        floor_line = [wet[i, j] ? bh[i, j] : 0.0 for (i, j) in cols]
        panels = ((Tm, zc, mline, mtl, "T model", :thermal, (5, 28)), (Tg, -dep[kg], gline, gtl, "T GLORYS", :thermal, (5, 28)),
                  (Sm, zc, mline, mtl, "S model", :haline, (31, 36.5)), (Sg, -dep[kg], gline, gtl, "S GLORYS", :haline, (31, 36.5)))
        for (col, (A, zz, ml, tl, name, cm, cr)) in enumerate(panels)
            ax = Axis(fig[row, col], title = "$(name), $(lbl)", xlabel = xl, ylabel = col == 1 ? "z (m)" : "")
            p = sortperm(zz)
            hm = heatmap!(ax, x, zz[p], A[:, p]; colormap = cm, colorrange = cr, nan_color = :gray85)
            lines!(ax, x, -ml; color = :white, linewidth = 2)
            lines!(ax, x, -tl; color = :white, linewidth = 1.5, linestyle = :dash)
            lines!(ax, x, floor_line; color = :black, linewidth = 1)
            ylims!(ax, -ZMAX, 0)
            (col == 2 || col == 4) && Colorbar(fig[row, 5 + (col == 4)], hm; label = col == 2 ? "°C" : "")
        end
    end
    out = joinpath(OUT, @sprintf("%s_sections_day%02d.png", TAG, round(Int, d)))
    save(out, fig; px_per_unit = 1.2); println("saved ", out)
end
