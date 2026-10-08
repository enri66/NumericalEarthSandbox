# Is the model's mixed-layer depth bias (model − GLORYS) concentrated at fronts? If so, a mixed-layer eddy
# restratification parameterization (Fox-Kemper et al. 2008, 2011) is the right fix; if the bias is uniform, it is not.
# For each column: model and GLORYS MLD (as in transects_vs_glorys.jl), the model's horizontal buoyancy gradient
# averaged over its own mixed layer |∇b|, and the Fox-Kemper restratifying buoyancy flux scale H²|∇b|²/|f|.
#   MAB_TAG=/abs/path/prefix MAB_DAYS=20,30 julia --project=<run env> scripts/mld_bias_vs_fronts.jl
using Oceananigans, NumericalEarth, Printf, Statistics, Dates, CairoMakie
using Oceananigans.Grids: znodes, λnodes, φnodes
const NCD = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "NCDatasets")]
const SWP = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "SeawaterPolynomials")]

const GDIR = get(ENV, "MAB_DATA_DIR", joinpath(homedir(), "Data", "NumericalEarth"))
const PREFIX = ENV["MAB_TAG"]
const start_date = DateTime(get(ENV, "MAB_START_DATE", "2019-04-01"))   # the run's start (its clock's zero)
const eos = SWP.TEOS10.TEOS10EquationOfState()
σ₀(T, S) = SWP.ρ(T, S, 0, eos) - 1000
const g, ρ₀, Ω, Rₑ = 9.81, 1026.0, 7.292e-5, 6.371e6

# ---- model output, combining per-rank files (see transects_vs_glorys.jl) ----
function global_grid(vol)
    JLD2 = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "JLD2")]
    stem = basename(splitext(vol)[1])
    paths = filter(f -> occursin(Regex("^" * stem * "_rank\\d+\\.jld2\$"), basename(f)), readdir(dirname(vol); join = true))
    grids = [JLD2.jldopen(f -> f["serialized/grid"], p) for p in paths]
    ugs = [gr.underlying_grid for gr in grids]
    order = sortperm([ug.λᶠᵃᵃ[1] for ug in ugs])
    bottom(r) = (b = grids[r].immersed_boundary.bottom_height; n = ugs[r].Nx; m = ugs[r].Ny;
                 ndims(b) == 2 ? collect(b[1:n, 1:m]) : collect(b[1:n, 1:m, 1]))
    bh = cat([bottom(r) for r in order]...; dims = 1)
    u1, u2 = ugs[order[1]], ugs[order[end]]
    gr = LatitudeLongitudeGrid(CPU(); size = (size(bh, 1), u1.Ny, u1.Nz), longitude = (u1.λᶠᵃᵃ[1], u2.λᶠᵃᵃ[u2.Nx + 1]),
                               latitude = (u1.φᵃᶠᵃ[1], u1.φᵃᶠᵃ[u1.Ny + 1]), z = collect(znodes(u1, Face())),
                               halo = Oceananigans.Grids.halo_size(u1))
    return ImmersedBoundaryGrid(gr, GridFittedBottom(bh))
end
vol = PREFIX * "_volume_daily.jld2"
grid_kw = isfile(vol) ? (;) : (; grid = global_grid(vol))
TV = FieldTimeSeries(vol, "T"; backend = OnDisk(), grid_kw...); SV = FieldTimeSeries(vol, "S"; backend = OnDisk(), grid_kw...)
grid = TV.grid; ug = grid.underlying_grid
Nx, Ny, Nz = size(ug)
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center())); zc = collect(znodes(ug, Center()))
B = grid.immersed_boundary.bottom_height
bh = [B[i, j, 1] for i in 1:Nx, j in 1:Ny]
wet = bh .< 0
wet3 = [zc[k] > bh[i, j] for i in 1:Nx, j in 1:Ny, k in 1:Nz]
tdays = TV.times ./ 86400
days = parse.(Float64, split(get(ENV, "MAB_DAYS", string(round(Int, tdays[end]))), ","))

# ---- GLORYS (as in transects_vs_glorys.jl) ----
gfile(var, date) = first(filter(x -> occursin("$(var)_GLORYSDaily_$(Dates.format(date, "yyyy-mm-dd"))T", x) &&
                                     endswith(x, "-76.0_-64.0_34.0_42.0.nc"), readdir(GDIR; join = true)))
function read_glorys(var, date)
    ds = NCD.Dataset(gfile(var, date))
    lon = Float64.(ds["longitude"][:]); lat = Float64.(ds["latitude"][:]); dep = Float64.(ds["depth"][:])
    A = [ismissing(x) ? NaN : Float64(x) for x in ds[var][:, :, :, 1]]
    close(ds)
    return lon, lat, dep, A
end
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
glorys_profile(A, lon, lat, dep, i, j) =
    [dep[k] < -bh[i, j] ? bilinear(view(A, :, :, k), lon, lat, λ[i], φ[j]) : NaN for k in eachindex(dep)]

# MLD: σ₀ exceeds its 10 m value by 0.03 kg/m³ (de Boyer Montégut et al. 2004); depth positive down, top first
function mixed_layer_depth(depth, T, S)
    ok = isfinite.(T) .& isfinite.(S)
    d = depth[ok]; σ = σ₀.(T[ok], S[ok])
    (length(d) < 3 || d[1] > 10 || d[end] < 10) && return NaN
    m = searchsortedlast(d, 10.0)
    σref = m == length(d) ? σ[m] : σ[m] + (σ[m+1] - σ[m]) * (10 - d[m]) / (d[m+1] - d[m])
    for k in m+1:length(d)
        σ[k] > σref + 0.03 && return max(10.0, d[k-1] + (d[k] - d[k-1]) * (σref + 0.03 - σ[k-1]) / (σ[k] - σ[k-1]))
    end
    return d[end]
end

region(i, j) = -bh[i, j] < 200 ? 1 : -bh[i, j] < 1000 ? 2 : 3
region_names = ("shelf (<200 m)", "slope (200-1000 m)", "deep (>1000 m)")
Δλ = deg2rad(λ[2] - λ[1]); Δφ = deg2rad(φ[2] - φ[1])
f = 2Ω .* sind.(φ)

for dwant in days
    n = argmin(abs.(tdays .- dwant)); d = tdays[n]
    T = Array(interior(TV[n])); S = Array(interior(SV[n]))
    T[.!wet3] .= NaN; S[.!wet3] .= NaN
    lon, lat, dep, GT = glorys_at("thetao", d); _, _, _, GS = glorys_at("so", d)
    depth_m = reverse(-zc)

    H = fill(NaN, Nx, Ny); Hg = fill(NaN, Nx, Ny); bml = fill(NaN, Nx, Ny)
    for i in 1:Nx, j in 1:Ny
        wet[i, j] || continue
        H[i, j]  = mixed_layer_depth(depth_m, reverse(T[i, j, :]), reverse(S[i, j, :]))
        Hg[i, j] = mixed_layer_depth(dep, glorys_profile(GT, lon, lat, dep, i, j), glorys_profile(GS, lon, lat, dep, i, j))
        isfinite(H[i, j]) || continue
        ks = [k for k in 1:Nz if wet3[i, j, k] && -zc[k] <= H[i, j]]
        isempty(ks) && (ks = [Nz])
        bml[i, j] = mean(-g / ρ₀ * σ₀(T[i, j, k], S[i, j, k]) for k in ks)   # buoyancy averaged over the mixed layer
    end
    # |∇b| of the mixed-layer buoyancy, centred differences where both neighbours are wet
    gradb = fill(NaN, Nx, Ny)
    for i in 2:Nx-1, j in 2:Ny-1
        all(isfinite, (bml[i-1, j], bml[i+1, j], bml[i, j-1], bml[i, j+1])) || continue
        bx = (bml[i+1, j] - bml[i-1, j]) / (2Rₑ * cosd(φ[j]) * Δλ)
        by = (bml[i, j+1] - bml[i, j-1]) / (2Rₑ * Δφ)
        gradb[i, j] = hypot(bx, by)
    end
    # Fox-Kemper et al. (2008) restratifying buoyancy flux scale, Cₑ = 0.06 (m² s⁻³)
    wb = [0.06 * H[i, j]^2 * gradb[i, j]^2 / abs(f[j]) for i in 1:Nx, j in 1:Ny]
    bias = H .- Hg

    @printf("\n==== %s, day %.0f (%s) ====\n", basename(PREFIX), d, Dates.format(start_date + Second(round(Int, d * 86400)), "yyyy-mm-dd"))
    println("MLD bias (model − GLORYS, m) binned by the model's mixed-layer |∇b| (quartiles within each region)")
    println("region               quartile   median |∇b| (s⁻²)   mean bias   median bias   cells")
    for r in 1:3
        c = [(i, j) for i in 1:Nx, j in 1:Ny if wet[i, j] && region(i, j) == r && isfinite(bias[i, j]) && isfinite(gradb[i, j])]
        isempty(c) && continue
        gb = [gradb[i, j] for (i, j) in c]; bs = [bias[i, j] for (i, j) in c]
        q = quantile(gb, [0, 0.25, 0.5, 0.75, 1])
        for k in 1:4
            sel = (gb .>= q[k]) .& (k == 4 ? gb .<= q[k+1] : gb .< q[k+1])
            @printf("%-20s   Q%d        %9.2e          %+7.1f      %+7.1f     %5d\n", region_names[r], k,
                    median(gb[sel]), mean(bs[sel]), median(bs[sel]), count(sel))
        end
        lgb = log10.(gb)
        @printf("%-20s   corr(bias, log|∇b|) = %+.2f   corr(bias, log wb_FK) = %+.2f\n", "", cor(bs, lgb),
                cor(bs, log10.([wb[i, j] for (i, j) in c])))
    end

    # Checks that position errors (e.g. a displaced Gulf Stream) cannot produce:
    # (1) deep water south of 36°N, away from the stream, where model and GLORYS water masses agree;
    south = [(i, j) for i in 1:Nx, j in 1:Ny if wet[i, j] && region(i, j) == 3 && φ[j] < 36 && isfinite(bias[i, j])]
    bs = [bias[i, j] for (i, j) in south]
    @printf("deep, south of 36°N:  model MLD median %.1f m, GLORYS %.1f m, bias mean %+.1f median %+.1f (%d cells)\n",
            median(H[i, j] for (i, j) in south), median(Hg[i, j] for (i, j) in south), mean(bs), median(bs), length(south))
    # (2) water-mass matched: median MLD in classes of surface temperature, wherever those columns are in each field
    sst_m = T[:, :, Nz]
    sst_g = [wet[i, j] ? bilinear(view(GT, :, :, 1), lon, lat, λ[i], φ[j]) : NaN for i in 1:Nx, j in 1:Ny]
    println("SST class (°C)    model MLD (median, cells)    GLORYS MLD (median, cells)    difference")
    for (lo, hi) in ((4, 8), (8, 12), (12, 16), (16, 20), (20, 24), (24, 30))
        m = [H[i, j] for i in 1:Nx, j in 1:Ny if wet[i, j] && lo <= sst_m[i, j] < hi && isfinite(H[i, j])]
        gl = [Hg[i, j] for i in 1:Nx, j in 1:Ny if wet[i, j] && lo <= sst_g[i, j] < hi && isfinite(Hg[i, j])]
        (length(m) < 20 || length(gl) < 20) && continue
        @printf("%5.0f-%-5.0f         %6.1f m  %6d              %6.1f m  %6d            %+6.1f m\n",
                lo, hi, median(m), length(m), median(gl), length(gl), median(m) - median(gl))
    end

    fig = Figure(size = (1700, 900), fontsize = 14)
    Label(fig[0, 1:4], @sprintf("%s day %.0f: is the MLD bias at fronts?", basename(PREFIX), d), fontsize = 18)
    for (col, (A, t, cm, cr)) in enumerate(((bias, "MLD bias, model − GLORYS (m)", :balance, (-40, 40)),
                                           (log10.(gradb), "log₁₀ mixed-layer |∇b| (s⁻²)", :viridis, (-9, -6.5)),
                                           (log10.(wb), "log₁₀ Fox-Kemper wb scale (m² s⁻³)", :viridis, (-11, -7))))
        ax = Axis(fig[1, col]; title = t, aspect = DataAspect())
        hm = heatmap!(ax, λ, φ, A; colormap = cm, colorrange = cr, nan_color = :gray85)
        Colorbar(fig[2, col], hm; vertical = false, flipaxis = false)
    end
    ax = Axis(fig[1, 4]; title = "MLD bias vs mixed-layer |∇b|", xlabel = "log₁₀ |∇b| (s⁻²)", ylabel = "bias (m)")
    for (r, color) in zip(1:3, (:seagreen, :orange, :royalblue))
        c = [(i, j) for i in 1:Nx, j in 1:Ny if wet[i, j] && region(i, j) == r && isfinite(bias[i, j]) && isfinite(gradb[i, j])]
        scatter!(ax, [log10(gradb[i, j]) for (i, j) in c], [bias[i, j] for (i, j) in c]; markersize = 3, color = (color, 0.4),
                 label = region_names[r])
    end
    hlines!(ax, 0; color = :black); ylims!(ax, -60, 100); axislegend(ax; position = :lt)
    out = PREFIX * @sprintf("_mld_bias_fronts_day%02d.png", round(Int, d))
    save(out, fig; px_per_unit = 1.2); println("saved ", out)
end
