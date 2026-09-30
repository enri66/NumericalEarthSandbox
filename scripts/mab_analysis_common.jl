# Shared pieces of the MAB analysis scripts (transects_vs_glorys.jl, seasonal_sections.jl, animate_surface_year.jl):
# opening a run's output (including script 05's per-rank files), GLORYS daily fields, and the mixed-layer and
# seasonal-thermocline depths of a profile.
using Oceananigans, NumericalEarth, Printf, Statistics, Dates
using Oceananigans.Grids: znodes, λnodes, φnodes
const NCD = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "NCDatasets")]           # loaded by NumericalEarth
const SWP = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "SeawaterPolynomials")]  # loaded by Oceananigans
const JLD2 = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "JLD2")]

const OUT  = joinpath(homedir(), "Data", "mab_glorys_obc")
const GDIR = get(ENV, "MAB_DATA_DIR", joinpath(homedir(), "Data", "NumericalEarth"))
const start_date = DateTime(2019, 4, 1)

const eos = SWP.TEOS10.TEOS10EquationOfState()
σ₀(T, S) = SWP.ρ(T, S, 0, eos) - 1000

# A run tag in ~/Data/mab_glorys_obc or, as the run scripts take it, an absolute path prefix.
run_prefix(tag) = isabspath(tag) ? tag : joinpath(OUT, tag)

# ---------------- run output ----------------
# Script 05 writes one file per MPI rank (`_rank<r>`), and Oceananigans combines them on an ImmersedBoundaryGrid only
# when given the global grid, so that grid is rebuilt from the ranks' own grids (slabs in x). Reading those grids needs
# the MPI environment the run used (e.g. /t0/workdir/enrique/mpi05_ib on triton).
function global_grid(path)
    stem = basename(splitext(path)[1])
    paths = filter(f -> occursin(Regex("^" * stem * "_rank\\d+\\.jld2\$"), basename(f)), readdir(dirname(path); join = true))
    grids = [JLD2.jldopen(f -> f["serialized/grid"], p) for p in paths]
    ugs = [g.underlying_grid for g in grids]
    order = sortperm([ug.λᶠᵃᵃ[1] for ug in ugs])
    function bottom(r)
        b = grids[r].immersed_boundary.bottom_height
        n, m = ugs[r].Nx, ugs[r].Ny
        return ndims(b) == 2 ? collect(b[1:n, 1:m]) : collect(b[1:n, 1:m, 1])
    end
    bh = cat([bottom(r) for r in order]...; dims = 1)
    first_ug, last_ug = ugs[order[1]], ugs[order[end]]
    Ny, Nz = first_ug.Ny, first_ug.Nz
    g = LatitudeLongitudeGrid(CPU(); size = (size(bh, 1), Ny, Nz),
                              longitude = (first_ug.λᶠᵃᵃ[1], last_ug.λᶠᵃᵃ[last_ug.Nx + 1]),
                              latitude = (first_ug.φᵃᶠᵃ[1], first_ug.φᵃᶠᵃ[Ny + 1]),
                              z = collect(znodes(first_ug, Face())), halo = Oceananigans.Grids.halo_size(first_ug))
    return ImmersedBoundaryGrid(g, GridFittedBottom(bh))
end

# A FieldTimeSeries from a run's output file, whether it was written as one file or one per rank
function open_series(path, name; kw...)
    grid_kw = isfile(path) ? (;) : (; grid = global_grid(path))
    return FieldTimeSeries(path, name; grid_kw..., kw...)
end

# All frames of a horizontal-slice output (the surface and SSH files) as an (x, y, frame) array, with their times.
# Oceananigans' combining reader expects full-depth fields, so script 05's per-rank slices are joined here instead:
# each rank's interior is cut out of its array, with the halo widths inferred from the array size (the split-explicit
# free surface keeps wider halos than the other fields).
isface(L) = L === Face || L isa Face

function surface_series(path, name)
    if isfile(path)
        fts = FieldTimeSeries(path, name)
        return cat([Array(interior(fts[n]))[:, :, 1] for n in eachindex(fts.times)]...; dims = 3), collect(fts.times)
    end
    stem = basename(splitext(path)[1])
    paths = filter(f -> occursin(Regex("^" * stem * "_rank\\d+\\.jld2\$"), basename(f)), readdir(dirname(path); join = true))
    files = [JLD2.jldopen(p, "r") for p in paths]
    ugs = [f["serialized/grid"].underlying_grid for f in files]
    order = sortperm([ug.λᶠᵃᵃ[1] for ug in ugs])
    files, ugs = files[order], ugs[order]
    loc = files[1]["timeseries/$name/serialized/location"]
    iterations = sort(parse.(Int, filter(!=("serialized"), collect(keys(files[1]["timeseries/$name"])))))
    times = [files[1]["timeseries/t/$it"] for it in iterations]
    frames = map(iterations) do it
        parts = map(eachindex(files)) do r
            a = files[r]["timeseries/$name/$it"]
            a = ndims(a) == 2 ? a : a[:, :, end]
            nx = ugs[r].Nx + (isface(loc[1]) && r == length(files) ? 1 : 0)
            ny = ugs[r].Ny + (isface(loc[2]) ? 1 : 0)
            hx = (size(a, 1) - nx) ÷ 2; hy = (size(a, 2) - ny) ÷ 2
            a[hx+1:hx+nx, hy+1:hy+ny]
        end
        cat(parts...; dims = 1)
    end
    close.(files)
    return cat(frames...; dims = 3), times
end

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

# GLORYS profile (native levels) at (x, y); levels deeper than H (the model's depth there) are dropped
glorys_profile(A, lon, lat, dep, x, y, H) =
    [dep[k] < H ? bilinear(view(A, :, :, k), lon, lat, x, y) : NaN for k in eachindex(dep)]

# ---------------- mixed layer and thermocline from one profile ----------------
# Mixed-layer depth: where σ₀ first exceeds its 10 m value by 0.03 kg/m³ (de Boyer Montégut et al. 2004).
# Seasonal thermocline depth: the midpoint of the largest temperature decrease between two adjacent levels, between the
# mixed-layer base and 300 m, where that decrease is at least 0.02 °C/m (otherwise there is none).
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
