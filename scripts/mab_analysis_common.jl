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
const start_date = DateTime(get(ENV, "MAB_START_DATE", "2019-04-01"))   # the run's start (its clock's zero)

const eos = SWP.TEOS10.TEOS10EquationOfState()
σ₀(T, S) = SWP.ρ(T, S, 0, eos) - 1000

# A run tag in ~/Data/mab_glorys_obc or, as the run scripts take it, an absolute path prefix.
run_prefix(tag) = isabspath(tag) ? tag : joinpath(OUT, tag)

# ---------------- run output ----------------
# Script 05 writes one file per MPI rank (`_rank<r>`), in a px × py layout (MAB_PARTITION_X, MAB_PARTITION_Y). The ranks'
# own grids tell where each piece sits: `rank_blocks` returns the file paths as a matrix [x-block, y-block], with the
# blocks ordered by the longitude and latitude of each piece's first face.
function rank_blocks(path)
    stem = basename(splitext(path)[1])
    paths = filter(f -> occursin(Regex("^" * stem * "_rank\\d+\\.jld2\$"), basename(f)), readdir(dirname(path); join = true))
    ugs = [JLD2.jldopen(f -> f["serialized/grid"], p).underlying_grid for p in paths]
    λ0 = [round(ug.λᶠᵃᵃ[1]; digits = 6) for ug in ugs]; φ0 = [round(ug.φᵃᶠᵃ[1]; digits = 6) for ug in ugs]
    λs = sort(unique(λ0)); φs = sort(unique(φ0))
    blocks = Matrix{String}(undef, length(λs), length(φs))
    for (n, p) in enumerate(paths)
        blocks[findfirst(==(λ0[n]), λs), findfirst(==(φ0[n]), φs)] = p
    end
    return blocks
end

# Global grid, rebuilt from the ranks' own grids (reading them needs the MPI environment the run used, e.g.
# /t0/workdir/enrique/mpi05_ib on triton)
function global_grid(path)
    blocks = rank_blocks(path)
    grids = Dict(b => JLD2.jldopen(f -> f["serialized/grid"], blocks[b]) for b in CartesianIndices(blocks))
    ug(b) = grids[b].underlying_grid
    function bottom(b)
        a = grids[b].immersed_boundary.bottom_height
        n, m = ug(b).Nx, ug(b).Ny
        return ndims(a) == 2 ? collect(a[1:n, 1:m]) : collect(a[1:n, 1:m, 1])
    end
    bh = reduce(vcat, [reduce(hcat, [bottom(CartesianIndex(a, b)) for b in 1:size(blocks, 2)]) for a in 1:size(blocks, 1)])
    first_ug, last_ug = ug(CartesianIndex(1, 1)), ug(CartesianIndex(size(blocks)...))
    Nz = first_ug.Nz
    g = LatitudeLongitudeGrid(CPU(); size = (size(bh, 1), size(bh, 2), Nz),
                              longitude = (first_ug.λᶠᵃᵃ[1], last_ug.λᶠᵃᵃ[last_ug.Nx + 1]),
                              latitude = (first_ug.φᵃᶠᵃ[1], last_ug.φᵃᶠᵃ[last_ug.Ny + 1]),
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
# free surface keeps wider halos than the other fields), and the pieces are laid out in the ranks' px × py arrangement.
isface(L) = L === Face || L isa Face

function surface_series(path, name)
    if isfile(path)
        fts = FieldTimeSeries(path, name)
        return cat([Array(interior(fts[n]))[:, :, 1] for n in eachindex(fts.times)]...; dims = 3), collect(fts.times)
    end
    blocks = rank_blocks(path)
    nbx, nby = size(blocks)
    files = Dict(b => JLD2.jldopen(blocks[b], "r") for b in CartesianIndices(blocks))
    ref = files[CartesianIndex(1, 1)]
    loc = ref["timeseries/$name/serialized/location"]
    iterations = sort(parse.(Int, filter(!=("serialized"), collect(keys(ref["timeseries/$name"])))))
    times = [ref["timeseries/t/$it"] for it in iterations]
    frames = map(iterations) do it
        rows = map(1:nbx) do a
            cols = map(1:nby) do b
                f = files[CartesianIndex(a, b)]; u = f["serialized/grid"].underlying_grid
                arr = f["timeseries/$name/$it"]
                arr = ndims(arr) == 2 ? arr : arr[:, :, end]
                nx = u.Nx + (isface(loc[1]) && a == nbx ? 1 : 0)
                ny = u.Ny + (isface(loc[2]) && b == nby ? 1 : 0)
                hx = (size(arr, 1) - nx) ÷ 2; hy = (size(arr, 2) - ny) ÷ 2
                arr[hx+1:hx+nx, hy+1:hy+ny]
            end
            reduce(hcat, cols)
        end
        reduce(vcat, rows)
    end
    close.(values(files))
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
