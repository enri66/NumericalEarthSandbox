# Response of the phytoplankton to the light supply: compares runs of bgc/mab_bgc.jl that differ only in MAB_BGC_PAR_FRACTION or
# MAB_BGC_WATER_BLUE_ATTENUATION (restarted from the same checkpoint). For the frames the runs have in common it prints, for the
# shelf, slope and deep regions and the CNSM and ISSM cells, the mean surface chlorophyll (P × 1.31 mg Chl/mmol N, the top level of
# the de-tided daily output), the phytoplankton biomass in the top 50 m (mmol N/m²) and, from the limitation files, the surface
# growth rate and PAR, with the change relative to the first run (the reference).
#
#   BGC_TAGS="reference=/t0/.../tagA,more light=/t0/.../tagB,..." julia --project=<env> bgc/light_sensitivity.jl
using Printf, Statistics, Dates
using Oceananigans
using Oceananigans.Grids: znodes
using NumericalEarth.DataWrangling: jldopen

const RUNS = [(split(kv, "=")[1], split(kv, "=")[2]) for kv in split(ENV["BGC_TAGS"], ",")]
const H, Nz, CPD = 7, parse(Int, get(ENV, "MAB_NZ", "50")), parse(Int, get(ENV, "MAB_CELLS_PER_DEGREE", "12"))
const λ_west, φ_south = -76 + 2 / 12, 34 + 2 / 12
const CHL_PER_N = 1.31
const MOORINGS = ("CNSM" => (-70.88, 40.10, 147.0), "ISSM" => (-70.88, 40.36, 92.0))

function exponential_scale(N, H, Δtop)
    top(h) = H * expm1(H / N / h) / expm1(H / h)
    lo, hi = H / 50, 100H
    for _ in 1:200; mid = (lo + hi) / 2; top(mid) > Δtop ? (hi = mid) : (lo = mid); end
    return (lo + hi) / 2
end
Δz_top, Z_BOTTOM = parse(Float64, get(ENV, "MAB_DZ_TOP", "2")), parse(Float64, get(ENV, "MAB_ZBOTTOM", "5500"))
zgrid = RectilinearGrid(size = (1, 1, Nz), x = (0, 1), y = (0, 1),
                        z = ExponentialDiscretization(Nz, -Z_BOTTOM, 0; scale = exponential_scale(Nz, Z_BOTTOM, Δz_top)))
const zf = collect(znodes(zgrid, Face())); const Δz = diff(zf)
const K50 = findall(k -> -zf[k+1] < 50, 1:Nz)             # levels (bottom first) whose upper face is shallower than 50 m

iters(f, name) = sort(parse.(Int, filter(k -> k != "serialized", keys(f["timeseries/$name"]))))
function layout(tag)
    files = sort(filter(f -> occursin(Regex("^" * basename(tag) * "_bottom_rank\\d+\\.jld2\$"), basename(f)), readdir(dirname(tag); join = true)))
    parts = map(files) do f
        jldopen(f) do file; (I = file["I_OFF"], J = file["J_OFF"], kb = file["k_bottom"], dzb = file["Δz_bottom"]); end
    end
    return parts
end
parts = layout(first(RUNS)[2]); NR = length(parts)
const Nx = maximum(p.I + size(p.kb, 1) for p in parts); const Ny = maximum(p.J + size(p.kb, 2) for p in parts)
kb = zeros(Int, Nx, Ny); depth = zeros(Nx, Ny)
for p in parts, j in axes(p.kb, 2), i in axes(p.kb, 1)
    k = p.kb[i, j]; kb[p.I + i, p.J + j] = k
    depth[p.I + i, p.J + j] = k > 0 ? sum(Δz[k+1:end]) + p.dzb[i, j] : 0.0
end
wet = kb .> 0
λm(i) = λ_west + (i - 0.5) / CPD; φm(j) = φ_south + (j - 0.5) / CPD
regions = (shelf = (depth .> 0) .& (depth .< 100), slope = (depth .≥ 100) .& (depth .< 1000), deep = depth .≥ 1000)
cells = Dict{String, Tuple{Int, Int}}()
for (name, (lon, lat, nominal)) in MOORINGS
    cand = [(i, j) for i in 1:Nx, j in 1:Ny if wet[i, j] && abs(λm(i) - lon) * cosd(lat) < 0.13 && abs(φm(j) - lat) < 0.13]
    cells[name] = cand[argmin([abs(depth[i, j] - nominal) / nominal + 0.5 * hypot((λm(i) - lon) * cosd(lat), φm(j) - lat) / 0.13 for (i, j) in cand])]
end

# per run: the frame times of the de-tided volume output and the limitation output
function read_run(tag)
    ranks = 0:NR-1
    vfiles = ["$(tag)_bgc_volume_daily_rank$(r).jld2" for r in ranks]
    lfiles = ["$(tag)_bgc_limitation_rank$(r).jld2" for r in ranks]
    vt = jldopen(f -> [f["timeseries/t/$it"] for it in iters(f, "P")], vfiles[1])
    lt = jldopen(f -> [f["timeseries/t/$it"] for it in iters(f, "PAR")], lfiles[1])
    return (; vfiles, lfiles, vt, lt)
end
runs = [(label, tag, read_run(tag)) for (label, tag) in RUNS]
common_v = sort(collect(intersect((Set(round.(Int, r.vt)) for (_, _, r) in runs)...)))
common_l = sort(collect(intersect((Set(round.(Int, r.lt)) for (_, _, r) in runs)...)))
println("frames in common: $(length(common_v)) de-tided volume frames, $(length(common_l)) limitation frames")
isempty(common_v) && error("no common volume frames")

function volume_P(run, t)
    n = argmin(abs.(run.vt .- t))
    A = zeros(Float32, Nx, Ny, Nz)
    for (p, file) in zip(parts, run.vfiles)
        jldopen(file) do f
            it = iters(f, "P")[n]
            B = f["timeseries/P/$it"][H+1:end-H, H+1:end-H, H+1:end-H]
            A[p.I+1:p.I+size(B, 1), p.J+1:p.J+size(B, 2), :] .= B
        end
    end
    return A
end
function limitation_surface(run, t, name)
    n = argmin(abs.(run.lt .- t))
    A = fill(NaN32, Nx, Ny)
    for (p, file) in zip(parts, run.lfiles)
        jldopen(file) do f
            it = iters(f, name)[n]
            B = f["timeseries/$name/$it"]
            size(B, 1) == size(p.kb, 1) + 2H && (B = B[H+1:end-H, H+1:end-H, :])
            A[p.I+1:p.I+size(B, 1), p.J+1:p.J+size(B, 2)] .= B[:, :, end]
        end
    end
    return A
end

# quantity series: label => region/mooring => mean over the common frames
function stats(run)
    chl = Dict(k => Float64[] for k in (keys(regions)..., Symbol.(keys(cells))...))
    bio = Dict(k => Float64[] for k in keys(chl)); growth = Dict(k => Float64[] for k in keys(chl)); par = Dict(k => Float64[] for k in keys(chl))
    function add!(d, k, A, mask) ; v = A[mask .& isfinite.(A)]; isempty(v) || push!(d[k], mean(v)); end
    for t in common_v
        P = volume_P(run, t)
        Pchl = P[:, :, Nz] .* CHL_PER_N
        Pint = sum(P[:, :, K50[i]] .* Δz[K50[i]] for i in eachindex(K50))
        G = t in common_l ? limitation_surface(run, t, "growth") : nothing
        R = t in common_l ? limitation_surface(run, t, "PAR") : nothing
        for k in keys(chl)
            mask = k in keys(regions) ? regions[k] : (m = falses(Nx, Ny); m[cells[String(k)]...] = true; m)
            add!(chl, k, Pchl, mask); add!(bio, k, Pint, mask)
            G === nothing || (add!(growth, k, G, mask); add!(par, k, R, mask))
        end
    end
    return (; chl, bio, growth, par)
end
results = [(label, stats(run)) for (label, _, run) in runs]
ref = first(results)[2]
for (title, field, fmt) in (("surface chlorophyll (mg/m³)", :chl, "%8.3f"), ("phytoplankton in the top 50 m (mmol N/m²)", :bio, "%8.2f"),
                             ("surface growth rate (1/day)", :growth, "%8.3f"), ("surface PAR (W/m²)", :par, "%8.1f"))
    println("\n", title, ": mean over the common frames (change relative to ", first(results)[1], ")")
    println(@sprintf("%-14s", ""), join([@sprintf("%22s", String(k)) for k in keys(ref.chl)]))
    for (label, r) in results
        row = @sprintf("%-14s", label)
        for k in keys(ref.chl)
            v = isempty(getfield(r, field)[k]) ? NaN : mean(getfield(r, field)[k]); v0 = isempty(getfield(ref, field)[k]) ? NaN : mean(getfield(ref, field)[k])
            row *= @sprintf("%11.3f (%+6.1f%%)", v, 100 * (v / v0 - 1))
        end
        println(row)
    end
end
