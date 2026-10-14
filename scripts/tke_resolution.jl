# Near-surface TKE at several horizontal resolutions, and the size of the submesoscale restratification sink proposed
# for CATKE (G. Wagner, "Unresolved restratification and the CATKE energy budget", 10 Oct 2026),
#   𝒥_b = −C_Δ Δ |∇ₕb| √e     (his eq. 2.10; Δ the grid spacing, C_Δ a dimensionless coefficient).
#
# Reads the checkpoints of script-05 runs (instantaneous T, S, u, v, e, κc, κu, written every MAB_CHECKPOINT_EVERY days)
# at the days found in all of them. For every wet column it takes the mixed-layer depth h of mab_analysis_common.jl and
# averages over two layers, the whole mixed layer (0 to h) and its lower half (h/2 to h):
#   e                       CATKE's TKE
#   G = Δ |∇ₕb| √e          the proposed sink per unit C_Δ (b = −g σ₀/ρ₀, ∇ₕb along model levels)
#   κc N²                   CATKE's own buoyancy sink (stable stratification)
#   P = κu S²               CATKE's shear production
#   ℓc = κc / √e            CATKE's tracer mixing length; the covariance bound (his eq. 2.4) allows at most
#                           C_Δ ≤ ℓc / (√12 h_χ), for a vertical scale h_χ of the unresolved buoyancy variance
# Columns mixed to within 10 % of the bottom are left out (their lower mixed layer is the bottom boundary layer), as is
# a band of MAB_EDGE_DEG degrees along the domain edges (sponges). Each run's column values are averaged onto the
# coarsest run's grid (blocks of cells), so all resolutions are compared cell by cell on one grid.
#
# Prints, by region and for all days together: the median TKE in both layers at each resolution; for the lower half,
# G against κc N² and P, ℓc and the C_Δ bound; and the coarse runs' TKE excess over the finest one against the coarse
# run's G (is the excess where the proposed sink would act?). Saves a figure and the block values. Usage:
#   MAB_TAGS=/t0/.../res_test/tke12,/t0/.../res_test/tke24,/t0/.../res_test/hr100_fx3 MAB_START_DATE=2019-08-29 \
#       julia -t 16 --project=/t0/workdir/enrique/mpi05_ib scripts/tke_resolution.jl
# MAB_DAYS_LIST: days to use (default: every checkpoint day after day 0 present in all runs)
# MAB_EDGE_DEG:  width of the edge band left out (default 0.75)
# MAB_OUT:       output prefix (default: the first tag with _tke_resolution)
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const TAGS = [run_prefix(t) for t in split(ENV["MAB_TAGS"], ",")]
const EDGE = parse(Float64, get(ENV, "MAB_EDGE_DEG", "0.75"))
const OUTP = get(ENV, "MAB_OUT", TAGS[1] * "_tke_resolution")
const g, ρ₀, R = 9.81, 1026.0, 6.371e6

# ---------------- reading checkpoints ----------------
rank_number(path) = parse(Int, match(r"_rank(\d+)\.jld2$", path).captures[1])
checkpoint_path(prefix, r, it) = "$(prefix)_checkpoint_rank$(r)_iteration$(it).jld2"
const MODEL = "simulation/model/ocean/model"

function checkpoint_days(prefix)
    files = filter(f -> occursin(Regex("^" * basename(prefix) * "_checkpoint_rank0_iteration\\d+\\.jld2\$"), f), readdir(dirname(prefix)))
    its = sort(parse.(Int, [match(r"iteration(\d+)\.jld2$", f).captures[1] for f in files]))
    days = map(its) do it
        JLD2.jldopen(f -> f["$MODEL/clock"].time, checkpoint_path(prefix, 0, it)) / 86400
    end
    return Dict(round(Int, d) => it for (d, it) in zip(days, its) if abs(d - round(d)) < 1e-6)
end

# Interior cell-centre arrays of T, S, e, u and v (averaged to centres) and face arrays of κc and κu, laid out globally
function read_state(prefix, it, blocks, ugs)
    nbx, nby = size(blocks)
    nx = [ugs[a, 1].Nx for a in 1:nbx]; ny = [ugs[1, b].Ny for b in 1:nby]; Nz = ugs[1, 1].Nz
    X0 = cumsum([0; nx[1:end-1]]); Y0 = cumsum([0; ny[1:end-1]])
    Nx, Ny = sum(nx), sum(ny)
    out = Dict(n => fill(NaN, Nx, Ny, Nz) for n in (:T, :S, :e, :u, :v))
    merge!(out, Dict(n => fill(NaN, Nx, Ny, Nz + 1) for n in (:κc, :κu)))
    for a in 1:nbx, b in 1:nby
        f = JLD2.jldopen(checkpoint_path(prefix, rank_number(blocks[a, b]), it), "r")
        n, m = nx[a], ny[b]
        I = X0[a] .+ (1:n); J = Y0[b] .+ (1:m)
        for (name, key) in ((:T, "tracers/T"), (:S, "tracers/S"), (:e, "tracers/e"), (:u, "velocities/u"), (:v, "velocities/v"))
            d = f["$MODEL/$key/data"]
            hx = (size(d, 1) - n) ÷ 2; hy = (size(d, 2) - m) ÷ 2; hz = (size(d, 3) - Nz) ÷ 2
            K = hz .+ (1:Nz)
            if name === :u      # x faces: centre i is between faces i and i + 1 (the last one is in the halo or on the boundary)
                out[name][I, J, :] .= (d[hx .+ (1:n), hy .+ (1:m), K] .+ d[hx .+ (2:n+1), hy .+ (1:m), K]) ./ 2
            elseif name === :v
                out[name][I, J, :] .= (d[hx .+ (1:n), hy .+ (1:m), K] .+ d[hx .+ (1:n), hy .+ (2:m+1), K]) ./ 2
            else
                out[name][I, J, :] .= d[hx .+ (1:n), hy .+ (1:m), K]
            end
        end
        for (name, key) in ((:κc, "closure_fields/κc"), (:κu, "closure_fields/κu"))
            d = f["$MODEL/$key/data"]
            hx = (size(d, 1) - n) ÷ 2; hy = (size(d, 2) - m) ÷ 2; hz = (size(d, 3) - Nz - 1) ÷ 2
            out[name][I, J, :] .= d[hx .+ (1:n), hy .+ (1:m), hz .+ (1:Nz+1)]
        end
        close(f)
    end
    return out
end

# ---------------- column diagnostics ----------------
# Mean of x over the levels whose depth lies in [lo, hi] (weights w), or the level nearest the middle if none does
function layer_mean(x, depth, w, lo, hi)
    s = 0.0; ws = 0.0
    for k in eachindex(x)
        (lo <= depth[k] <= hi && isfinite(x[k])) || continue
        s += w[k] * x[k]; ws += w[k]
    end
    ws > 0 && return s / ws
    k = argmin(abs.(depth .- (lo + hi) / 2))
    return x[k]
end

const QUANTITIES = (:h, :e_ml, :e_low, :G_low, :B_low, :P_low, :ℓ_low, :∇b_low)

function column_diagnostics(prefix, it, blocks, ugs, grid)
    st = read_state(prefix, it, blocks, ugs)
    ug = grid.underlying_grid
    Nx, Ny, Nz = size(ug)
    zc = collect(znodes(ug, Center())); zf = collect(znodes(ug, Face()))
    φ = collect(φnodes(ug, Center())); λ = collect(λnodes(ug, Center()))
    Δλ = deg2rad(λ[2] - λ[1]); Δφ = deg2rad(φ[2] - φ[1])
    dx = [R * cosd(φ[j]) * Δλ for j in 1:Ny]; dy = R * Δφ
    bh = [grid.immersed_boundary.bottom_height[i, j, 1] for i in 1:Nx, j in 1:Ny]
    for k in 1:Nz, j in 1:Ny, i in 1:Nx            # below the bottom: no data
        zc[k] < bh[i, j] && (st[:T][i, j, k] = NaN; st[:S][i, j, k] = NaN; st[:e][i, j, k] = NaN)
    end
    b = similar(st[:T])
    Threads.@threads for k in 1:Nz
        for j in 1:Ny, i in 1:Nx
            t, s = st[:T][i, j, k], st[:S][i, j, k]
            b[i, j, k] = isfinite(t) && isfinite(s) ? -g * σ₀(t, s) / ρ₀ : NaN
        end
    end
    # |∇ₕb| at cell centres along levels: centred differences, one-sided next to land
    ∇b = fill(NaN, Nx, Ny, Nz)
    Threads.@threads for k in 1:Nz
        for j in 1:Ny, i in 1:Nx
            isfinite(b[i, j, k]) || continue
            ok(p, q) = 1 <= p <= Nx && 1 <= q <= Ny && isfinite(b[p, q, k])
            gx = ok(i - 1, j) && ok(i + 1, j) ? (b[i+1, j, k] - b[i-1, j, k]) / 2dx[j] :
                 ok(i + 1, j) ? (b[i+1, j, k] - b[i, j, k]) / dx[j] : ok(i - 1, j) ? (b[i, j, k] - b[i-1, j, k]) / dx[j] : NaN
            gy = ok(i, j - 1) && ok(i, j + 1) ? (b[i, j+1, k] - b[i, j-1, k]) / 2dy :
                 ok(i, j + 1) ? (b[i, j+1, k] - b[i, j, k]) / dy : ok(i, j - 1) ? (b[i, j, k] - b[i, j-1, k]) / dy : NaN
            ∇b[i, j, k] = hypot(gx, gy)
        end
    end
    depth_c = reverse(-zc)                          # top to bottom
    Δzc = reverse(diff(zf))
    out = Dict(q => fill(NaN, Nx, Ny) for q in QUANTITIES)
    Threads.@threads for j in 1:Ny
        Δ = sqrt(dx[j] * dy)
        for i in 1:Nx
            bh[i, j] < 0 || continue
            (λ[i] - λ[1] > EDGE && λ[end] - λ[i] > EDGE && φ[j] - φ[1] > EDGE && φ[end] - φ[j] > EDGE) || continue
            T = reverse(st[:T][i, j, :]); S = reverse(st[:S][i, j, :])
            h, _ = column_depths(depth_c, T, S)
            (isfinite(h) && h < 0.9 * -bh[i, j]) || continue
            e = reverse(st[:e][i, j, :])
            G = [Δ * ∇b[i, j, k] * sqrt(max(st[:e][i, j, k], 0)) for k in Nz:-1:1]
            ∇bc = reverse(∇b[i, j, :])
            # faces k = 2..Nz (between centres k - 1 and k), from the top down
            fk = Nz:-1:2
            depth_f = [-zf[k] for k in fk]
            Δzf = [zc[k] - zc[k-1] for k in fk]
            N² = [(b[i, j, k] - b[i, j, k-1]) / (zc[k] - zc[k-1]) for k in fk]
            S² = [((st[:u][i, j, k] - st[:u][i, j, k-1])^2 + (st[:v][i, j, k] - st[:v][i, j, k-1])^2) / (zc[k] - zc[k-1])^2 for k in fk]
            ef = [(st[:e][i, j, k] + st[:e][i, j, k-1]) / 2 for k in fk]
            κc = [st[:κc][i, j, k] for k in fk]; κu = [st[:κu][i, j, k] for k in fk]
            B = κc .* N²; P = κu .* S²; ℓ = κc ./ sqrt.(max.(ef, 1e-12))
            out[:h][i, j] = h
            out[:e_ml][i, j] = layer_mean(e, depth_c, Δzc, 0.0, h)
            out[:e_low][i, j] = layer_mean(e, depth_c, Δzc, h / 2, h)
            out[:G_low][i, j] = layer_mean(G, depth_c, Δzc, h / 2, h)
            out[:∇b_low][i, j] = layer_mean(∇bc, depth_c, Δzc, h / 2, h)
            out[:B_low][i, j] = layer_mean(B, depth_f, Δzf, h / 2, h)
            out[:P_low][i, j] = layer_mean(P, depth_f, Δzf, h / 2, h)
            out[:ℓ_low][i, j] = layer_mean(ℓ, depth_f, Δzf, h / 2, h)
        end
    end
    return out
end

# Average over blocks of r × r cells (finite values only)
function coarsen(A, r)
    Nx, Ny = size(A) .÷ r
    C = fill(NaN, Nx, Ny)
    for I in 1:Nx, J in 1:Ny
        v = filter(isfinite, A[(I-1)*r+1:I*r, (J-1)*r+1:J*r])
        isempty(v) || (C[I, J] = mean(v))
    end
    return C
end

# ---------------- all runs ----------------
runs = map(TAGS) do prefix
    blocks = rank_blocks(prefix * ".jld2")
    ugs = map(p -> JLD2.jldopen(f -> f["serialized/grid"], p).underlying_grid, blocks)
    grid = global_grid(prefix * ".jld2")
    (; prefix, name = basename(prefix), blocks, ugs, grid, days = checkpoint_days(prefix), Nx = grid.underlying_grid.Nx)
end
days = haskey(ENV, "MAB_DAYS_LIST") ? parse.(Int, split(ENV["MAB_DAYS_LIST"], ",")) :
       sort(collect(intersect([Set(k for k in keys(r.days) if k > 0) for r in runs]...)))
coarse = runs[argmin([r.Nx for r in runs])]
Ncx = coarse.Nx
factor(r) = r.Nx ÷ Ncx
all(r -> r.Nx == factor(r) * Ncx, runs) || error("grid sizes are not multiples of the coarsest: $([r.Nx for r in runs])")
println("runs: ", join(["$(r.name) ($(r.Nx) columns in x, $(factor(r))× the coarsest)" for r in runs], ", "))
println("days: ", join(days, ", "))

cb = coarse.grid.immersed_boundary.bottom_height
cbh = [cb[i, j, 1] for i in 1:Ncx, j in 1:coarse.grid.underlying_grid.Ny]
region_of(d) = d < 200 ? 1 : d < 1000 ? 2 : 3
regions = ("shelf (<200 m)", "slope (200-1000 m)", "deep (>1000 m)")
creg = [cbh[i, j] < 0 ? region_of(-cbh[i, j]) : 0 for i in axes(cbh, 1), j in axes(cbh, 2)]

# block values: blk[run][quantity] is (Ncx, Ncy, ndays)
blk = map(runs) do r
    d = Dict(q => fill(NaN, size(cbh)..., length(days)) for q in QUANTITIES)
    for (n, day) in enumerate(days)
        t0 = time()
        c = column_diagnostics(r.prefix, r.days[day], r.blocks, r.ugs, r.grid)
        for q in QUANTITIES
            d[q][:, :, n] .= coarsen(c[q], factor(r))
        end
        @printf("  %s day %d done (%.0f s)\n", r.name, day, time() - t0)
        flush(stdout)
    end
    d
end

# compare only blocks with values in every run on the same day
common = trues(size(cbh)..., length(days))
for d in blk, q in (:e_ml, :e_low)
    common .&= isfinite.(d[q])
end
med(x) = (y = filter(isfinite, x); isempty(y) ? NaN : median(y))
ranks(x) = (r = similar(x); r[sortperm(x)] = 1:length(x); r)
rankcor(x, y) = cor(ranks(x), ranks(y))       # Spearman: a few extreme blocks do not set it
sel(A, reg) = [A[i, j, n] for i in axes(A, 1), j in axes(A, 2), n in axes(A, 3) if common[i, j, n] && creg[i, j] == reg]

println("\n1. Median TKE (10⁻⁵ m²/s²) over common blocks and days, and ratio to $(runs[1].name)")
@printf("%-20s %6s", "region", "blocks")
for r in runs
    @printf("  %22s", r.name)
end
println()
for (reg, rname) in enumerate(regions)
    for (q, label) in ((:e_ml, "whole ML"), (:e_low, "lower half"))
        @printf("%-20s %6d", "$rname $label"[1:min(end, 20)], length(sel(blk[1][q], reg)))
        base = med(sel(blk[1][q], reg))
        for d in blk
            m = med(sel(d[q], reg))
            @printf("  %12.2f (%5.2f×)", 1e5 * m, m / base)
        end
        println()
    end
end
println("median mixed-layer depth (m): ", join(["$(r.name) " * join([@sprintf("%.1f", med(sel(blk[k][:h], reg))) for reg in 1:3], "/") for (k, r) in enumerate(runs)], ";  "), "  (shelf/slope/deep)")

println("\n   By day, deep water, lower half of the ML (10⁻⁵ m²/s²):")
for (n, day) in enumerate(days)
    @printf("   day %2d", day)
    for (k, d) in enumerate(blk)
        v = [d[:e_low][i, j, n] for i in axes(cbh, 1), j in axes(cbh, 2) if common[i, j, n] && creg[i, j] == 3]
        @printf("  %s %6.2f", runs[k].name, 1e5 * med(v))
    end
    println()
end

println("\n2. Lower half of the mixed layer: the proposed sink per unit C_Δ, G = Δ|∇ₕb|√e, against CATKE's terms (medians)")
@printf("%-20s %-12s %10s %10s %10s %8s %8s %10s %12s\n", "region", "run", "G 1e-9", "κcN² 1e-9", "κuS² 1e-9", "G/κcN²", "G/κuS²",
        "ℓc (m)", "C_Δ max h10")
for (reg, rname) in enumerate(regions), (k, r) in enumerate(runs)
    G = sel(blk[k][:G_low], reg); Bv = sel(blk[k][:B_low], reg); Pv = sel(blk[k][:P_low], reg); ℓ = sel(blk[k][:ℓ_low], reg)
    @printf("%-20s %-12s %10.2f %10.2f %10.2f %8.3f %8.3f %10.3f %12.4f\n", rname[1:min(end, 20)], r.name, 1e9 * med(G), 1e9 * med(Bv),
            1e9 * med(Pv), med(G ./ Bv), med(G ./ Pv), med(ℓ), med(ℓ) / (sqrt(12) * 10))
end
println("   (G at each run's own resolution, averaged onto the common blocks; C_Δ max h10 = ℓc/(√12·10 m), the eq. 2.4 bound for h_χ = 10 m)")

finest = argmax([r.Nx for r in runs])
println("\n3. TKE excess of each coarser run over $(runs[finest].name), lower half of the ML, against the coarser run's G")
for (k, r) in enumerate(runs)
    k == finest && continue
    for (reg, rname) in enumerate(regions)
        ΔE = sel(blk[k][:e_low], reg) .- sel(blk[finest][:e_low], reg)
        G = sel(blk[k][:G_low], reg); E = sel(blk[k][:e_low], reg); grad = sel(blk[k][:∇b_low], reg)
        ok = isfinite.(ΔE) .& isfinite.(G)
        count(ok) < 20 && continue
        q = quantile(grad[ok], (0.25, 0.5, 0.75)); edges = (0.0, q..., Inf)
        @printf("%-12s %-20s rank corr(ΔE, G) %5.2f   by quartile of |∇ₕb|: median ΔE/E", r.name, rname[1:min(end, 20)], rankcor(ΔE[ok], G[ok]))
        for qb in 1:4
            m = ok .& (edges[qb] .<= grad .< edges[qb+1])
            @printf("  Q%d %+5.2f", qb, med(ΔE[m] ./ E[m]))
        end
        println()
    end
end

# ---------------- figure: time-mean lower-ML TKE at each resolution and the ratio to the finest ----------------
λf = collect(λnodes(coarse.grid.underlying_grid, Face())); φf = collect(φnodes(coarse.grid.underlying_grid, Face()))
tmean(A) = [(v = [A[i, j, n] for n in axes(A, 3) if common[i, j, n]]; isempty(v) ? NaN : mean(v)) for i in axes(A, 1), j in axes(A, 2)]
fig = Figure(size = (520 * (length(runs) + 1), 900), fontsize = 16)
ref = tmean(blk[finest][:e_low])
for (k, r) in enumerate(runs)
    E = tmean(blk[k][:e_low])
    ax = Axis(fig[1, k], title = "$(r.name): lower-ML e (log₁₀ m²/s²)", aspect = DataAspect())
    hm = heatmap!(ax, λf, φf, log10.(E); colormap = :viridis, colorrange = (-6, -3), nan_color = :gray80)
    k == length(runs) && Colorbar(fig[1, k+1], hm)
    if k != finest
        ax2 = Axis(fig[2, k], title = "$(r.name) / $(runs[finest].name)", aspect = DataAspect())
        hm2 = heatmap!(ax2, λf, φf, log2.(E ./ ref); colormap = :balance, colorrange = (-2, 2), nan_color = :gray80)
        Colorbar(fig[2, length(runs)+1], hm2; label = "log₂ ratio")
    end
end
Label(fig[0, :], "Time mean over days $(first(days))-$(last(days)), common blocks only", fontsize = 18)
save(OUTP * ".png", fig; px_per_unit = 1)
JLD2.jldsave(OUTP * ".jld2"; days, names = [r.name for r in runs], blocks = blk, common, region = creg)
println("\nsaved ", OUTP, ".png and .jld2")
