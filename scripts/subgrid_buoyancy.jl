# The unresolved buoyancy variance and the diffusivity-stratification covariance, measured in a fine run by averaging it
# over boxes of r × r cells, against the estimates of G. Wagner, "Unresolved restratification and the CATKE energy
# budget" (10 Oct 2026):
#   χ = ½ (b′)‾  within each box, against the gradient estimate (his eq. 2.7) from the box-mean buoyancy b̄,
#       χ_g = (Δx² (∂x b̄)² + Δy² (∂y b̄)²) / 24,  Δx, Δy the box widths;
#   𝒥_b = −(κ′γ′)‾ = −[(κc N²)‾ − κ̄c N̄²]  at the faces (γ = N²), his covariance correction to the averaged CATKE flux;
#       his closure (eq. 2.10) is 𝒥_b = −C_Δ Δ |∇ₕb̄| √ē, so C_Δ = −𝒥_b / (Δ |∇ₕb̄| √ē) is what the run implies, and
#       C_χ = −𝒥_b / √(2χ ē) with the measured χ (eq. 2.9).
# Bars are means over the r × r cells of a box at one level (boxes with any dry cell at that level are left out).
# Values are averaged over the lower half of the box's mixed layer (h̄/2 to h̄, h̄ the mean mixed-layer depth of its
# columns, mab_analysis_common.jl) and binned by region and by quartile of |∇ₕb̄|. Also the resolved horizontal
# buoyancy gradient |∇ₕb| in the mixed layer at the run's own grid and on each box grid. Reads the run's checkpoints
# (instantaneous T, S, e, κc). Usage:
#   MAB_TAG=/t0/.../res_test/w60 MAB_START_DATE=2020-01-13 MAB_DAYS_LIST=4,5,6 MAB_BOX=2,4,5,10 \
#       julia -t 32 --project=/t0/workdir/enrique/mpi05_ib scripts/subgrid_buoyancy.jl
# MAB_BOX:      box sizes r in cells (default 2,4,5,10; 5 boxes a 1/60° run into 1/12° cells)
# MAB_EDGE_DEG: width of the edge band left out (default 0.75)
include(joinpath(@__DIR__, "mab_analysis_common.jl"))
include(joinpath(@__DIR__, "mab_checkpoints.jl"))

const PREFIX = run_prefix(ENV["MAB_TAG"])
const BOXES = parse.(Int, split(get(ENV, "MAB_BOX", "2,4,5,10"), ","))
const EDGE = parse(Float64, get(ENV, "MAB_EDGE_DEG", "0.75"))
const g, ρ₀, R = 9.81, 1026.0, 6.371e6

blocks = rank_blocks(PREFIX * ".jld2")
ugs = map(p -> JLD2.jldopen(f -> f["serialized/grid"], p).underlying_grid, blocks)
grid = global_grid(PREFIX * ".jld2"); ug = grid.underlying_grid
Nx, Ny, Nz = size(ug)
zc = collect(znodes(ug, Center())); zf = collect(znodes(ug, Face()))
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center()))
Δλ = deg2rad(λ[2] - λ[1]); Δφ = deg2rad(φ[2] - φ[1])
dx = [R * cosd(φ[j]) * Δλ for j in 1:Ny]; dy = R * Δφ
bh = [grid.immersed_boundary.bottom_height[i, j, 1] for i in 1:Nx, j in 1:Ny]
ckdays = checkpoint_days(PREFIX)
days = haskey(ENV, "MAB_DAYS_LIST") ? parse.(Int, split(ENV["MAB_DAYS_LIST"], ",")) : sort(filter(>(0), collect(keys(ckdays))))
depth_c = reverse(-zc)
regions = ("shelf (<200 m)", "slope (200-1000 m)", "deep (>1000 m)")
region_of(d) = d < 200 ? 1 : d < 1000 ? 2 : 3
med(x) = (y = filter(isfinite, x); isempty(y) ? NaN : median(y))
q90(x) = (y = filter(isfinite, x); isempty(y) ? NaN : quantile(y, 0.9))
ranks(x) = (r = similar(x); r[sortperm(x)] = 1:length(x); r)
rankcor(x, y) = (k = isfinite.(x) .& isfinite.(y); count(k) < 10 ? NaN : cor(ranks(x[k]), ranks(y[k])))

# |∇ₕb| along levels on a grid of spacing (ddx[j], ddy), centred differences, one-sided where a neighbour is missing
function gradient_magnitude(B, ddx, ddy)
    nx, ny, nz = size(B)
    G = fill(NaN, nx, ny, nz)
    Threads.@threads for k in 1:nz
        for j in 1:ny, i in 1:nx
            isfinite(B[i, j, k]) || continue
            ok(p, q) = 1 <= p <= nx && 1 <= q <= ny && isfinite(B[p, q, k])
            gx = ok(i - 1, j) && ok(i + 1, j) ? (B[i+1, j, k] - B[i-1, j, k]) / 2ddx[j] :
                 ok(i + 1, j) ? (B[i+1, j, k] - B[i, j, k]) / ddx[j] : ok(i - 1, j) ? (B[i, j, k] - B[i-1, j, k]) / ddx[j] : NaN
            gy = ok(i, j - 1) && ok(i, j + 1) ? (B[i, j+1, k] - B[i, j-1, k]) / 2ddy :
                 ok(i, j + 1) ? (B[i, j+1, k] - B[i, j, k]) / ddy : ok(i, j - 1) ? (B[i, j, k] - B[i, j-1, k]) / ddy : NaN
            G[i, j, k] = hypot(gx, gy)
        end
    end
    return G
end

# rows[r] collects one named tuple per box and day
rows = Dict(r => Any[] for r in BOXES)
native_grad = Float64[]                                   # |∇ₕb| in the mixed layer at the run's own grid
coarse_grad = Dict(r => Float64[] for r in BOXES)
for day in days
    t0 = time()
    st = read_state(PREFIX, ckdays[day], blocks, ugs)
    b = fill(NaN, Nx, Ny, Nz)
    Threads.@threads for k in 1:Nz
        for j in 1:Ny, i in 1:Nx
            (zc[k] > bh[i, j] && isfinite(st[:T][i, j, k])) || continue
            b[i, j, k] = -g * σ₀(st[:T][i, j, k], st[:S][i, j, k]) / ρ₀
        end
    end
    # mixed-layer depth of every column (NaN where mixed to within 10 % of the bottom)
    h = fill(NaN, Nx, Ny)
    Threads.@threads for j in 1:Ny
        for i in 1:Nx
            bh[i, j] < 0 || continue
            m, _ = column_depths(depth_c, reverse([zc[k] > bh[i, j] ? st[:T][i, j, k] : NaN for k in 1:Nz]),
                                         reverse([zc[k] > bh[i, j] ? st[:S][i, j, k] : NaN for k in 1:Nz]))
            isfinite(m) && m < 0.9 * -bh[i, j] && (h[i, j] = m)
        end
    end
    # N² and the κc N² product at faces k = 2..Nz (face k between centres k - 1 and k)
    N² = fill(NaN, Nx, Ny, Nz + 1)
    for k in 2:Nz
        N²[:, :, k] .= (b[:, :, k] .- b[:, :, k-1]) ./ (zc[k] - zc[k-1])
    end
    κN² = st[:κc] .* N²

    ∇b = gradient_magnitude(b, dx, dy)
    for j in 1:4:Ny, i in 1:4:Nx, k in 1:Nz             # a sample of the native mixed-layer gradients
        isfinite(h[i, j]) && -zc[k] <= h[i, j] && isfinite(∇b[i, j, k]) && push!(native_grad, ∇b[i, j, k])
    end

    for r in BOXES
        NX, NY = Nx ÷ r, Ny ÷ r
        box(A, I, J, k) = view(A, (I-1)*r+1:I*r, (J-1)*r+1:J*r, k)
        b̄ = fill(NaN, NX, NY, Nz); ē = fill(NaN, NX, NY, Nz); χ = fill(NaN, NX, NY, Nz)
        Threads.@threads for k in 1:Nz
            for J in 1:NY, I in 1:NX
                v = box(b, I, J, k)
                all(isfinite, v) || continue
                m = mean(v); b̄[I, J, k] = m; χ[I, J, k] = mean((v .- m) .^ 2) / 2
                ē[I, J, k] = mean(box(st[:e], I, J, k))
            end
        end
        Dx = [r * dx[(J-1)*r + (r+1)÷2] for J in 1:NY]; Dy = r * dy
        ∇b̄ = gradient_magnitude(b̄, Dx, Dy)
        # gradient estimate of χ (eq. 2.7), with the x and y parts of the gradient taken separately
        χg = fill(NaN, NX, NY, Nz)
        for k in 1:Nz, J in 1:NY, I in 1:NX
            isfinite(b̄[I, J, k]) || continue
            ok(p, q) = 1 <= p <= NX && 1 <= q <= NY && isfinite(b̄[p, q, k])
            (ok(I - 1, J) && ok(I + 1, J) && ok(I, J - 1) && ok(I, J + 1)) || continue
            gx = (b̄[I+1, J, k] - b̄[I-1, J, k]) / 2Dx[J]; gy = (b̄[I, J+1, k] - b̄[I, J-1, k]) / 2Dy
            χg[I, J, k] = (Dx[J]^2 * gx^2 + Dy^2 * gy^2) / 24
        end
        for J in 1:NY, I in 1:NX
            jc = (J-1)*r + (r+1)÷2; ic = (I-1)*r + (r+1)÷2
            (λ[ic] - λ[1] > EDGE && λ[end] - λ[ic] > EDGE && φ[jc] - φ[1] > EDGE && φ[end] - φ[jc] > EDGE) || continue
            hs = view(h, (I-1)*r+1:I*r, (J-1)*r+1:J*r)
            all(isfinite, hs) || continue
            h̄ = mean(hs)
            depth = mean(view(bh, (I-1)*r+1:I*r, (J-1)*r+1:J*r))
            Δ = sqrt(Dx[J] * Dy)
            χs = Float64[]; χgs = Float64[]; es = Float64[]; gs = Float64[]; 𝒥s = Float64[]; Bs = Float64[]
            for k in 1:Nz
                d = -zc[k]
                if h̄ / 2 <= d <= h̄ && isfinite(χ[I, J, k])
                    push!(χs, χ[I, J, k]); push!(χgs, χg[I, J, k]); push!(es, ē[I, J, k]); push!(gs, ∇b̄[I, J, k])
                end
                d <= h̄ && isfinite(∇b̄[I, J, k]) && push!(coarse_grad[r], ∇b̄[I, J, k])
            end
            for k in 2:Nz
                d = -zf[k]
                h̄ / 2 <= d <= h̄ || continue
                κ = box(st[:κc], I, J, k); n = box(N², I, J, k); κn = box(κN², I, J, k)
                all(isfinite, n) || continue
                push!(𝒥s, -(mean(κn) - mean(κ) * mean(n))); push!(Bs, mean(κ) * mean(n))
            end
            (isempty(χs) || isempty(𝒥s)) && continue
            e = max(mean(es), 0.0); G = mean(filter(isfinite, gs))
            push!(rows[r], (; day, region = region_of(-depth), h̄, Δ, χ = mean(χs), χg = mean(filter(isfinite, χgs)),
                            ∇b̄ = G, ē = e, 𝒥 = mean(𝒥s), B = mean(Bs),
                            CΔ = -mean(𝒥s) / (Δ * G * sqrt(e)), Cχ = -mean(𝒥s) / sqrt(2 * mean(χs) * e)))
        end
    end
    @printf("day %d done (%.0f s)\n", day, time() - t0); flush(stdout)
end

col(rs, q) = [getproperty(x, q) for x in rs]
println("\n$(basename(PREFIX)), days $(join(days, ","))")
@printf("resolved |∇ₕb| in the mixed layer (10⁻⁸ s⁻²): native grid median %.3f, 90th percentile %.3f\n",
        1e8 * med(native_grad), 1e8 * q90(native_grad))
for r in BOXES
    @printf("   boxes of %2d × %-2d cells: median %.3f, 90th percentile %.3f\n", r, r, 1e8 * med(coarse_grad[r]), 1e8 * q90(coarse_grad[r]))
end

println("\nLower half of the mixed layer, per box (medians; χ in 10⁻¹⁰ m²/s⁴, 𝒥_b and κ̄N̄² in 10⁻⁹ m²/s³)")
@printf("%-6s %-20s %6s %8s %9s %9s %8s %9s %9s %9s %8s %9s %9s %9s\n", "box", "region", "n", "Δ (km)", "χ", "χ_grad",
        "χg/χ", "rc(χg,χ)", "𝒥_b", "κ̄N̄²", "𝒥_b<0", "C_Δ", "C_χ", "rc(𝒥,G)")
for r in BOXES, (reg, rname) in enumerate(regions)
    rs = filter(x -> x.region == reg, rows[r])
    length(rs) < 20 && continue
    χ = col(rs, :χ); χg = col(rs, :χg); 𝒥 = col(rs, :𝒥)
    G = col(rs, :Δ) .* col(rs, :∇b̄) .* sqrt.(col(rs, :ē))
    @printf("%2d×%-3d %-20s %6d %8.1f %9.3f %9.3f %8.2f %9.2f %9.3f %9.3f %8.2f %9.4f %9.4f %9.2f\n", r, r, rname[1:min(end, 20)],
            length(rs), med(col(rs, :Δ)) / 1e3, 1e10 * med(χ), 1e10 * med(χg), med(χg ./ χ), rankcor(χg, χ), 1e9 * med(𝒥),
            1e9 * med(col(rs, :B)), count(<(0), 𝒥) / length(𝒥), med(col(rs, :CΔ)), med(col(rs, :Cχ)), rankcor(-𝒥, G))
end
println("(χg/χ: gradient estimate over measured variance; 𝒥_b<0: share of boxes with the sign of a TKE sink; C_Δ, C_χ: the",
        " coefficients the run implies; rc: rank correlations, rc(𝒥,G) between −𝒥_b and Δ|∇ₕb̄|√ē)")

println("\nBy quartile of |∇ₕb̄| (deep water and slope together)")
for r in BOXES
    rs = filter(x -> x.region >= 2, rows[r])
    length(rs) < 40 && continue
    gq = quantile(col(rs, :∇b̄), (0.25, 0.5, 0.75)); edges = (0.0, gq..., Inf)
    for qb in 1:4
        rb = filter(x -> edges[qb] <= x.∇b̄ < edges[qb+1], rs)
        @printf("  %2d×%-2d Q%d |∇ₕb̄| %7.3f  χ %8.3f  χg/χ %6.2f  𝒥_b %8.3f  κ̄N̄² %8.3f  𝒥_b<0 %5.2f  C_Δ %8.4f\n", r, r, qb,
                1e8 * med(col(rb, :∇b̄)), 1e10 * med(col(rb, :χ)), med(col(rb, :χg) ./ col(rb, :χ)), 1e9 * med(col(rb, :𝒥)),
                1e9 * med(col(rb, :B)), count(x -> x.𝒥 < 0, rb) / length(rb), med(col(rb, :CΔ)))
    end
end
