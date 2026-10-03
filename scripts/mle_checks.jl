# Two checks before writing a mixed-layer eddy (Fox-Kemper) restratification scheme, on a script-05 run's de-tided
# daily T and S:
#
# A. Is the mixed-layer depth bias against Argo larger at stronger mixed-layer fronts, where the scheme acts? Each Argo
#    profile (argo_common.jl) is paired with the model's mixed-layer horizontal buoyancy gradient |∇b| at its cell and
#    day, and the biases (model − Argo, GLORYS − Argo) are binned by quartiles of |∇b|.
#
# B. How strong would the restratification be? For deep-water columns (every other column) on every day, the peak
#    eddy buoyancy flux of Fox-Kemper et al. (2008), w'b' = Cₑ H² |∇b|² / |f| with Cₑ = 0.06, at mid-depth of the
#    mixed layer, and the same with the resolution factor of Fox-Kemper et al. (2011), Δs / L_f, for grid spacing Δs
#    and a front width L_f = 5 km. With CATKE's surface buoyancy flux Jᵇ in the run (MAB_CATKE_OUTPUT=true), both are
#    compared with it over the columns that are losing buoyancy (Jᵇ > 0).
#
# H is the mixed-layer depth of mab_analysis_common.jl, b = -g σ₀ / ρ₀, and ∇b is the centred difference of the
# mixed-layer mean b between neighbouring wet columns (one-sided next to land). Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/catke_fall/ctl ARGO_DIR=/t0/workdir/enrique/Data/Argo/mab \
#       julia --project=. scripts/mle_checks.jl
include(joinpath(@__DIR__, "argo_common.jl"))

const TAG      = get(ENV, "MAB_TAG", "ctl")
const PREFIX   = run_prefix(TAG)
const ARGO_DIR = get(ENV, "ARGO_DIR", joinpath(homedir(), "Data", "Argo", "mab"))
const g, ρ₀, Ω, R = 9.81, 1026.0, 7.292e-5, 6.371e6
const Cₑ, L_f = 0.06, 5e3

vol = PREFIX * "_volume_daily.jld2"
TV = open_series(vol, "T"; backend = OnDisk()); SV = open_series(vol, "S"; backend = OnDisk())
grid = TV.grid; ug = grid.underlying_grid
Nx, Ny, Nz = size(ug)
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center())); zc = collect(znodes(ug, Center()))
B = grid.immersed_boundary.bottom_height
bh = [B[i, j, 1] for i in 1:Nx, j in 1:Ny]
wet = bh .< 0
Δλ = deg2rad(λ[2] - λ[1]); Δφ = deg2rad(φ[2] - φ[1])
dx = [R * cosd(φ[j]) * Δλ for j in 1:Ny]; dy = R * Δφ
f = [2Ω * sind(φ[j]) for j in 1:Ny]
depth_c = reverse(-zc)

# Mixed-layer depth and mixed-layer mean buoyancy of every wet column, and |∇b| of the latter, for one frame
function mixed_layer_fields(n)
    T = Array(interior(TV[n])); S = Array(interior(SV[n]))
    H = fill(NaN, Nx, Ny); bml = fill(NaN, Nx, Ny)
    for i in 1:Nx, j in 1:Ny
        wet[i, j] || continue
        t = reverse(T[i, j, :]); s = reverse(S[i, j, :])
        below = depth_c .> -bh[i, j]; t[below] .= NaN; s[below] .= NaN
        h, _ = column_depths(depth_c, t, s)
        isfinite(h) || continue
        k = findall(m -> depth_c[m] <= h && isfinite(t[m]) && isfinite(s[m]), eachindex(depth_c))
        isempty(k) && continue
        H[i, j] = h; bml[i, j] = -g * mean(σ₀.(t[k], s[k])) / ρ₀
    end
    ∂(a, b, Δ) = isfinite(a) && isfinite(b) ? (b - a) / Δ : NaN
    G = fill(NaN, Nx, Ny)
    for i in 1:Nx, j in 1:Ny
        isfinite(bml[i, j]) || continue
        gx = i > 1 && i < Nx && isfinite(bml[i-1, j]) && isfinite(bml[i+1, j]) ? ∂(bml[i-1, j], bml[i+1, j], 2dx[j]) :
             i < Nx && isfinite(bml[i+1, j]) ? ∂(bml[i, j], bml[i+1, j], dx[j]) : i > 1 ? ∂(bml[i-1, j], bml[i, j], dx[j]) : NaN
        gy = j > 1 && j < Ny && isfinite(bml[i, j-1]) && isfinite(bml[i, j+1]) ? ∂(bml[i, j-1], bml[i, j+1], 2dy) :
             j < Ny && isfinite(bml[i, j+1]) ? ∂(bml[i, j], bml[i, j+1], dy) : j > 1 ? ∂(bml[i, j-1], bml[i, j], dy) : NaN
        G[i, j] = hypot(gx, gy)
    end
    return H, G
end

fk08(H, G, j) = Cₑ * H^2 * G^2 / abs(f[j])
fk11(H, G, j) = fk08(H, G, j) * sqrt(dx[j] * dy) / L_f

# ---------------- A: bias against Argo by front strength ----------------
matches = filter(m -> isfinite(first(column_depths(m.argo...))), matched_argo_profiles(PREFIX, ARGO_DIR))
cache = Dict{Int, Any}()
rows = map(matches) do m
    H, G = get!(() -> mixed_layer_fields(m.frame), cache, m.frame)
    length(cache) > 3 && delete!(cache, minimum(keys(cache)))
    (; m.month, m.region, argo = first(column_depths(m.argo...)), model = first(column_depths(m.model...)),
       glorys = first(column_depths(m.glorys...)), G = G[m.i, m.j], fk = fk08(H[m.i, m.j], G[m.i, m.j], m.j))
end
rows = filter(r -> isfinite(r.G) && isfinite(r.model), rows)
meddiff(x, y) = (k = isfinite.(x) .& isfinite.(y); count(k) == 0 ? NaN : median(x[k] .- y[k]))
col(rs, q) = [getproperty(r, q) for r in rs]

println("\nA. Mixed-layer depth bias against Argo by quartile of the model's mixed-layer |∇b| at the profile")
for (label, sel) in (("all months", r -> true), [(mo, r -> r.month == mo) for mo in sort(unique(col(rows, :month)))]...)
    rs = filter(sel, rows)
    length(rs) < 12 && continue
    q = quantile(col(rs, :G), (0.25, 0.5, 0.75))
    @printf("%-11s (%3d profiles)  quartile: |∇b| range (10⁻⁸ s⁻²)   median model − Argo   median GLORYS − Argo   median FK08 w'b' (10⁻⁹ m²/s³)\n",
            label, length(rs))
    edges = (0.0, q..., Inf)
    for b in 1:4
        rb = filter(r -> edges[b] <= r.G < edges[b+1], rs)
        @printf("              Q%d  %6.2f - %-8.2f %4d   %+8.1f m   %+8.1f m   %8.2f\n", b, 1e8 * edges[b], 1e8 * min(edges[b+1], maximum(col(rs, :G))),
                length(rb), meddiff(col(rb, :model), col(rb, :argo)), meddiff(col(rb, :glorys), col(rb, :argo)), 1e9 * median(col(rb, :fk)))
    end
end

# ---------------- B: size of the restratification ----------------
catke_file = PREFIX * "_catke_daily.jld2"
have_Jᵇ = isfile(catke_file) || any(f -> startswith(basename(f), basename(PREFIX) * "_catke_daily_rank"), readdir(dirname(PREFIX)))
if have_Jᵇ
    Jᵇ_frames, Jᵇ_times = surface_series(catke_file, "Jᵇ")
end
deep = [(i, j) for i in 1:2:Nx, j in 1:2:Ny if bh[i, j] < -1000]
println("\nB. Fox-Kemper peak eddy buoyancy flux in $(length(deep)) deep-water columns (Cₑ = $Cₑ; 2011 factor Δs/L_f with L_f = $(L_f / 1e3) km, Δs ≈ $(round(sqrt(dx[Ny ÷ 2] * dy) / 1e3, digits = 1)) km)")
@printf("%-12s %7s %12s %14s %14s %14s %16s %16s\n", "date", "H", "|∇b| 1e-8", "FK08 1e-9", "FK11 1e-9", "Jᵇ>0 1e-9",
        "ΣFK08/ΣJᵇ⁺", "ΣFK11/ΣJᵇ⁺")
monthly = Dict{String, Vector{NTuple{3, Float64}}}()
for (n, t) in enumerate(TV.times)
    H, G = mixed_layer_fields(n)
    date = start_date + Second(round(Int, t))
    J = fill(NaN, Nx, Ny)
    if have_Jᵇ
        m₁ = findfirst(τ -> abs(τ - t) < 3600, Jᵇ_times); m₂ = findfirst(τ -> abs(τ - t - 86400) < 3600, Jᵇ_times)
        (isnothing(m₁) || isnothing(m₂)) || (J = (Jᵇ_frames[:, :, m₁] .+ Jᵇ_frames[:, :, m₂]) ./ 2)
    end
    vals = [(H[i, j], G[i, j], fk08(H[i, j], G[i, j], j), fk11(H[i, j], G[i, j], j), J[i, j]) for (i, j) in deep if isfinite(G[i, j])]
    isempty(vals) && continue
    cooling = filter(v -> isfinite(v[5]) && v[5] > 0, vals)
    r08 = isempty(cooling) ? NaN : sum(getindex.(cooling, 3)) / sum(getindex.(cooling, 5))
    r11 = isempty(cooling) ? NaN : sum(getindex.(cooling, 4)) / sum(getindex.(cooling, 5))
    jm = isempty(cooling) ? NaN : median(getindex.(cooling, 5))
    if n % 3 == 0
        @printf("%-12s %7.1f %12.2f %14.2f %14.2f %14.2f %16.2f %16.2f\n", Dates.format(date, "yyyy-mm-dd"), median(getindex.(vals, 1)),
                1e8 * median(getindex.(vals, 2)), 1e9 * median(getindex.(vals, 3)), 1e9 * median(getindex.(vals, 4)), 1e9 * jm, r08, r11)
    end
    push!(get!(monthly, Dates.format(date, "yyyy-mm"), NTuple{3, Float64}[]), (r08, r11, 1e9 * median(getindex.(vals, 3))))
end
println("\nMonthly medians of the daily values: ΣFK08/ΣJᵇ⁺, ΣFK11/ΣJᵇ⁺ (over cooling columns), median FK08 w'b' (10⁻⁹ m²/s³)")
for mo in sort(collect(keys(monthly)))
    v = monthly[mo]
    med(x) = (y = filter(isfinite, x); isempty(y) ? NaN : median(y))
    @printf("  %s  %6.2f  %6.2f  %8.2f\n", mo, med(getindex.(v, 1)), med(getindex.(v, 2)), med(getindex.(v, 3)))
end
