# Where would the proposed CATKE restratification sink act, compared with where the mixed layer is too deep? The sink
# (G. Wagner, "Unresolved restratification and the CATKE energy budget", 10 Oct 2026, eq. 2.10) is
#   𝒥_b = −C_Δ Δ |∇ₕb| √e,
# Δ the grid spacing. For a script-05 run made with MAB_CATKE_OUTPUT=true, each Argo profile (argo_common.jl) is paired
# with the model's daily column there, and over the lower half of the model's mixed layer (h/2 to h) we take
#   G = Δ |∇ₕb| √e        (the sink per unit C_Δ; b = −g σ₀/ρ₀ from the de-tided daily T and S, ∇ₕb along levels)
#   κc N²                  CATKE's buoyancy sink
#   κc S²                  a lower bound on CATKE's shear production κu S² (κu is not in the daily output)
#   ℓc = κc / √e           with the covariance bound (eq. 2.4) C_Δ ≤ ℓc / (√12 h_χ)
# from the daily means of catke_diagnostics.jl. κc N², κc S² and ℓc are products of daily means, which overstate them
# where κc is episodic (night-time convection, storms); tke_resolution.jl has the instantaneous values from checkpoints. The profiles are binned by quartile of the lower-mixed-layer |∇ₕb|, per
# month, with the median mixed-layer depth bias against Argo in each bin. Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/catke_fall/diag ARGO_DIR=/t0/workdir/enrique/Data/Argo/mab \
#       julia --project=/t0/workdir/enrique/mpi05_ib scripts/greg_sink_vs_argo.jl
include(joinpath(@__DIR__, "argo_common.jl"))

const TAG      = get(ENV, "MAB_TAG", "diag")
const PREFIX   = run_prefix(TAG)
const ARGO_DIR = get(ENV, "ARGO_DIR", joinpath(homedir(), "Data", "Argo", "mab"))
const g, ρ₀, R = 9.81, 1026.0, 6.371e6

vol = PREFIX * "_volume_daily.jld2"; catke_file = PREFIX * "_catke_daily.jld2"
TV = open_series(vol, "T"; backend = OnDisk()); SV = open_series(vol, "S"; backend = OnDisk())
K = Dict(n => open_series(catke_file, n; backend = OnDisk()) for n in ("κc", "N²", "S²", "e"))
grid = TV.grid; ug = grid.underlying_grid
Nx, Ny, Nz = size(ug)
λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center()))
zc = collect(znodes(ug, Center())); zf = collect(znodes(ug, Face()))
B = grid.immersed_boundary.bottom_height
bh = [B[i, j, 1] for i in 1:Nx, j in 1:Ny]
Δλ = deg2rad(λ[2] - λ[1]); Δφ = deg2rad(φ[2] - φ[1])
dx = [R * cosd(φ[j]) * Δλ for j in 1:Ny]; dy = R * Δφ
catke_times = K["κc"].times

# The daily fields of frame n: T, S (centred on t) and the CATKE means over the two days centred on t
cache = Dict{Int, Any}()
function frame_fields(n)
    haskey(cache, n) && return cache[n]
    length(cache) > 2 && delete!(cache, minimum(keys(cache)))
    t = TV.times[n]
    m₁ = findfirst(τ -> abs(τ - t) < 3600, catke_times); m₂ = findfirst(τ -> abs(τ - t - 86400) < 3600, catke_times)
    (isnothing(m₁) || isnothing(m₂)) && return (cache[n] = nothing)
    T = Array(interior(TV[n])); S = Array(interior(SV[n]))
    A = Dict(name => (Array(interior(K[name][m₁])) .+ Array(interior(K[name][m₂]))) ./ 2 for name in keys(K))
    b = [zc[k] > bh[i, j] ? -g * σ₀(T[i, j, k], S[i, j, k]) / ρ₀ : NaN for i in 1:Nx, j in 1:Ny, k in 1:Nz]
    return (cache[n] = (; T, S, A, b))
end

function gradb(b, i, j, k)
    ok(p, q) = 1 <= p <= Nx && 1 <= q <= Ny && isfinite(b[p, q, k])
    ok(i, j) || return NaN
    gx = ok(i - 1, j) && ok(i + 1, j) ? (b[i+1, j, k] - b[i-1, j, k]) / 2dx[j] :
         ok(i + 1, j) ? (b[i+1, j, k] - b[i, j, k]) / dx[j] : ok(i - 1, j) ? (b[i, j, k] - b[i-1, j, k]) / dx[j] : NaN
    gy = ok(i, j - 1) && ok(i, j + 1) ? (b[i, j+1, k] - b[i, j-1, k]) / 2dy :
         ok(i, j + 1) ? (b[i, j+1, k] - b[i, j, k]) / dy : ok(i, j - 1) ? (b[i, j, k] - b[i, j-1, k]) / dy : NaN
    return hypot(gx, gy)
end

# Means over the lower half of the mixed layer, h/2 to h (cell centres for G and |∇b|, faces for the CATKE terms)
function lower_ml(F, i, j, h)
    Δ = sqrt(dx[j] * dy)
    G = Float64[]; ∇ = Float64[]; Bk = Float64[]; Pc = Float64[]; ℓ = Float64[]
    for k in 1:Nz
        d = -zc[k]
        (h / 2 <= d <= h && zc[k] > bh[i, j]) || continue
        gb = gradb(F.b, i, j, k); e = max(F.A["e"][i, j, k], 0)
        push!(G, Δ * gb * sqrt(e)); push!(∇, gb)
    end
    for k in 2:Nz
        d = -zf[k]
        (h / 2 <= d <= h && zf[k] > bh[i, j]) || continue
        κ = F.A["κc"][i, j, k]; ef = max((F.A["e"][i, j, k] + F.A["e"][i, j, k-1]) / 2, 1e-12)
        push!(Bk, κ * F.A["N²"][i, j, k]); push!(Pc, κ * F.A["S²"][i, j, k]); push!(ℓ, κ / sqrt(ef))
    end
    m(x) = (y = filter(isfinite, x); isempty(y) ? NaN : mean(y))
    return (G = m(G), ∇b = m(∇), B = m(Bk), Pc = m(Pc), ℓ = m(ℓ))
end

matches = filter(m -> isfinite(first(column_depths(m.argo...))), matched_argo_profiles(PREFIX, ARGO_DIR))
rows = []
for m in matches
    F = frame_fields(m.frame)
    isnothing(F) && continue
    h = first(column_depths(m.model...))
    (isfinite(h) && h < 0.9 * -bh[m.i, m.j]) || continue
    v = lower_ml(F, m.i, m.j, h)
    isfinite(v.G) || continue
    push!(rows, (; m.month, m.region, h, argo = first(column_depths(m.argo...)), glorys = first(column_depths(m.glorys...)), v...))
end
col(rs, q) = [getproperty(r, q) for r in rs]
med(x) = (y = filter(isfinite, x); isempty(y) ? NaN : median(y))

println("Lower half of the model's mixed layer at Argo profiles, by quartile of |∇ₕb| there (medians; G, κcN², κcS² in 10⁻⁹ m²/s³;")
println("κcN², κcS² and ℓc from products of daily means, so too large where κc is episodic)")
for (label, sel) in (("all months", r -> true), [(mo, r -> r.month == mo) for mo in sort(unique(col(rows, :month)))]...)
    rs = filter(sel, rows)
    length(rs) < 12 && continue
    q = quantile(col(rs, :∇b), (0.25, 0.5, 0.75)); edges = (0.0, q..., Inf)
    @printf("\n%-10s (%3d profiles)\n   %-26s %4s %9s %9s %8s %9s %9s %8s %8s %8s %11s\n", label, length(rs), "|∇ₕb| (10⁻⁸ s⁻²)", "n",
            "mod−Argo", "GLO−Argo", "h (m)", "G", "κcN²", "κcS²", "G/κcN²", "ℓc (m)", "C_Δmax h10")
    for b in 1:4
        rb = filter(r -> edges[b] <= r.∇b < edges[b+1], rs)
        @printf("   Q%d %8.3f - %-12.3f %4d %+8.1f m %+8.1f m %8.1f %9.3f %9.3f %8.3f %8.3f %8.3f %11.5f\n", b, 1e8 * edges[b],
                1e8 * min(edges[b+1], maximum(col(rs, :∇b))), length(rb),
                med(col(rb, :h) .- col(rb, :argo)), med(col(rb, :glorys) .- col(rb, :argo)), med(col(rb, :h)),
                1e9 * med(col(rb, :G)), 1e9 * med(col(rb, :B)), 1e9 * med(col(rb, :Pc)), med(col(rb, :G) ./ col(rb, :B)),
                med(col(rb, :ℓ)), med(col(rb, :ℓ)) / (sqrt(12) * 10))
    end
end

println("\nBy size of the bias (all months): terciles of model − Argo mixed-layer depth")
bias = col(rows, :h) .- col(rows, :argo)
t = quantile(filter(isfinite, bias), (1 / 3, 2 / 3)); edges = (-Inf, t..., Inf)
for b in 1:3
    rb = rows[(edges[b] .<= bias .< edges[b+1])]
    @printf("   bias %+7.1f to %+7.1f m: %4d profiles, median |∇ₕb| %.3f×10⁻⁸, G %.3f, κcN² %.3f, G/κcN² %.3f, ℓc %.3f m\n",
            max(edges[b], minimum(filter(isfinite, bias))), min(edges[b+1], maximum(filter(isfinite, bias))), length(rb),
            1e8 * med(col(rb, :∇b)), 1e9 * med(col(rb, :G)), 1e9 * med(col(rb, :B)), med(col(rb, :G) ./ col(rb, :B)), med(col(rb, :ℓ)))
end
