# Upper-ocean stratification of a run and of GLORYS against Argo, profile by profile (argo_common.jl): is the model's
# water column easier to mix down than the real one, so that ordinary convective deepening goes further? For each
# profile, from σ₀ on a 1 m grid (constant above the shallowest level):
#   N² below the mixed layer: g/ρ₀ Δσ₀/Δz over the 20 m below each source's own mixed-layer depth (s⁻²);
#   the cooling needed to mix the column down to H = 100 and 200 m: the buoyancy SI(H) = g/ρ₀ ∫₀ᴴ [σ₀(H) − σ₀(z)] dz
#   (m²/s²), expressed as heat ρ₀ cₚ SI / (g α) (MJ/m², with α of the Argo profile's 10 m water), i.e. how much heat
#   the surface must lose before the mixed layer reaches H, if salinity does not change.
# Tables by month (median of each source and median model − Argo), and a figure of the monthly median σ₀(z) − σ₀(10 m)
# profiles. Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/sponge_year_mpi2/spy ARGO_DIR=/t0/workdir/enrique/Data/Argo/mab \
#       julia --project=. scripts/stratification_vs_argo.jl
using CairoMakie
include(joinpath(@__DIR__, "argo_common.jl"))

const TAG      = get(ENV, "MAB_TAG", "mab_H90")
const PREFIX   = run_prefix(TAG)
const ARGO_DIR = get(ENV, "ARGO_DIR", joinpath(homedir(), "Data", "Argo", "mab"))
const g, ρ₀, cₚ = 9.81, 1026.0, 3991.9
const ZGRID = 0.0:1.0:300.0

# σ₀ on ZGRID: linear between levels, constant above the shallowest (if it is no deeper than 10 m), NaN below the deepest
function σ_on_grid(depth, T, S)
    ok = isfinite.(T) .& isfinite.(S)
    d = depth[ok]; σ = σ₀.(T[ok], S[ok])
    out = fill(NaN, length(ZGRID))
    (length(d) < 3 || d[1] > 10) && return out
    for (n, z) in enumerate(ZGRID)
        z > d[end] && break
        k = searchsortedlast(d, z)
        out[n] = k == 0 ? σ[1] : k == length(d) ? σ[k] : σ[k] + (σ[k+1] - σ[k]) * (z - d[k]) / (d[k+1] - d[k])
    end
    return out
end

σ_at(s, z) = (n = round(Int, z) + 1; 1 <= n <= length(s) ? s[n] : NaN)
N²_below(s, mld) = isfinite(mld) ? g / ρ₀ * (σ_at(s, mld + 20) - σ_at(s, mld)) / 20 : NaN
SI(s, H) = (n = round(Int, H) + 1; isfinite(s[n]) ? g / ρ₀ * sum(s[n] .- s[1:n]) * step(ZGRID) : NaN)
α(T, S) = -(σ₀(T + 0.01, S) - σ₀(T - 0.01, S)) / 0.02 / (1000 + σ₀(T, S))

function at10(depth, T, S)
    ok = isfinite.(T) .& isfinite.(S)
    k = findfirst(>=(10), depth[ok])
    isnothing(k) ? (NaN, NaN) : (T[ok][k], S[ok][k])
end

rows = []
for m in matched_argo_profiles(PREFIX, ARGO_DIR)
    a_mld, _ = column_depths(m.argo...)
    isfinite(a_mld) || continue
    T10, S10 = at10(m.argo...)
    heat = ρ₀ * cₚ / (g * α(T10, S10)) / 1e6            # MJ/m² per m²/s² of SI
    r = (; m.time, m.region, m.month, heat)
    for (name, prof) in ((:argo, m.argo), (:model, m.model), (:glorys, m.glorys))
        s = σ_on_grid(prof...)
        mld, _ = column_depths(prof...)
        r = merge(r, NamedTuple{(Symbol(name, :_σ), Symbol(name, :_mld), Symbol(name, :_N²), Symbol(name, :_Q100), Symbol(name, :_Q200))}(
                  (s, mld, N²_below(s, mld), heat * SI(s, 100), heat * SI(s, 200))))
    end
    push!(rows, r)
end

col(rs, name) = [getproperty(r, name) for r in rs]
med(x) = (f = filter(isfinite, x); isempty(f) ? NaN : median(f))
meddiff(x, y) = (f = isfinite.(x) .& isfinite.(y); count(f) == 0 ? NaN : median(x[f] .- y[f]))
months = sort(unique(col(rows, :month)))

for (q, label, scale, unit) in ((:N², "N² in the 20 m below the mixed layer", 1e5, "10⁻⁵ s⁻²"),
                                (:Q100, "cooling needed to mix to 100 m", 1.0, "MJ/m²"),
                                (:Q200, "cooling needed to mix to 200 m", 1.0, "MJ/m²"))
    println("\n$(label) ($(unit)), median by month; model − Argo and GLORYS − Argo are medians of profile differences")
    @printf("%-9s %5s %9s %9s %9s %14s %14s\n", "month", "n", "Argo", "model", "GLORYS", "model − Argo", "GLORYS − Argo")
    for mo in months
        rs = filter(x -> x.month == mo, rows)
        a = scale .* col(rs, Symbol(:argo_, q)); mm = scale .* col(rs, Symbol(:model_, q)); gg = scale .* col(rs, Symbol(:glorys_, q))
        @printf("%-9s %5d %9.2f %9.2f %9.2f %+14.2f %+14.2f\n", mo, count(isfinite, a), med(a), med(mm), med(gg),
                meddiff(mm, a), meddiff(gg, a))
    end
end

# ---------------- figure: monthly median σ₀(z) − σ₀(10 m) ----------------
anomaly(s) = s .- σ_at(s, 10)
ncol = min(length(months), 6)
nrow = cld(length(months), ncol)
fig = Figure(size = (300 * ncol + 100, 520 * nrow), fontsize = 14)
Label(fig[0, 1:ncol], "$(basename(PREFIX)): monthly median σ₀(z) − σ₀(10 m) at Argo profiles (deep water)", fontsize = 18)
for (n, mo) in enumerate(months)
    rs = filter(x -> x.month == mo && x.region == 3, rows)
    ax = Axis(fig[cld(n, ncol), mod1(n, ncol)], title = "$(mo) ($(length(rs)))", xlabel = "Δσ₀ (kg/m³)",
              ylabel = mod1(n, ncol) == 1 ? "depth (m)" : "", yreversed = true, limits = (-0.05, 1.6, 0, 300))
    for (name, color) in ((:argo, :black), (:model, :firebrick), (:glorys, :royalblue))
        isempty(rs) && continue
        A = reduce(hcat, [anomaly(getproperty(r, Symbol(name, :_σ))) for r in rs])
        lines!(ax, [med(A[k, :]) for k in axes(A, 1)], collect(ZGRID); color, label = String(name))
        hlines!(ax, med(col(rs, Symbol(name, :_mld))); color, linestyle = :dash, linewidth = 1)
    end
    vlines!(ax, 0.03; color = :gray, linestyle = :dot)
    n == 1 && axislegend(ax; position = :rb)
end
out = PREFIX * "_stratification_vs_argo.png"
save(out, fig; px_per_unit = 1.2)
println("saved ", out)
