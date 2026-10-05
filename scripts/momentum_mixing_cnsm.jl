# How do the 3D model and a single CATKE column at the Pioneer central mooring (CNSM) take up storm momentum? Compares
# the 3D run's hourly mooring column (MAB_TAG, with κu and e in its mooring output, 29a2092 or later) with the columns
# saved by column_cnsm.jl (COLUMN_TAG_<variant>.jld2), over a window (MIX_START, MIX_END):
#   CATKE's momentum viscosity κu: depth-time sections, and mean profiles in calm and storm periods;
#   the mixed-layer depth; the near-inertial (16-22 h) velocity at each depth, as mean kinetic-energy profiles;
#   the depth where κu falls below 10⁻³ m²/s (how deep surface momentum is mixed within hours);
#   the squared shear of the near-inertial velocity across the mixed-layer base.
# Usage:
#   MAB_TAG=/t0/.../res_test/r12k COLUMN_TAG=/t0/.../res_test/r12 julia --project=. scripts/momentum_mixing_cnsm.jl
using CairoMakie
include(joinpath(@__DIR__, "ooi_common.jl"))

const TAG        = get(ENV, "MAB_TAG", "r12k")
const COLUMN_TAG = get(ENV, "COLUMN_TAG", "r12")
const VARIANTS   = split(get(ENV, "COLUMN_VARIANTS", "catke,catke_relative"), ",")
const W0, W1     = DateTime(get(ENV, "MIX_START", "2019-10-09")), DateTime(get(ENV, "MIX_END", "2019-10-23"))
const STORM      = (DateTime(2019, 10, 16, 12), DateTime(2019, 10, 20))
const CALM       = (DateTime(2019, 10, 11), DateTime(2019, 10, 15))
const axis       = collect(W0:Hour(1):W1)

# ---------------- 3D column ----------------
function read_3d(prefix, name = "CNSM")
    for f in filter(f -> occursin(r"_moorings_rank\d+\.jld2$", f) && startswith(basename(f), basename(prefix) * "_moorings"),
                    readdir(dirname(prefix); join = true))
        out = JLD2.jldopen(f) do file
            haskey(file, name) || return nothing
            t0 = DateTime(file["start_date"])
            t = [round(t0 + Second(round(Int, s)), Hour) for s in file["$name/time"]]
            keep = [true; diff(t) .> Hour(0)]
            depth = file["$name/cell"][5]
            zc = file["z_center"]; zf = file["z_face"]
            kc = findall(z -> z > -depth, zc); kf = findall(z -> z >= -depth - 1e-6, zf)
            (; t = t[keep], zc = zc[kc], zf = zf[kf], u = file["$name/u"][kc, keep], v = file["$name/v"][kc, keep],
               T = file["$name/T"][kc, keep], S = file["$name/S"][kc, keep], κu = file["$name/κu"][kf, keep])
        end
        isnothing(out) || return out
    end
    error("no $name column with κu in $prefix")
end

function read_column(variant)
    JLD2.jldopen(COLUMN_TAG * "_column_cnsm_$(variant).jld2") do file
        t0 = DateTime(file["start_date"])
        zf = file["z_face"]; zc = (zf[1:end-1] .+ zf[2:end]) ./ 2
        (; t = [round(t0 + Second(round(Int, s)), Hour) for s in file["time"]], zc, zf, u = file["u"], v = file["v"],
           T = file["T"], S = file["S"], κu = file["κu"])
    end
end

on_axis(src, A) = (ti = Dict(t => n for (n, t) in enumerate(src.t)); [haskey(ti, t) ? A[k, ti[t]] : NaN for k in axes(A, 1), t in axis])
mld(src, n) = first(column_depths(reverse(-src.zc), reverse(src.T[:, n]), reverse(src.S[:, n])))
win(r) = findall(t -> r[1] <= t < r[2], axis)
medrow(A, cols) = [nanmean(A[k, cols]) for k in axes(A, 1)]

sources = [("3D " * basename(TAG), read_3d(TAG)); [("column " * v, read_column(v)) for v in VARIANTS]]
println("CNSM, $(W0) to $(W1); storm window $(STORM[1]) to $(STORM[2]), calm $(CALM[1]) to $(CALM[2])")
@printf("%-22s %12s %12s %16s %16s %18s %18s\n", "source", "MLD calm", "MLD storm", "κu ML calm", "κu ML storm",
        "κu<1e-3 depth st.", "NI KE top 30 m st.")
results = Dict{String, Any}()
for (label, src) in sources
    κu = on_axis(src, src.κu); u = on_axis(src, src.u); v = on_axis(src, src.v)
    ti = Dict(t => n for (n, t) in enumerate(src.t))
    h = [haskey(ti, t) ? mld(src, ti[t]) : NaN for t in axis]
    ub = reduce(vcat, [bandpass(u[k, :])' for k in axes(u, 1)]); vb = reduce(vcat, [bandpass(v[k, :])' for k in axes(v, 1)])
    ke = 0.5 .* (ub .^ 2 .+ vb .^ 2)
    depth_f = -src.zf; depth_c = -src.zc
    # mean κu inside the mixed layer at each hour
    κml = [nanmean(κu[[d < h[n] for d in depth_f], n]) for n in eachindex(axis)]
    # deepest face where κu ≥ 1e-3
    pen = [(k = findall(x -> isfinite(x) && x >= 1e-3, κu[:, n]); isempty(k) ? 0.0 : maximum(depth_f[k])) for n in eachindex(axis)]
    top = [d <= 30 for d in depth_c]
    keTop = [nanmean(ke[top, n]) for n in eachindex(axis)]
    sw, cw = win(STORM), win(CALM)
    @printf("%-22s %10.1f m %10.1f m %14.3f %16.3f %16.1f m %16.2e\n", label, nanmean(h[cw]), nanmean(h[sw]), nanmean(κml[cw]),
            nanmean(κml[sw]), nanmean(pen[sw]), nanmean(keTop[sw]))
    results[label] = (; depth_f, depth_c, κu, ke, h, κml, pen, keTop, storm_κ = medrow(κu, sw), calm_κ = medrow(κu, cw),
                       storm_ke = medrow(ke, sw))
end

# ---------------- figure ----------------
labels = first.(sources)
td = [Dates.value(t - axis[1]) / 86_400_000 for t in axis]
fig = Figure(size = (1500, 950), fontsize = 13)
Label(fig[0, 1:length(labels)], "CNSM: CATKE momentum viscosity κu, 3D model against single columns", fontsize = 17)
hm = nothing
for (c, label) in enumerate(labels)
    r = results[label]
    ax = Axis(fig[1, c], title = label, yreversed = true, ylabel = c == 1 ? "depth (m)" : "", limits = (nothing, (0, 120)),
              xlabel = "days since $(Dates.format(axis[1], "yyyy-mm-dd"))")
    global hm = heatmap!(ax, td, r.depth_f, log10.(max.(r.κu, 1e-6))'; colormap = :viridis, colorrange = (-5, 0))
    lines!(ax, td, r.h; color = :white, linewidth = 1.5)
end
Colorbar(fig[1, length(labels) + 1], hm; label = "log₁₀ κu (m²/s); white: MLD")
ax1 = Axis(fig[2, 1], title = "mean κu, calm (dashed) and storm (solid)", xscale = log10, yreversed = true, ylabel = "depth (m)",
           xlabel = "κu (m²/s)", limits = ((1e-6, 2), (0, 120)))
ax2 = Axis(fig[2, 2], title = "near-inertial KE during the storm window", xscale = log10, yreversed = true, xlabel = "KE (J/kg)",
           limits = (nothing, (0, 120)))
ax3 = Axis(fig[2, 3], title = "mean κu in the mixed layer (m²/s)", xlabel = "days", yscale = log10)
for (label, color) in zip(labels, (:black, :firebrick, :royalblue, :darkorange))
    r = results[label]
    lines!(ax1, max.(r.storm_κ, 1e-6), r.depth_f; color, label); lines!(ax1, max.(r.calm_κ, 1e-6), r.depth_f; color, linestyle = :dash)
    lines!(ax2, max.(r.storm_ke, 1e-7), r.depth_c; color)
    lines!(ax3, td, max.(r.κml, 1e-5); color)
end
axislegend(ax1; position = :rb)
out = TAG * "_momentum_mixing_cnsm.png"
save(out, fig; px_per_unit = 1.2)
println("saved ", out)
