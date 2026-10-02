# Which of CATKE's mixing-length branches mixes the top of the pycnocline? Reads a script-05 run made with
# MAB_CATKE_OUTPUT=true (daily means of κc and its shear and convective parts, N², S², see catke_diagnostics.jl) and
# the de-tided daily T, S, and, for deep-water columns (bottom deeper than 1000 m), sorts the faces of each column by
# depth relative to that day's mixed-layer depth h (column_depths in mab_analysis_common.jl):
#   upper mixed layer: 10 m to h/2;   base: h/2 to h;   transition layer: h to h + 20 m;   below: h + 20 to h + 60 m.
# For each layer: medians of κc, κc_shear, κc_convective, the mean of κc (episodic mixing shows there), the share of faces where the convective part is the larger
# (and so sets κc), and Ri = N²/S² from the daily means. Daily table and monthly summary, and a figure of the median
# profiles against depth relative to h, one panel per month. Usage:
#   MAB_TAG=/t0/workdir/enrique/runs/catke_fall/diag julia --project=. scripts/catke_branches.jl
using CairoMakie
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

const TAG    = get(ENV, "MAB_TAG", "diag")
const PREFIX = run_prefix(TAG)
const STRIDE = parse(Int, get(ENV, "CATKE_COLUMN_STRIDE", "2"))   # every STRIDE-th column in x and y

vol = PREFIX * "_volume_daily.jld2"; catke_file = PREFIX * "_catke_daily.jld2"
TV = open_series(vol, "T"; backend = OnDisk()); SV = open_series(vol, "S"; backend = OnDisk())
K = Dict(n => open_series(catke_file, n; backend = OnDisk()) for n in ("κc", "κc_shear", "κc_convective", "N²", "S²"))
grid = TV.grid; ug = grid.underlying_grid
Nx, Ny, Nz = size(ug)
zc = collect(znodes(ug, Center())); zf = collect(znodes(ug, Face()))
B = grid.immersed_boundary.bottom_height
bh = [B[i, j, 1] for i in 1:Nx, j in 1:Ny]
columns = [(i, j) for i in 1:STRIDE:Nx, j in 1:STRIDE:Ny if bh[i, j] < -1000]
depth_c = reverse(-zc); depth_f = -zf                      # depth_f[k] belongs to face k (bottom to top)
layers = (("upper ML", (h -> 10.0, h -> h / 2)), ("ML base", (h -> h / 2, h -> h)),
          ("transition", (h -> h, h -> h + 20)), ("below", (h -> h + 20, h -> h + 60)))
ζbins = -60.0:4.0:60.0                                     # (depth − h) bins for the profiles (m)

# A catke_daily frame at t is the mean over the day before t, and a volume_daily frame is centred on its time t, so a
# volume frame is paired with the mean of the catke frames at t and t + 1 day (the two days centred on t)
catke_times = K["κc"].times
rows = []
profiles = Dict{String, Any}()
for (n, t) in enumerate(TV.times)
    m₁ = findfirst(τ -> abs(τ - t) < 3600, catke_times); m₂ = findfirst(τ -> abs(τ - t - 86400) < 3600, catke_times)
    (isnothing(m₁) || isnothing(m₂)) && continue
    T = Array(interior(TV[n])); S = Array(interior(SV[n]))
    A = Dict(name => (Array(interior(K[name][m₁])) .+ Array(interior(K[name][m₂]))) ./ 2 for name in keys(K))
    date = start_date + Second(round(Int, t))
    month = Dates.format(date, "yyyy-mm")
    acc = Dict(l => Dict(q => Float64[] for q in (:κc, :sh, :cv, :conv_wins, :Ri)) for (l, _) in layers)
    mlds = Float64[]
    prof = get!(profiles, month) do
        Dict(q => [Float64[] for _ in ζbins] for q in (:sh, :cv, :Ri))
    end
    for (i, j) in columns
        h, _ = column_depths(depth_c, reverse(T[i, j, :]), reverse(S[i, j, :]))
        isfinite(h) || continue
        push!(mlds, h)
        for k in 2:Nz                                       # interior faces
            zf[k] > bh[i, j] || continue
            d = depth_f[k]
            κ, sh, cv = A["κc"][i, j, k], A["κc_shear"][i, j, k], A["κc_convective"][i, j, k]
            Ri = A["N²"][i, j, k] / max(A["S²"][i, j, k], 1e-14)
            for (l, (lo, hi)) in layers
                if lo(h) <= d < hi(h)
                    push!(acc[l][:κc], κ); push!(acc[l][:sh], sh); push!(acc[l][:cv], cv)
                    push!(acc[l][:conv_wins], cv > sh); push!(acc[l][:Ri], Ri)
                end
            end
            b = searchsortedlast(ζbins, d - h)
            if 1 <= b <= length(ζbins)
                push!(prof[:sh][b], sh); push!(prof[:cv][b], cv); push!(prof[:Ri][b], Ri)
            end
        end
    end
    push!(rows, (; date, month, mld = median(mlds), acc))
end

med(x) = isempty(x) ? NaN : median(x)
println("$(basename(PREFIX)): $(length(columns)) deep-water columns, $(length(rows)) days")
println("\nDaily: median MLD and, in the transition layer (h to h + 20 m), median κc / κc_shear / κc_convective (m²/s),")
println("the share of faces where the convective part sets κc, and median Ri")
@printf("%-12s %7s %10s %10s %10s %8s %9s\n", "date", "MLD", "κc", "shear", "convect.", "conv %", "Ri")
for r in rows
    a = r.acc["transition"]
    @printf("%-12s %7.1f %10.2e %10.2e %10.2e %8.0f %9.2f\n", Dates.format(r.date, "yyyy-mm-dd"), r.mld,
            med(a[:κc]), med(a[:sh]), med(a[:cv]), 100 * mean(a[:conv_wins]), med(a[:Ri]))
end

println("\nMonthly, by layer: median κc / shear / convective, mean κc (m²/s), convective share (%), median Ri")
for month in sort(unique([r.month for r in rows]))
    rs = filter(r -> r.month == month, rows)
    @printf("%s (%d days, median MLD %.1f m)\n", month, length(rs), median([r.mld for r in rs]))
    for (l, _) in layers
        pool(q) = reduce(vcat, [r.acc[l][q] for r in rs])
        @printf("   %-11s %10.2e %10.2e %10.2e %10.2e %6.0f %9.2f\n", l, med(pool(:κc)), med(pool(:sh)), med(pool(:cv)),
                mean(pool(:κc)), 100 * mean(pool(:conv_wins)), med(pool(:Ri)))
    end
end

# ---------------- figure ----------------
months = sort(collect(keys(profiles)))
fig = Figure(size = (420 * length(months) + 80, 620), fontsize = 14)
Label(fig[0, 1:length(months)], "$(basename(PREFIX)): CATKE tracer diffusivity by branch, median over deep-water faces, vs depth below the mixed-layer base",
      fontsize = 16)
ζc = ζbins .+ step(ζbins) / 2
for (c, month) in enumerate(months)
    p = profiles[month]
    ax = Axis(fig[1, c], title = month, xlabel = "κ (m²/s)", ylabel = c == 1 ? "depth − MLD (m)" : "", xscale = log10,
              yreversed = true, limits = (1e-6, 1e0, ζbins[1], ζbins[end]))
    lines!(ax, [max(med(x), 1e-7) for x in p[:sh]], ζc; color = :royalblue, label = "shear branch")
    lines!(ax, [max(med(x), 1e-7) for x in p[:cv]], ζc; color = :firebrick, label = "convective branch")
    hlines!(ax, 0; color = :black, linestyle = :dash)
    hlines!(ax, 20; color = :gray, linestyle = :dot)
    axr = Axis(fig[1, c], xaxisposition = :top, xlabel = "Ri", xscale = log10, yreversed = true,
               limits = (1e-1, 1e3, ζbins[1], ζbins[end]), yticksvisible = false, yticklabelsvisible = false)
    lines!(axr, [clamp(med(x), 1e-1, 1e3) for x in p[:Ri]], ζc; color = :black, linewidth = 1)
    c == 1 && axislegend(ax; position = :lb)
end
out = PREFIX * "_catke_branches.png"
save(out, fig; px_per_unit = 1.2)
println("saved ", out)
