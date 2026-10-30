# Two summary figures of the BGC work of 2026-10-29/30, made from the text output of the analysis scripts:
#   final_depth_chl.png  surface chlorophyll by water depth, satellite against the model (bgc/validate_bgc.jl, "by water depth"), and
#                        the model/satellite ratio of the reference run and the plankton-rate sensitivity runs
#   final_nbudget.png    nitrogen budget of the shelf (bgc/nbudget_analysis.jl): the unexplained remainder X per day before and after
#                        the sediment-thickness fix, and the mean terms after it
#   julia --project=<env with CairoMakie> bgc/final_summary_figures.jl
# Environment: DEPTH_OUTS = "label=validate_bgc output,label=..." (the first one is the reference run), NB_BEFORE and NB_AFTER
# (outputs of nbudget_analysis.jl), FIG_OUT (directory).
using CairoMakie, Printf, Statistics

const DEPTH_OUTS = [(String(first(split(p, "=", limit = 2))), String(last(split(p, "=", limit = 2)))) for p in split(ENV["DEPTH_OUTS"], ",")]
const NB_BEFORE, NB_AFTER = ENV["NB_BEFORE"], ENV["NB_AFTER"]
const OUT = get(ENV, "FIG_OUT", ".")
mkpath(OUT)

function depth_table(path)
    rows = NamedTuple[]
    inside = false
    for line in eachline(path)
        occursin("by water depth", line) && (inside = true; continue)
        inside || continue
        m = match(r"^\s+(\d+)-(\S+)\s+m\s+columns\s+(\d+)\s+cell-days\s+(\d+)\s+satellite\s+([\d.]+)\s+model\s+([\d.]+)\s+model/satellite\s+([\d.]+)", line)
        isnothing(m) && (isempty(rows) || break; continue)
        hi = m[2] == "Inf" ? Inf : parse(Float64, m[2])
        push!(rows, (lo = parse(Float64, m[1]), hi = hi, sat = parse(Float64, m[5]), model = parse(Float64, m[6]), ratio = parse(Float64, m[7])))
    end
    return rows
end
bin_label(r) = isinf(r.hi) ? ">$(Int(r.lo))" : "$(Int(r.lo))-$(Int(r.hi))"

function budget_table(path)
    days = Float64[]; vals = Vector{Float64}[]
    names = ("dW/dt", "river", "lateral", "return", "sinking", "sponge", "ScaleNeg", "X")
    for line in eachline(path)
        m = match(r"^\s+(\d+\.\d)\s+([+-].*)$", line)
        isnothing(m) && continue
        push!(days, parse(Float64, m[1])); push!(vals, parse.(Float64, split(m[2])))
    end
    return days, reduce(hcat, vals), names
end

# ---- chlorophyll by depth ----
tables = [(label, depth_table(path)) for (label, path) in DEPTH_OUTS]
ref_label, ref = first(tables)
fig = Figure(size = (1200, 480))
labels = bin_label.(ref)
ax = Axis(fig[1, 1], title = "surface chlorophyll by water depth ($ref_label)", xlabel = "water depth (m)", ylabel = "mg Chl/m³",
          xticks = (1:length(ref), labels), xticklabelrotation = π / 6, yscale = log10)
barplot!(ax, (1:length(ref)) .- 0.2, getfield.(ref, :sat); width = 0.4, color = :gray45, label = "satellite (Copernicus L3, 4 km)")
barplot!(ax, (1:length(ref)) .+ 0.2, getfield.(ref, :model); width = 0.4, color = :steelblue, label = "model")
axislegend(ax, position = :rt)
ax2 = Axis(fig[1, 2], title = "model / satellite", xlabel = "water depth (m)", ylabel = "ratio", xticks = (1:length(ref), labels),
           xticklabelrotation = π / 6, yscale = log10)
hlines!(ax2, [1.0]; color = :black, linestyle = :dash)
for (n, (label, rows)) in enumerate(tables)
    scatterlines!(ax2, 1:length(rows), getfield.(rows, :ratio); label, color = Makie.wong_colors()[n], markersize = 9)
end
axislegend(ax2, position = :rb)
save(joinpath(OUT, "final_depth_chl.png"), fig)

# ---- nitrogen budget ----
dB, VB, names = budget_table(NB_BEFORE)
dA, VA, _ = budget_table(NB_AFTER)
iX = findfirst(==("X"), names)
fig = Figure(size = (1200, 480))
ax = Axis(fig[1, 1], title = "unexplained nitrogen change of the shelf water", xlabel = "model day", ylabel = "kmol N/day")
hlines!(ax, [0.0]; color = :black, linewidth = 0.8)
scatterlines!(ax, dB, VB[iX, :]; color = :tomato, label = "before the thickness fix (mean $(@sprintf("%.2e", mean(VB[iX, 3:end]))))")
scatterlines!(ax, dA, VA[iX, :]; color = :steelblue, label = "with the thickness fix (mean $(@sprintf("%.2e", mean(VA[iX, 3:end]))))")
axislegend(ax, position = :rc)
ax2 = Axis(fig[1, 2], title = "mean terms of the budget with the fix (days 59-65)", ylabel = "kmol N/day", xticks = (1:length(names), collect(names)), xticklabelrotation = π / 6)
barplot!(ax2, 1:length(names), vec(mean(VA[:, 3:end]; dims = 2)); color = [v < 0 ? :tomato : :steelblue for v in vec(mean(VA[:, 3:end]; dims = 2))])
save(joinpath(OUT, "final_nbudget.png"), fig)
println("wrote ", joinpath(OUT, "final_depth_chl.png"), " and final_nbudget.png")
