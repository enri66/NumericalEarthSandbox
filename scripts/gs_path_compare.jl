# Gulf Stream path of two runs against the satellite, from the daily paths ring_census.jl writes
# (<RING_OUT>/gulf_stream_paths.csv: date, source, longitude, latitude).
#   julia --project=. scripts/gs_path_compare.jl <a.csv> <label a> <b.csv> <label b> <output prefix>
# Over the days and longitudes both runs cover: the monthly mean path latitude and its day-to-day spread for each run, the
# satellite and GLORYS (from run a's file), and each run's daily offset from the satellite, mean and rms over those longitudes.
using CairoMakie, Dates, Printf, Statistics

a_csv, a_label, b_csv, b_label, out = ARGS[1:5]

function read_paths(file)
    P = Dict{Tuple{String, Date, Float64}, Float64}()
    for line in Iterators.drop(eachline(file), 1)
        d, s, λ, φ = split(line, ",")
        P[(String(s), Date(d), round(parse(Float64, λ); digits = 4))] = parse(Float64, φ)
    end
    return P
end
A, B = read_paths(a_csv), read_paths(b_csv)
keyset(P, s) = Set((d, λ) for (src, d, λ) in keys(P) if src == s)
days = sort(collect(intersect(Set(first.(keyset(A, "model"))), Set(first.(keyset(B, "model"))))))
lons = sort(collect(intersect(Set(last.(keyset(A, "model"))), Set(last.(keyset(B, "model"))))))
@printf("%d common days %s..%s, %d longitudes %.2f..%.2f\n", length(days), days[1], days[end], length(lons), lons[1], lons[end])

table(P, s) = [get(P, (s, d, λ), NaN) for d in days, λ in lons]           # day × longitude
runs = [(a_label, table(A, "model")), (b_label, table(B, "model"))]
sat, glorys, satb = table(A, "satellite"), table(A, "GLORYS"), table(B, "satellite")
nanmean(x) = (f = filter(isfinite, x); isempty(f) ? NaN : mean(f))
nanstd(x) = (f = filter(isfinite, x); length(f) < 2 ? NaN : std(f))
@printf("satellite path in run a's box minus run b's box: mean %+.3f°\n", nanmean(sat .- satb))

months = unique(Date.(year.(days), month.(days), 1))
println("\nmonthly mean path latitude (mean / day-to-day std averaged over longitude) and offset from the satellite:")
for m in months
    rows = findall(d -> Date(year(d), month(d), 1) == m, days)
    line = @sprintf("  %s (%2d d)  satellite %.2f/%.2f  GLORYS %.2f", Dates.format(m, "yyyy-mm"), length(rows),
                    nanmean(sat[rows, :]), nanmean([nanstd(sat[rows, k]) for k in eachindex(lons)]), nanmean(glorys[rows, :]))
    for (label, P) in runs
        D = P[rows, :] .- sat[rows, :]
        line *= @sprintf("  | %s %.2f/%.2f offset %+.2f° rms %.2f°", label, nanmean(P[rows, :]),
                         nanmean([nanstd(P[rows, k]) for k in eachindex(lons)]), nanmean(D), sqrt(nanmean(D .^ 2)))
    end
    println(line)
end

colors = Dict(a_label => :firebrick, b_label => :darkorange, "satellite" => :black, "GLORYS" => :royalblue)
fig = Figure(size = (430 * length(months) + 60, 820), fontsize = 14)
for (k, m) in enumerate(months)
    rows = findall(d -> Date(year(d), month(d), 1) == m, days)
    ax = Axis(fig[1, k], title = "$(Dates.format(m, "yyyy-mm")) ($(length(rows)) days)", xlabel = "longitude",
              ylabel = k == 1 ? "Gulf Stream path latitude (°N)" : "")
    for (label, P) in [("satellite", sat), ("GLORYS", glorys), runs...]
        μ = [nanmean(P[rows, j]) for j in eachindex(lons)]; σ = [nanstd(P[rows, j]) for j in eachindex(lons)]
        label == "GLORYS" || band!(ax, lons, μ .- σ, μ .+ σ; color = (colors[label], 0.12))
        lines!(ax, lons, μ; color = colors[label], linewidth = label == "GLORYS" ? 1.5 : 2.5, label = label)
    end
    k == 1 && axislegend(ax; position = :lt, labelsize = 11)
end
ax = Axis(fig[2, 1:length(months)], ylabel = "mean path offset (°)",
          title = @sprintf("daily mean offset from the satellite path, %.1f..%.1f°E", lons[1], lons[end]))
for (label, P) in runs
    lines!(ax, 1:length(days), [nanmean(P[n, :] .- sat[n, :]) for n in eachindex(days)]; color = colors[label], label = "$label − satellite")
end
hlines!(ax, 0; color = :gray)
ticks = [findfirst(d -> Date(year(d), month(d), 1) == m, days) for m in months]
ax.xticks = (ticks, Dates.format.(months, "u d"))
axislegend(ax; position = :lb)
save(out * ".png", fig; px_per_unit = 1.3)
println("saved ", out * ".png")
