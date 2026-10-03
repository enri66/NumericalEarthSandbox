# Argo profiles in the MAB box, and the matching model and GLORYS columns, for mld_vs_argo.jl and
# stratification_vs_argo.jl. The profiles come from the Ifremer ERDDAP ArgoFloats dataset as monthly CSV files (one
# header row and one units row). Argo levels are used only with good quality flags (1 or 2); delayed-mode and adjusted
# profiles use the adjusted values. Each profile is matched with the model's de-tided daily column and GLORYS's column
# at the nearest model cell and the nearest daily frame.
include(joinpath(@__DIR__, "mab_analysis_common.jl"))

good(flag) = flag in ("1", "2")
parsefloat(s) = (x = tryparse(Float64, s); isnothing(x) ? NaN : x)

function read_argo(dir)
    profiles = Dict{Tuple{String, String, String}, Any}()
    for f in sort(filter(endswith(".csv"), readdir(dir; join = true)))
        lines = readlines(f)
        length(lines) < 3 && continue
        header = split(lines[1], ",")
        col = Dict(name => i for (i, name) in enumerate(header))
        for line in lines[3:end]
            v = split(line, ",")
            good(v[col["position_qc"]]) || continue
            adjusted = v[col["data_mode"]] in ("D", "A")
            p, pq = adjusted ? (v[col["pres_adjusted"]], v[col["pres_adjusted_qc"]]) : (v[col["pres"]], v[col["pres_qc"]])
            T, Tq = adjusted ? (v[col["temp_adjusted"]], v[col["temp_adjusted_qc"]]) : (v[col["temp"]], v[col["temp_qc"]])
            S, Sq = adjusted ? (v[col["psal_adjusted"]], v[col["psal_adjusted_qc"]]) : (v[col["psal"]], v[col["psal_qc"]])
            (good(pq) && good(Tq) && good(Sq)) || continue
            key = (v[col["platform_number"]], v[col["cycle_number"]], v[col["direction"]])
            prof = get!(profiles, key) do
                (time = DateTime(v[col["time"]][1:19]), lat = parsefloat(v[col["latitude"]]), lon = parsefloat(v[col["longitude"]]),
                 p = Float64[], T = Float64[], S = Float64[])
            end
            push!(prof.p, parsefloat(p)); push!(prof.T, parsefloat(T)); push!(prof.S, parsefloat(S))
        end
    end
    out = []
    for (key, prof) in profiles
        o = sortperm(prof.p)
        push!(out, (; key, prof.time, prof.lat, prof.lon, p = prof.p[o], T = prof.T[o], S = prof.S[o]))
    end
    return sort(out; by = x -> x.time)
end

const region_names = ("shelf (<200 m)", "slope (200-1000 m)", "deep (>1000 m)")

# Every Argo profile inside the run's grid and within a day of one of its daily frames, with the model's and GLORYS's
# columns there: depths positive downward, T and S NaN below the bottom. `region` indexes region_names; `i`, `j` and
# `frame` are the model cell and the volume_daily frame used.
function matched_argo_profiles(prefix, argo_dir)
    argo = read_argo(argo_dir)
    @printf("%d Argo profiles with good data, %s to %s\n", length(argo), first(argo).time, last(argo).time)

    vol = prefix * "_volume_daily.jld2"
    TV = open_series(vol, "T"; backend = OnDisk()); SV = open_series(vol, "S"; backend = OnDisk())
    grid = TV.grid; ug = grid.underlying_grid
    Nx, Ny, Nz = size(ug)
    λ = collect(λnodes(ug, Center())); φ = collect(φnodes(ug, Center())); zc = collect(znodes(ug, Center()))
    B = grid.immersed_boundary.bottom_height
    bh = [B[i, j, 1] for i in 1:Nx, j in 1:Ny]
    depth_m = reverse(-zc)
    tdays = TV.times ./ 86400
    region(i, j) = -bh[i, j] < 200 ? 1 : -bh[i, j] < 1000 ? 2 : 3

    out = []
    # The profiles are in time order, so only the current day's model and GLORYS fields are kept
    current = Dict{Symbol, Any}(:n => 0, :day => -1)
    for prof in argo
        (λ[1] <= prof.lon <= λ[end] && φ[1] <= prof.lat <= φ[end]) || continue
        i = argmin(abs.(λ .- prof.lon)); j = argmin(abs.(φ .- prof.lat))
        bh[i, j] < 0 || continue
        d = Dates.value(prof.time - start_date) / 86_400_000          # days since the start, as a real number
        n = argmin(abs.(tdays .- d))
        abs(tdays[n] - d) <= 1 || continue                            # outside the run's daily frames
        if current[:n] != n
            current[:T] = Array(interior(TV[n])); current[:S] = Array(interior(SV[n])); current[:n] = n
        end
        T = current[:T][i, j, :]; S = current[:S][i, j, :]
        below = [zc[k] <= bh[i, j] for k in 1:Nz]
        T[below] .= NaN; S[below] .= NaN
        day = round(Int, tdays[n])
        if current[:day] != day
            lon, lat, dep, GT = glorys_at("thetao", day); _, _, _, GS = glorys_at("so", day)
            current[:glorys] = (; lon, lat, dep, GT, GS); current[:day] = day
        end
        g = current[:glorys]
        gT = glorys_profile(g.GT, g.lon, g.lat, g.dep, prof.lon, prof.lat, -bh[i, j])
        gS = glorys_profile(g.GS, g.lon, g.lat, g.dep, prof.lon, prof.lat, -bh[i, j])
        push!(out, (; prof.time, prof.lat, prof.lon, i, j, frame = n, region = region(i, j), month = Dates.format(prof.time, "yyyy-mm"),
                    argo = (depth = prof.p, T = prof.T, S = prof.S),
                    model = (depth = depth_m, T = reverse(T), S = reverse(S)),
                    glorys = (depth = g.dep, T = gT, S = gS)))
    end
    @printf("%d profiles inside the run's grid and period\n\n", length(out))
    return out
end
