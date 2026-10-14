# Reading script 05's per-rank checkpoints (instantaneous model state) as global arrays, for tke_resolution.jl and
# subgrid_buoyancy.jl. Needs mab_analysis_common.jl.
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
