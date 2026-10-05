# Along the open boundaries, how do the model's water depths compare with TPXO10's? The tidal transports enter the
# model unconverted (TPXO's m²/s), while the GLORYS subtidal transport is GLORYS's depth-mean velocity times the
# model's own depth. Where the model is shallower than TPXO the tide therefore enters with a larger velocity, by
# H_TPXO / H_model. For each boundary: the model depth of each wet boundary cell, TPXO's depth there, the TPXO M2
# normal transport, and the M2-transport-weighted mean of H_TPXO / H_model, with the cells that carry most of the
# tidal transport listed.
# Usage:
#   MAB_TAG=/t0/.../res_test/cd03 TPXO_DIR=/t0/workdir/enrique/Data/TPXO10_atlas_v2_nc julia --project=. scripts/tide_boundary_depths.jl
include(joinpath(@__DIR__, "mab_analysis_common.jl"))
ENV["TPXO_DIR"] = get(ENV, "TPXO_DIR", joinpath(homedir(), "Data", "TPXO10_atlas_v2_nc"))
include(joinpath(@__DIR__, "tpxo.jl"))
using NumericalEarth: TPXO10Atlas
using Oceananigans: tidal_atlas_constants

const TAG    = get(ENV, "MAB_TAG", "cd03")
const PREFIX = run_prefix(TAG)

grid = global_grid(PREFIX * "_barotropic.jld2"); ug = grid.underlying_grid
Nx, Ny = size(ug, 1), size(ug, 2)
λc = collect(λnodes(ug, Center())); φc = collect(φnodes(ug, Center()))
λf = collect(λnodes(ug, Face()));   φf = collect(φnodes(ug, Face()))
B = grid.immersed_boundary.bottom_height; bh = [B[i, j, 1] for i in 1:Nx, j in 1:Ny]

window = load_tpxo((:m2,); λ_bounds = (λf[1] - 1, λf[end] + 1), φ_bounds = (φf[1] - 1, φf[end] + 1), dir = ENV["TPXO_DIR"], verbose = false)

# (name, boundary cells (i, j), boundary face positions, which TPXO transport is normal there)
boundaries = (("north", [(i, Ny) for i in 1:Nx], [(λc[i], φf[end]) for i in 1:Nx], :northward_transport),
              ("south", [(i, 1) for i in 1:Nx],  [(λc[i], φf[1]) for i in 1:Nx],   :northward_transport),
              ("east",  [(Nx, j) for j in 1:Ny], [(λf[end], φc[j]) for j in 1:Ny], :eastward_transport),
              ("west",  [(1, j) for j in 1:Ny],  [(λf[1], φc[j]) for j in 1:Ny],   :eastward_transport))

for (name, cells, faces, component) in boundaries
    wet = findall(c -> bh[c...] < 0, cells)
    isempty(wet) && continue
    U = getproperty(tidal_atlas_constants(TPXO10Atlas(), faces[wet], :M2; dir = ENV["TPXO_DIR"]), component)
    Hm = [-bh[cells[k]...] for k in wet]
    Ht = [first(tpxo_depth(window, faces[k]...)) for k in wet]
    w = abs.(U)
    ok = isfinite.(Ht) .& (Ht .> 0)
    r = Ht[ok] ./ Hm[ok]
    @printf("\n== %s boundary: %d wet cells; M2-transport-weighted mean H_TPXO / H_model = %.2f (median of cells %.2f)\n",
            name, length(wet), sum(w[ok] .* r) / sum(w[ok]), median(r))
    @printf("   share of the M2 transport where the model is shallower than TPXO by more than 20%%: %.0f%%\n",
            100 * sum(w[ok][r .> 1.2]) / sum(w[ok]))
    order = sortperm(w; rev = true)
    println("   largest M2 transports:   position          H_model   H_TPXO   |U_M2| (m²/s)   ratio")
    for k in order[1:min(8, end)]
        λ, φ = faces[wet[k]]
        @printf("                          %7.2f, %5.2f   %6.0f m  %6.0f m   %10.2f   %6.2f\n", λ, φ, Hm[k], Ht[k], w[k], Ht[k] / Hm[k])
    end
end
