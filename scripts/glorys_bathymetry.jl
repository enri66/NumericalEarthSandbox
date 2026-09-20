# ==================================================================
# Match the model bathymetry to GLORYS near the open boundaries.
#
# The barotropic transport through a boundary is U = ∫u dz over the full water
# column, so it scales with the local depth H. If the model's H at a boundary
# cell differs from GLORYS's H there, the transport we drive in is wrong even
# when the velocity profile is right — and on the shelf/slope a 10% depth error
# is a 10% transport error, which is enough to spin up a spurious boundary jet.
# (Observed: max|u| = 3.5 m/s pinned to the SW corner with unmatched bathymetry.)
#
# `glorys_deptho_on_grid` reads GLORYS's STATIC `deptho` field (the real sea-floor
# depth — NOT the NaN mask of an inpainted tracer, which gets filled during
# inpainting) and bilinearly interpolates it onto the model cell centres.
#
# `match_boundary_bathymetry!` overwrites `bottom_height` in an `n_match`-cell
# band adjacent to each requested boundary with the GLORYS depth, blending
# smoothly (raised cosine) back to the model bathymetry over the band. Where
# GLORYS is land the model value is kept. `minimum_depth` is respected.
# ==================================================================

using NumericalEarth
using NumericalEarth.DataWrangling: Metadata, BoundingBox
using NumericalEarth.DataWrangling.GLORYS: GLORYSStatic
const NCD = NumericalEarth.DataWrangling.NCDatasets
using Oceananigans
using Oceananigans.Grids: λnodes, φnodes
using Downloads: download
using Printf

"""
    glorys_deptho_on_grid(grid, data_dir; region) -> Matrix

GLORYS static sea-floor depth (positive, metres) bilinearly interpolated to the
`Center, Center` points of `grid`. `NaN` where GLORYS has no ocean.
"""
function glorys_deptho_on_grid(grid, data_dir; region)
    m = Metadata(:depth; dataset = GLORYSStatic(), dir = data_dir, region)
    fname = NumericalEarth.DataWrangling.metadata_filename(m.dataset, m.name, m.dates, m.region)
    path  = joinpath(data_dir, fname)
    isfile(path) || download(m)
    isfile(path) || (p = NumericalEarth.DataWrangling.metadata_path(m); path = p isa AbstractVector ? first(p) : p)

    ds = NCD.Dataset(path)
    dg = Array{Float64}(replace(ds["deptho"][:, :], missing => NaN))
    λg = Array{Float64}(ds["longitude"][:])
    φg = Array{Float64}(ds["latitude"][:])
    close(ds)

    ug = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
    λm = collect(λnodes(ug, Center())); φm = collect(φnodes(ug, Center()))

    @inline function bilerp(A, xs, ys, x, y)
        i = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
        j = clamp(searchsortedlast(ys, y), 1, length(ys) - 1)
        tx = (x - xs[i]) / (xs[i+1] - xs[i]); ty = (y - ys[j]) / (ys[j+1] - ys[j])
        v = ((1-tx)*(1-ty)*A[i,j] + tx*(1-ty)*A[i+1,j] + (1-tx)*ty*A[i,j+1] + tx*ty*A[i+1,j+1])
        # if any corner is land, fall back to nearest finite corner
        if !isfinite(v)
            cs = ((A[i,j],(1-tx)*(1-ty)), (A[i+1,j],tx*(1-ty)), (A[i,j+1],(1-tx)*ty), (A[i+1,j+1],tx*ty))
            fin = filter(c -> isfinite(c[1]), cs)
            v = isempty(fin) ? NaN : argmax(c -> c[2], fin)[1]
        end
        return v
    end

    return [bilerp(dg, λg, φg, x, y) for x in λm, y in φm]
end

"""
    match_boundary_bathymetry!(bottom_height, grid, data_dir;
                               region, sides, n_match = 4, minimum_depth = 10)

Blend the model `bottom_height` toward GLORYS `deptho` within `n_match` cells of
each boundary in `sides`. Returns the GLORYS depth field for inspection.
"""
function match_boundary_bathymetry!(bottom_height, grid, data_dir;
                                    region, sides = (:west, :east, :south, :north),
                                    n_match = 4, minimum_depth = 10)
    Hg = glorys_deptho_on_grid(grid, data_dir; region)
    b  = interior(bottom_height)                    # OffsetArray view, positive-up (z of floor)
    Nx, Ny, _ = size(grid)

    # weight: 1 exactly on the boundary, 0 at n_match cells in (raised cosine)
    w(d) = d >= n_match ? 0.0 : 0.5 * (1 + cos(π * d / n_match))

    edge_dist(i, j) = begin
        d = Inf
        :west  in sides && (d = min(d, i - 1))
        :east  in sides && (d = min(d, Nx - i))
        :south in sides && (d = min(d, j - 1))
        :north in sides && (d = min(d, Ny - j))
        d
    end

    nmod = 0
    for i in 1:Nx, j in 1:Ny
        Hgij = Hg[i, j]
        isfinite(Hgij) || continue                  # GLORYS land: keep model value
        α = w(edge_dist(i, j))
        α == 0 && continue
        Hm = -b[i, j]                                # model depth, positive down
        Hm > 0 || continue                          # model land: keep
        Hnew = max(minimum_depth, (1 - α) * Hm + α * Hgij)
        b[i, j] = -Hnew
        nmod += 1
    end
    @printf("  bathymetry matched to GLORYS in %d cells (%d-cell band, sides=%s)\n",
            nmod, n_match, join(sides, ","))
    return Hg
end
