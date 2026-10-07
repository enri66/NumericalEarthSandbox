# Clean up a wet/land mask on a regional C-grid, and list the cells that need a decision.
#
# The mask is a Bool matrix `wet[i, j]` (true = ocean). `clean_mask` applies, in order:
#   1. water not connected to the open ocean through shared cell faces becomes land (lakes, inland bays, ponds); the
#      open ocean is the wet cells on the boundary sides in `open_sides`, and a connection through a corner does not count;
#   2. dead ends one cell wide (a wet cell with at most one wet neighbour through a face) become land, repeatedly, so a
#      bay one cell wide and several cells long is removed from its closed end;
#   3. forced decisions from an overrides table (`land` or `wet` for a cell or a lon/lat box).
# `mask_candidates` lists what is left that a person should look at: channels one cell wide (land on both sides across the
# channel), one-cell islands (land with no land neighbour), and cells joined to the rest of the water only through a corner.
# `write_candidates` writes them with indices and positions to a CSV for review, and the overrides file uses the same layout:
#   action,i,j          action = land | wet, one line per cell
#   action,lon0,lon1,lat0,lat1   a box (cells whose centres are inside it)
using Printf

const NEIGHBOURS = ((1, 0), (-1, 0), (0, 1), (0, -1))

inside(wet, i, j) = 1 <= i <= size(wet, 1) && 1 <= j <= size(wet, 2)
wet_at(wet, i, j) = inside(wet, i, j) && wet[i, j]

"Number of wet neighbours of (i, j) through a face."
face_neighbours(wet, i, j) = count(d -> wet_at(wet, i + d[1], j + d[2]), NEIGHBOURS)

"Wet cells connected to the open ocean (boundary wet cells on `open_sides`) through faces."
function ocean_connected(wet; open_sides = (:west, :east, :south, :north))
    Nx, Ny = size(wet)
    seen = falses(Nx, Ny)
    stack = Tuple{Int, Int}[]
    seed!(i, j) = (wet[i, j] && !seen[i, j]) && (seen[i, j] = true; push!(stack, (i, j)))
    :west  in open_sides && foreach(j -> seed!(1, j), 1:Ny)
    :east  in open_sides && foreach(j -> seed!(Nx, j), 1:Ny)
    :south in open_sides && foreach(i -> seed!(i, 1), 1:Nx)
    :north in open_sides && foreach(i -> seed!(i, Ny), 1:Nx)
    while !isempty(stack)
        i, j = pop!(stack)
        for d in NEIGHBOURS
            a, b = i + d[1], j + d[2]
            inside(wet, a, b) && wet[a, b] && !seen[a, b] && (seen[a, b] = true; push!(stack, (a, b)))
        end
    end
    return seen
end

"Remove dead-end wet cells (at most one wet face neighbour) until none are left. Returns the number removed."
function remove_dead_ends!(wet; open_cells = falses(size(wet)))
    removed = 0
    while true
        found = [(i, j) for i in axes(wet, 1), j in axes(wet, 2) if wet[i, j] && !open_cells[i, j] && face_neighbours(wet, i, j) <= 1]
        isempty(found) && break
        for (i, j) in found; wet[i, j] = false; end
        removed += length(found)
    end
    return removed
end

"Apply an overrides table (see the header). Returns the number of cells changed."
function apply_overrides!(wet, overrides, λ, φ)
    changed = 0
    for row in overrides
        value = row.action == "wet"
        cells = haskey(row, :i) ? [(row.i, row.j)] :
                [(i, j) for i in axes(wet, 1), j in axes(wet, 2) if row.lon0 <= λ[i] <= row.lon1 && row.lat0 <= φ[j] <= row.lat1]
        for (i, j) in cells
            wet[i, j] != value && (wet[i, j] = value; changed += 1)
        end
    end
    return changed
end

"""
    clean_mask(wet; open_sides, overrides = [], λ = nothing, φ = nothing, say = println) -> cleaned wet mask

Steps 1-3 of the header; reports how many cells each step turned to land.
"""
function clean_mask(wet0; open_sides = (:west, :east, :south, :north), overrides = [], λ = nothing, φ = nothing, say = println)
    wet = copy(wet0)
    connected = ocean_connected(wet; open_sides)
    isolated = count(wet .& .!connected)
    wet .&= connected
    say(@sprintf("mask cleanup: %d wet cells not connected to the open ocean through faces turned to land", isolated))
    onboundary = falses(size(wet))
    :west  in open_sides && (onboundary[1, :] .= true);  :east  in open_sides && (onboundary[end, :] .= true)
    :south in open_sides && (onboundary[:, 1] .= true);  :north in open_sides && (onboundary[:, end] .= true)
    dead = remove_dead_ends!(wet; open_cells = onboundary)
    say(@sprintf("mask cleanup: %d dead-end cells (one wet face neighbour) turned to land", dead))
    isempty(overrides) || say(@sprintf("mask cleanup: %d cells changed by the overrides table", apply_overrides!(wet, overrides, λ, φ)))
    say(@sprintf("mask cleanup: wet cells %d → %d", count(wet0), count(wet)))
    return wet
end

"""
    mask_candidates(wet) -> NamedTuple of vectors of (i, j)

`thin_channels`: wet cells with land on both sides across the channel (west and east, or south and north) and wet cells
along it. `islands`: land cells with no land neighbour through a face. `corner_links`: wet cells whose only link to
another wet cell is through a corner.
"""
function mask_candidates(wet)
    Nx, Ny = size(wet)
    thin = Tuple{Int, Int}[]; islands = Tuple{Int, Int}[]; corner = Tuple{Int, Int}[]
    for j in 1:Ny, i in 1:Nx
        if wet[i, j]
            (!wet_at(wet, i-1, j) && !wet_at(wet, i+1, j) && (wet_at(wet, i, j-1) || wet_at(wet, i, j+1)) && inside(wet, i-1, j) && inside(wet, i+1, j)) && push!(thin, (i, j))
            (!wet_at(wet, i, j-1) && !wet_at(wet, i, j+1) && (wet_at(wet, i-1, j) || wet_at(wet, i+1, j)) && inside(wet, i, j-1) && inside(wet, i, j+1)) && push!(thin, (i, j))
            face = face_neighbours(wet, i, j)
            diag = count(d -> wet_at(wet, i + d[1], j + d[2]), ((1, 1), (1, -1), (-1, 1), (-1, -1)))
            (face == 0 && diag > 0) && push!(corner, (i, j))
        else
            1 < i < Nx && 1 < j < Ny && all(d -> wet[i + d[1], j + d[2]], NEIGHBOURS) && push!(islands, (i, j))
        end
    end
    return (; thin_channels = unique(thin), islands, corner_links = corner)
end

"Write the candidates as a CSV: kind,i,j,lon,lat"
function write_candidates(path, candidates, λ, φ)
    open(path, "w") do io
        println(io, "kind,i,j,lon,lat")
        for (kind, cells) in pairs(candidates), (i, j) in cells
            @printf(io, "%s,%d,%d,%.4f,%.4f\n", kind, i, j, λ[i], φ[j])
        end
    end
    return path
end
