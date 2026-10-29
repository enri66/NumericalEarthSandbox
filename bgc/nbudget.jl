# In-run nitrogen budget of the shelf columns (water depth < MAB_BGC_NBUDGET_DEPTH, default 100 m), included by mab_bgc.jl with
# MAB_BGC_NBUDGET=true. Every rank sums over its own shelf columns; each hour it appends one line to <TAG>_nbudget_rank<r>.txt:
#   time (s), N inventory W (mmol N), and four cumulative integrals (mmol N, since the start of this run):
#     lateral   the advective flux of N out through the shelf's lateral faces (the 100 m isobath, rank-internal faces count
#               only if one side is not shelf, and the open edges), from the velocity and the N concentration at the end of every
#               time step: it includes the tidal part that the daily de-tided output lacks (outward positive)
#     sponge    N added by the restoring sponges of the BGC tracers (positive = gain)
#     scaleneg  N added by ScaleNegativeTracers (the change of the N inventory across the modifier, summed over every stage;
#               positive = gain, which is what zeroing a cell with a negative total does)
# N = NO₃ + NH₄ + P + Z + DOM + sPOM + bPOM. The river input, the sediment and the sinking flux are taken from the log and the
# sediment output by bgc/nitrogen_budget.jl. The totals over the ranks are the sums of the files.
using Oceananigans.Operators: Vᶜᶜᶜ, Axᶠᶜᶜ, Ayᶜᶠᶜ
import Oceananigans.Biogeochemistry: update_biogeochemical_state!

const NB_TRACERS = (:NO₃, :NH₄, :P, :Z, :DOM, :sPOM, :bPOM)
const NB_DEPTH = parse(Float64, get(ENV, "MAB_BGC_NBUDGET_DEPTH", "100"))

mutable struct NBudget
    shelf :: BitMatrix          # (Nx + 2, Ny + 2): this rank's columns and a one-cell border; false outside the domain
    domain :: BitMatrix         # true where the cell is inside the domain
    before :: Float64
    scaleneg :: Float64
    lateral :: Float64
    sponge :: Float64
end

function NBudget(Nxl, Nyl)
    shelf = falses(Nxl + 2, Nyl + 2); domain = falses(Nxl + 2, Nyl + 2)
    for j in 0:Nyl+1, i in 0:Nxl+1
        I, J = I_OFF + i, J_OFF + j
        (1 ≤ I ≤ size(bh, 1) && 1 ≤ J ≤ size(bh, 2)) || continue
        domain[i+1, j+1] = true
        shelf[i+1, j+1] = bh[I, J] < 0 && bh[I, J] > -NB_DEPTH
    end
    return NBudget(shelf, domain, 0.0, 0.0, 0.0, 0.0)
end

@inline function ntot(tr, i, j, k)
    @inbounds tr.NO₃[i, j, k] + tr.NH₄[i, j, k] + tr.P[i, j, k] + tr.Z[i, j, k] + tr.DOM[i, j, k] + tr.sPOM[i, j, k] + tr.bPOM[i, j, k]
end

# N inventory of this rank's shelf columns (mmol N)
function nb_inventory(model, nb::NBudget)
    g = model.grid; Nx, Ny, Nz = size(g); tr = model.tracers
    W = 0.0
    for j in 1:Ny, i in 1:Nx
        nb.shelf[i+1, j+1] || continue
        for k in 1:Nz
            inactive_cell(i, j, k, g) && continue
            W += Vᶜᶜᶜ(i, j, k, g) * ntot(tr, i, j, k)
        end
    end
    return W
end

# net advective flux of N out of the shelf through its lateral faces (mmol N/s, outward positive)
function nb_lateral(model, nb::NBudget)
    g = model.grid; Nx, Ny, Nz = size(g); tr = model.tracers
    u, v = model.velocities.u, model.velocities.v
    out = 0.0
    for j in 1:Ny, i in 1:Nx
        nb.shelf[i+1, j+1] || continue
        for k in 1:Nz
            fe, fw, fn, fs = !nb.shelf[i+2, j+1], !nb.shelf[i, j+1], !nb.shelf[i+1, j+2], !nb.shelf[i+1, j]
            (fe | fw | fn | fs) || continue
            inactive_cell(i, j, k, g) && continue
            c0 = ntot(tr, i, j, k)
            if fe       # east
                if !nb.domain[i+2, j+1]
                    out += u[i+1, j, k] * Axᶠᶜᶜ(i+1, j, k, g) * c0
                elseif !inactive_cell(i+1, j, k, g)
                    out += u[i+1, j, k] * Axᶠᶜᶜ(i+1, j, k, g) * 0.5 * (c0 + ntot(tr, i+1, j, k))
                end
            end
            if fw         # west
                if !nb.domain[i, j+1]
                    out -= u[i, j, k] * Axᶠᶜᶜ(i, j, k, g) * c0
                elseif !inactive_cell(i-1, j, k, g)
                    out -= u[i, j, k] * Axᶠᶜᶜ(i, j, k, g) * 0.5 * (c0 + ntot(tr, i-1, j, k))
                end
            end
            if fn         # north
                if !nb.domain[i+1, j+2]
                    out += v[i, j+1, k] * Ayᶜᶠᶜ(i, j+1, k, g) * c0
                elseif !inactive_cell(i, j+1, k, g)
                    out += v[i, j+1, k] * Ayᶜᶠᶜ(i, j+1, k, g) * 0.5 * (c0 + ntot(tr, i, j+1, k))
                end
            end
            if fs         # south
                if !nb.domain[i+1, j]
                    out -= v[i, j, k] * Ayᶜᶠᶜ(i, j, k, g) * c0
                elseif !inactive_cell(i, j-1, k, g)
                    out -= v[i, j, k] * Ayᶜᶠᶜ(i, j, k, g) * 0.5 * (c0 + ntot(tr, i, j-1, k))
                end
            end
        end
    end
    return out
end

# N added per second by the restoring sponges of the BGC tracers in the shelf columns (mmol N/s)
function nb_sponge(model, nb::NBudget)
    g = model.grid; Nx, Ny, Nz = size(g)
    fields = merge(model.velocities, model.tracers)
    forcings = [(name, getproperty(model.forcing, name)) for name in NB_TRACERS if hasproperty(model.forcing, name)]
    isempty(forcings) && return 0.0
    total = 0.0
    for j in 1:Ny, i in 1:Nx
        nb.shelf[i+1, j+1] || continue
        for k in 1:Nz
            inactive_cell(i, j, k, g) && continue
            s = 0.0
            for (_, f) in forcings
                s += f(i, j, k, g, model.clock, fields)
            end
            total += Vᶜᶜᶜ(i, j, k, g) * s
        end
    end
    return total
end

# modifiers around ScaleNegativeTracers: the inventory before and after it
struct NProbe
    nb :: NBudget
    kind :: Symbol
end
function update_biogeochemical_state!(model, p::NProbe)
    W = nb_inventory(model, p.nb)
    p.kind == :before ? (p.nb.before = W) : (p.nb.scaleneg += W - p.nb.before)
    return nothing
end

function nb_step!(sim, nb::NBudget)
    m = sim.model.ocean.model
    Δt = sim.Δt
    nb.lateral += Δt * nb_lateral(m, nb)
    return nothing
end

# the sponge is evaluated once an hour (its targets and rates change slowly) and counted for the hour
function nb_write!(sim, nb::NBudget, file, interval)
    m = sim.model.ocean.model
    nb.sponge += interval * nb_sponge(m, nb)
    open(file, "a") do io
        println(io, m.clock.time, " ", nb_inventory(m, nb), " ", nb.lateral, " ", nb.sponge, " ", nb.scaleneg)
    end
    return nothing
end
