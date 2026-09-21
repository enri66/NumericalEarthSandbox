# Longitude-depth (x-z) temperature transects at three fixed latitudes (south/mid/north
# of the domain), from a run's final checkpoint — the daily/surface output writers only
# ever saved the top level, so a full-depth snapshot has to come from the checkpoint's
# complete model state instead of a JLD2Writer file.
# Usage:
#   MAB_TAG=mab_3months julia --project=. scripts/plot_temperature_transects.jl
using NumericalEarth
using NumericalEarth.NestedModels: Interpolated
using CopernicusMarine
using CopernicusClimateDataStore
using Oceananigans
using Oceananigans.OutputWriters: load_checkpoint_state
using Oceananigans.Grids: znodes, λnodes, φnodes
using CairoMakie, Printf

include(joinpath(@__DIR__, "glorys_bathymetry.jl"))

const DATA_DIR   = get(ENV, "MAB_DATA_DIR", joinpath(homedir(), "Data", "NumericalEarth"))
const DATA       = joinpath(homedir(), "Data", "mab_glorys_obc")
const TAG        = get(ENV, "MAB_TAG", "mab_3months")
const CHECKPOINT = get(ENV, "MAB_CHECKPOINT", "") # empty = the tag's highest-iteration checkpoint
const resolution = 1 / 12
const Nz = 40

const n_pad    = 2
const data_λ   = (-76.0, -64.0)
const data_φ   = ( 34.0,  42.0)
const λ_bounds = (data_λ[1] + n_pad*resolution, data_λ[2] - n_pad*resolution)
const φ_bounds = (data_φ[1] + n_pad*resolution, data_φ[2] - n_pad*resolution)
const Nλ = round(Int, (λ_bounds[2] - λ_bounds[1]) / resolution)
const Nφ = round(Int, (φ_bounds[2] - φ_bounds[1]) / resolution)

# Identical to 04_mab_glorys_tides_reservoirs.jl's grid (same cache, so this is fast).
z = ExponentialDiscretization(Nz, -4000, 0; scale = 1400)
grid = LatitudeLongitudeGrid(CPU(); size = (Nλ, Nφ, Nz),
                             longitude = λ_bounds, latitude = φ_bounds, z, halo = (7, 7, 7))
bottom_height = regrid_bathymetry(grid; dataset = ETOPO2022(), height_above_water = 1,
                                  minimum_depth = 10, major_basins = 1, interpolation_passes = 10)
match_boundary_bathymetry!(bottom_height, grid, DATA_DIR;
                           region = BoundingBox(longitude = data_λ, latitude = data_φ),
                           n_match = 4, minimum_depth = 10)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

checkpoint_path = if isempty(CHECKPOINT)
    candidates = filter(f -> startswith(f, basename(TAG) * "_checkpoint_iteration"), readdir(DATA))
    iterations = [parse(Int, match(r"iteration(\d+)\.jld2$", f).captures[1]) for f in candidates]
    joinpath(DATA, candidates[argmax(iterations)])
else
    CHECKPOINT
end
@info "reading" checkpoint_path

state = load_checkpoint_state(checkpoint_path)
day = round(state.model.ocean.model.clock.time / 86400, digits = 1)
T_with_halos = state.model.ocean.model.tracers.T.data
Hx, Hy, Hz = grid.underlying_grid.Hx, grid.underlying_grid.Hy, grid.underlying_grid.Hz
Nλg, Nφg, Nzg = size(grid)
T = T_with_halos[Hx+1:Hx+Nλg, Hy+1:Hy+Nφg, Hz+1:Hz+Nzg]

λc = collect(λnodes(grid.underlying_grid, Center()))
φc = collect(φnodes(grid.underlying_grid, Center()))
zc = collect(znodes(grid.underlying_grid, Center()))
bh = Array(interior(bottom_height))[:, :, 1]

# South / middle / north of the MODEL domain (not the full GLORYS box it's padded from).
φ_targets = (φ_bounds[1] + 0.1*(φ_bounds[2]-φ_bounds[1]),
             (φ_bounds[1] + φ_bounds[2]) / 2,
             φ_bounds[1] + 0.9*(φ_bounds[2]-φ_bounds[1]))
labels = ("south", "middle", "north")

function transect_panel(fig, row, φ_target, label; zlims = nothing)
    j = argmin(abs.(φc .- φ_target))
    Tslice = [zc[k] > bh[i, j] ? T[i, j, k] : NaN for i in 1:Nλg, k in 1:Nzg]

    ax = Axis(fig[row, 1], title = @sprintf("T at φ=%.2f°N (%s), day %.1f", φc[j], label, day),
              xlabel = "longitude (°E)", ylabel = "depth (m)")
    isnothing(zlims) || ylims!(ax, zlims)
    hm = heatmap!(ax, λc, zc, Tslice, colormap = :thermal)
    contour!(ax, λc, zc, Tslice, levels = 0:2:30, color = (:black, 0.3), linewidth = 0.5)
    lines!(ax, λc, [bh[i, j] for i in 1:Nλg], color = :black, linewidth = 2)
    Colorbar(fig[row, 2], hm, label = "°C")
end

fig = Figure(size = (1400, 1000), fontsize = 16)
for (row, (φ_target, label)) in enumerate(zip(φ_targets, labels))
    transect_panel(fig, row, φ_target, label)
end
out = joinpath(DATA, TAG * "_temperature_transects.png")
save(out, fig)
println("saved ", out)

fig_top = Figure(size = (1400, 1000), fontsize = 16)
for (row, (φ_target, label)) in enumerate(zip(φ_targets, labels))
    transect_panel(fig_top, row, φ_target, label; zlims = (-500, 0))
end
out_top = joinpath(DATA, TAG * "_temperature_transects_top500m.png")
save(out_top, fig_top)
println("saved ", out_top)
