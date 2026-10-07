# Does NumericalEarth PR #736 (TimeInterpolatedPotential, filled by `update_barotropic_potential!`) read the atmospheric
# pressure correctly when the pressure is a DISK-BACKED, WINDOWED FieldTimeSeries (the ERA5 mechanism)?
#
# The coupled `time_step!` steps the atmosphere first (its clock ticks to t + Δt and `update_field_time_series!` moves
# the in-memory window to that time), and the hook then interpolates the pressure at t (the snapshot called `previous`)
# and at t + Δt (`next`). If the window has just moved, the snapshot below t's interval can lie one index before the new
# window start. Here the pressure has 11 slices 60 s apart, only 3 in memory, and Δt = 60 s, so the window moves on
# every second step (PWT_DT = 45 puts the evaluation times between snapshots). After each step the potential is compared with one interpolated from a fully in-memory copy.
# Run with and without `--check-bounds=yes` (an out-of-window read is then a BoundsError instead of garbage).
using Oceananigans
using Oceananigans.Units
using Oceananigans.OutputReaders: OnDisk, InMemory
using NumericalEarth
using NumericalEarth.EarthSystemModels: OceanOnlyModel
using NumericalEarth.Atmospheres: PrescribedAtmosphere
using NumericalEarth.Radiations: PrescribedRadiation
using NumericalEarth.Oceans: forcing_barotropic_potential
using Printf

const DIR = mktempdir()
const GRID = LatitudeLongitudeGrid(size = (12, 10, 8), halo = (7, 7, 7), longitude = (-71, -70), latitude = (37, 38), z = (-100, 0))
const TIMES = collect(0:60:600.0)
const Δt = parse(Float64, get(ENV, "PWT_DT", "60"))     # s; 60 puts every evaluation exactly on a snapshot, 45 does not
const STEPS = floor(Int, 560 / Δt)
const PATH = joinpath(DIR, "p_series.jld2")

λc = collect(Oceananigans.Grids.λnodes(GRID, Center())); φc = collect(Oceananigans.Grids.φnodes(GRID, Center()))
slice(n) = [101325 + 100n + 30 * sinpi(2 * (λ + 71)) * cospi(φ - 37) for λ in λc, φ in φc]
P = [slice(n) for n in eachindex(TIMES)]

let f_tmp = Field{Center, Center, Nothing}(GRID)
    f = FieldTimeSeries{Center, Center, Nothing}(GRID, TIMES; backend = OnDisk(), path = PATH, name = "p")
    for n in eachindex(TIMES)
        set!(f_tmp, (λ, φ) -> 101325 + 100n + 30 * sinpi(2 * (λ + 71)) * cospi(φ - 37))
        set!(f, f_tmp, n)
    end
end

p_fts = FieldTimeSeries(PATH, "p"; backend = InMemory(3), architecture = CPU())
ocean = ocean_simulation(GRID; warn = false)
atmosphere = PrescribedAtmosphere(GRID; pressure = p_fts)
set!(atmosphere; u = -3.0, v = 2.0, T = 290.0, q = 0.008)
radiation = PrescribedRadiation(GRID)
set!(radiation; downwelling_shortwave = 100.0, downwelling_longwave = 300.0)
model = OceanOnlyModel(ocean; atmosphere, radiation)

ρᵒᶜ = model.interfaces.ocean_properties.reference_density
Φ = forcing_barotropic_potential(model.ocean)
println("barotropic potential: ", typeof(Φ).name.name, "; ρᵒᶜ = ", ρᵒᶜ)

# Reference: the same pressure interpolated linearly in time from the full in-memory array, divided by ρᵒᶜ
function reference(t)
    n = clamp(searchsortedlast(TIMES, t), 1, length(TIMES) - 1)
    χ = (t - TIMES[n]) / (TIMES[n+1] - TIMES[n])
    return ((1 - χ) .* P[n] .+ χ .* P[n+1]) ./ ρᵒᶜ
end
interior_values(f) = Array(interior(f))[:, :, 1]

@printf("%4s %8s %14s %16s %16s %12s\n", "step", "t (s)", "window start", "max|prev − ref|", "max|next − ref|", "max|u| (m/s)")
all_ok = true
for step in 1:STEPS
    time_step!(model, Δt)
    t₁, t₂ = Float64.(Array(Φ.times))
    e₁ = maximum(abs, interior_values(Φ.previous) .- reference(t₁))
    e₂ = maximum(abs, interior_values(Φ.next) .- reference(t₂))
    ok = max(e₁, e₂) <= 1e-6                      # false for NaN too
    global all_ok &= ok
    u = Array(interior(model.ocean.model.velocities.u))
    @printf("%4d %8.0f %14d %16.3e %16.3e %12.3e%s\n", step, t₁, p_fts.backend.start, e₁, e₂, maximum(abs, u), ok ? "" : "   <-- MISMATCH")
end
println(all_ok ? "RESULT: the potential matches the in-memory reference at every step" : "RESULT: the potential does NOT match the in-memory reference at some steps")
