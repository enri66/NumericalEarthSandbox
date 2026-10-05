# Runs the open-boundary probe configurations for 10 steps from rest (no restart) and serializes the prognostic
# fields, to compare two code versions:  OBC_TAG=before julia --project=<worktree> obc_continuous_state.jl
include(joinpath(@__DIR__, "obc_probe_models.jl"))
using Serialization
results = Dict{String, Any}()
for ts in (:QuasiAdamsBashforth2, :SplitRungeKutta3), (name, vs, trs) in configurations
    model = probe_model(ts, vs, trs)
    run!(Simulation(model, Δt=10.0, stop_iteration=10, verbose=false))
    results["$ts $name"] = Dict(string(k) => copy(Array(interior(f))) for (k, f) in pairs(Oceananigans.prognostic_fields(model)))
end
serialize("/tmp/obc_continuous_$(get(ENV, "OBC_TAG", "run")).jls", results)
