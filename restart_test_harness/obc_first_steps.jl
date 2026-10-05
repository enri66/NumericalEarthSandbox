include(joinpath(@__DIR__, "obc_probe_models.jl"))
using Serialization
out = Dict{String, Any}()
name, vs, trs = configurations[1]
model = probe_model(:SplitRungeKutta3, vs, trs)
sim = Simulation(model, Δt=10.0, stop_iteration=0, verbose=false)
for n in 1:4
    sim.stop_iteration = n; run!(sim)
    out["u$n"] = copy(Array(interior(model.velocities.u))); out["b$n"] = copy(Array(interior(model.tracers.b)))
end
serialize("/tmp/obc_first_$(get(ENV, "OBC_TAG", "run")).jls", out)
