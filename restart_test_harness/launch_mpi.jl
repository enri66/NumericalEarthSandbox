# Run a script on n MPI ranks with MPI.jl's own mpiexec, e.g. on a laptop:
#   julia --project=<env with Oceananigans and MPI> launch_mpi.jl 2 split_explicit_open_y_distributed.jl
using MPI
n = parse(Int, ARGS[1])
run(`$(mpiexec()) -n $n $(Base.julia_cmd()) --project=$(Base.active_project()) $(joinpath(@__DIR__, ARGS[2]))`)
