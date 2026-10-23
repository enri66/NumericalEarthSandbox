using MPI
MPI.Init()
comm = MPI.COMM_WORLD; r = MPI.Comm_rank(comm); n = MPI.Comm_size(comm)
x = fill(Float64(r), 10^7); y = similar(x)           # 80 MB
MPI.Barrier(comm)
t = @elapsed begin
    for _ in 1:3
        MPI.Sendrecv!(x, y, comm; dest = mod(r + 1, n), source = mod(r - 1, n))
    end
end
s = MPI.Allreduce(sum(y), +, comm)
ok = y[1] == Float64(mod(r - 1, n))
MPI.Barrier(comm)
bw = 3 * 80 / t
worst = MPI.Allreduce(bw, min, comm)
r == 0 && println("RESULT n=$n allreduce=$(s) neighbours ok on rank0=$ok  slowest rank bandwidth = $(round(worst, digits=1)) MB/s")
MPI.Finalize()
