#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --exclusive
#SBATCH --time=03:00:00
# NE-Ice: does the coupled ocean + ClimaSeaIce (EVP) model give the serial answer on 2 ranks? (restart_test_harness/distributed_sea_ice.jl)
# EVP_EXTEND=false (exported) exchanges the EVP velocity halos every substep; needed for 120-240 substeps on regional partitions.
# CASES: ";"-separated "NP PX PY EVP_SUBSTEPS FEATURES NSTEPS DT", FEATURES joined with + (sbatch --export splits on commas)
. /t0/workdir/enrique/mpi05_ib/ib_env.sh
export JULIA_PKG_GC_AUTO=false
cd /t0/workdir/enrique/Models/NumericalEarthSandbox
echo "=== node $SLURM_JOB_NODELIST  commit $(git log --oneline -1)"
IFS=";" read -ra cases <<< "${CASES:-2 2 1 8 immersed 4 300;2 1 2 8 immersed 4 300}"
for c in "${cases[@]}"; do
    read -r NP PX PY EVP_SUBSTEPS FEATURES NSTEPS DT <<< "$c"
    FEATURES=${FEATURES//+/,}; export PX PY EVP_SUBSTEPS FEATURES NSTEPS DT
    echo "=== case: $c"
    mpirun -np $NP --bind-to none -x LD_LIBRARY_PATH -x PATH -x JULIA_DEPOT_PATH -x MPITRAMPOLINE_LIB -x MPITRAMPOLINE_MPIEXEC -x JULIA_PKG_GC_AUTO \
           -x PX -x PY -x EVP_EXTEND -x EVP_SUBSTEPS -x FEATURES -x NSTEPS -x DT \
           julia -t 4 --project=/t0/workdir/enrique/mpi05_ib restart_test_harness/distributed_sea_ice.jl 2>&1 | grep -v -e "HostCPUFeatures" -e "CPU info is out-of-date" -e "incorrect CPU feature"
done
