#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --exclusive
#SBATCH --time=03:00:00
# distributed_sea_ice.jl with a chosen environment (JENV) and harness copy (HARNESS); CASES as in dsi.sbatch
. /t0/workdir/enrique/mpi05_ib/ib_env.sh
export JULIA_PKG_GC_AUTO=false
JENV=${JENV:-/t0/workdir/enrique/mpi05_ib}
HARNESS=${HARNESS:-/t0/workdir/enrique/runs/seaice/icefix/distributed_sea_ice.jl}
echo "=== node $SLURM_JOB_NODELIST  env $JENV  harness $HARNESS"
grep -A3 '^\[\[deps.NumericalEarth\]\]' $JENV/Manifest.toml | tail -2
IFS=";" read -ra cases <<< "${CASES:-2 2 1 8 immersed+nodynamics 1 300}"
for c in "${cases[@]}"; do
    read -r NP PX PY EVP_SUBSTEPS FEATURES NSTEPS DT <<< "$c"
    FEATURES=${FEATURES//+/,}; export PX PY EVP_SUBSTEPS FEATURES NSTEPS DT
    echo "=== case: $c  FILLICE=${FILLICE:-false}"
    mpirun -np $NP --bind-to none -x LD_LIBRARY_PATH -x PATH -x JULIA_DEPOT_PATH -x MPITRAMPOLINE_LIB -x MPITRAMPOLINE_MPIEXEC -x JULIA_PKG_GC_AUTO \
           -x PX -x PY -x EVP_SUBSTEPS -x FEATURES -x NSTEPS -x DT -x FILLICE -x SHOWDIFF \
           julia -t 4 --project=$JENV $HARNESS 2>&1 | grep -v -e "HostCPUFeatures" -e "CPU info is out-of-date" -e "incorrect CPU feature" -e "mpool.c" -e "incorrect CPU name" -e "Will continue"
done
