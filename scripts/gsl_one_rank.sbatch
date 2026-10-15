#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
# NE-Ice: script 05 in the GSL + Newfoundland-Labrador region (MAB_REGION=gsl) on one rank, production settings
# (z*, partial cells, bounds-preserving WENO-7 T/S, SSS restoring, T,S,u,v sponge). From sbatch --export: RUN_TAG,
# STAGE (bathy|model|run), DAYS, THREADS, plus any MAB_* overrides (exported as they are).
. /t0/workdir/enrique/mpi05_ib/ib_env.sh
export JULIA_PKG_GC_AUTO=false
cd /t0/workdir/enrique/Models/NESandbox_ice
echo "=== node $SLURM_JOB_NODELIST  tag $RUN_TAG  stage ${STAGE:-run}"
export MAB_REGION=${REGION:-gsl} MAB_START_DATE=${START:-2013-12-01} MAB_DAYS=${DAYS:-5} MAB_STAGE=${STAGE:-run} \
       MAB_ZSTAR=true MAB_BOTTOM=partial MAB_TRACER_BOUNDS=T:-2:31,S:0:37.2 MAB_SSS_PISTON=${SSS_PISTON:-0.5} \
       MAB_SPONGE_VARS=${SPONGE_VARS:-T,S,u,v} MAB_PROGRESS_MINUTES=${PROGRESS:-60} MAB_TAG=/t0/workdir/enrique/runs/seaice/gsl/$RUN_TAG
# launched through mpirun even on one rank, as production does (a bare julia start segfaulted, 2026-10-15)
mpirun -np 1 --bind-to none -x LD_LIBRARY_PATH -x PATH -x JULIA_DEPOT_PATH -x MPITRAMPOLINE_LIB -x MPITRAMPOLINE_MPIEXEC -x JULIA_PKG_GC_AUTO \
    $(env | grep -o "^MAB_[A-Z_0-9]*" | sed "s/^/-x /" | tr "\n" " ") \
    julia ${JULIA_FLAGS} -t ${THREADS:-48} --project=/t0/workdir/enrique/mpi05_ib scripts/05_mab_glorys_tides_reservoirs_distributed.jl 2>&1 | grep -v -e HostCPUFeatures -e "CPU info is out-of-date" -e "incorrect CPU feature"
