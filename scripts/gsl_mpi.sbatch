#!/bin/bash
#SBATCH --exclusive
# NE-Ice: script 05 in the GSL region with sea ice on NP ranks (--nodes/--ntasks from the sbatch line), the environment
# JENV (the sea-ice fork branches of HANDOFF_SEA_ICE.md for more than one rank). From sbatch --export: RUN_TAG, DAYS, STAGE, SEA_ICE,
# THREADS, PX, PY, JENV, SANDBOX, plus any MAB_* overrides (exported as they are).
. /t0/workdir/enrique/mpi05_ib/ib_env.sh
export JULIA_PKG_GC_AUTO=false
cd ${SANDBOX:-/t0/workdir/enrique/Models/NESandbox_icempi}
NP=${SLURM_NTASKS:-1}
echo "=== nodes $SLURM_JOB_NODELIST  tag $RUN_TAG  ranks $NP (${PX:-$NP} × ${PY:-1})  threads ${THREADS:-48}  env $JENV"
export MAB_REGION=gsl MAB_SEA_ICE=${SEA_ICE:-true} MAB_START_DATE=${START:-2013-12-01} MAB_DAYS=${DAYS:-1} MAB_STAGE=${STAGE:-run} \
       MAB_ZSTAR=true MAB_BOTTOM=partial MAB_TRACER_BOUNDS=T:-2:31,S:0:37.2 MAB_SSS_PISTON=${SSS_PISTON:-0.5} \
       MAB_SPONGE_VARS=${SPONGE_VARS:-T,S,u,v} MAB_PROGRESS_MINUTES=${PROGRESS:-60} \
       MAB_PARTITION_X=${PX:-$NP} MAB_PARTITION_Y=${PY:-1} MAB_TAG=/t0/workdir/enrique/runs/seaice/gslmpi/$RUN_TAG
mpirun -np $NP --map-by ppr:$((NP / SLURM_NNODES)):node --bind-to none -x LD_LIBRARY_PATH -x PATH -x JULIA_DEPOT_PATH \
    -x MPITRAMPOLINE_LIB -x MPITRAMPOLINE_MPIEXEC -x JULIA_PKG_GC_AUTO \
    $(env | grep -o "^MAB_[A-Z_0-9]*" | sed "s/^/-x /" | tr "\n" " ") \
    julia ${JULIA_FLAGS} -t ${THREADS:-48} --project=$JENV scripts/05_mab_glorys_tides_reservoirs_distributed.jl 2>&1 | grep --line-buffered -v -e HostCPUFeatures -e "CPU info is out-of-date" -e "incorrect CPU feature"
