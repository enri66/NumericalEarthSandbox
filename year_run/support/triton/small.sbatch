#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
# 1/12 degree, 2x2 ranks, 1 day: equivalence test of the sponge change. SCRIPTFILE and TAGNAME from --export.
. /t0/workdir/enrique/year_run/env/ib_env.sh
cd /t0/workdir/enrique/year_run/sandbox
export MAB_START_DATE=2019-09-01 MAB_DAYS=1 MAB_CELLS_PER_DEGREE=12 MAB_NZ=50 MAB_DZ_TOP=2 MAB_DT=300 MAB_PARTITION_X=2 MAB_PARTITION_Y=2
export MAB_LAND_FRACTION=0.5 MAB_MASK_CLEANUP=true MAB_ZSTAR=true MAB_BOTTOM=partial MAB_TRACER_BOUNDS="T:-2:31,S:0:37.2"
export MAB_SPONGE_VARS=T,S,u,v MAB_SPONGE_WIDTH=8 MAB_N_MATCH=4 MAB_CHECKPOINT_EVERY=100 MAB_PROGRESS_MINUTES=360
export MAB_TAG=/t0/workdir/enrique/year_run/test_sponge/$TAGNAME
MPIX=""; for v in $(compgen -e | grep "^MAB_"); do MPIX="$MPIX -x $v"; done
mpirun -np $SLURM_NNODES --map-by ppr:1:node --bind-to none -x LD_LIBRARY_PATH -x PATH -x JULIA_DEPOT_PATH -x MPITRAMPOLINE_LIB -x MPITRAMPOLINE_MPIEXEC $MPIX julia -t 16 --project=/t0/workdir/enrique/year_run/env $SCRIPTFILE
