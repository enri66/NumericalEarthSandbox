#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
# 1/12 degree, 1 rank, production settings, 10 timed steps after 3 warm-up steps. ENVD, THR from --export.
. /t0/workdir/enrique/year_run/env/ib_env.sh
cd /t0/workdir/enrique/year_run/sandbox
export MAB_START_DATE=2019-09-01 MAB_DAYS=1 MAB_CELLS_PER_DEGREE=12 MAB_NZ=100 MAB_DZ_TOP=1 MAB_DT=300 MAB_PARTITION_X=1 MAB_PARTITION_Y=1
export MAB_LAND_FRACTION=0.5 MAB_MASK_CLEANUP=true MAB_ZSTAR=true MAB_BOTTOM=partial MAB_TRACER_BOUNDS="T:-2:31,S:0:37.2"
export MAB_SPONGE_VARS=T,S,u,v MAB_SPONGE_WIDTH=8 MAB_N_MATCH=4 MAB_STAGE=profile MAB_STEPS=10 MAB_WARMUP=3 MAB_PROFILE_STEPS=1
export MAB_TAG=/t0/workdir/enrique/year_run/bench/$TAGNAME
echo "node: $(hostname -s) $(grep -m1 model name /proc/cpuinfo | cut -d: -f2) threads $THR env $ENVD"
julia -t $THR --project=$ENVD scripts/05_mab_glorys_tides_reservoirs_distributed.jl
