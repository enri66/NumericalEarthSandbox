#!/bin/bash
. /t0/workdir/enrique/year_run/env/ib_env.sh
export JULIA_PKG_GC_AUTO=false MAB_CELLS_PER_DEGREE=60 MAB_LAND_FRACTION=0.5 MAB_MASK_CLEANUP=true RIVER_DATE=$1
cd /t0/workdir/enrique/year_run/sandbox
julia -t 8 --project=/t0/workdir/enrique/year_run/env scripts/river_check.jl
