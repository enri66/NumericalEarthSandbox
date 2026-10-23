#!/bin/bash
. /t0/workdir/enrique/year_run/env/ib_env.sh
export JULIA_PKG_GC_AUTO=false WARM_PARTS=16
for p in $(seq 1 16); do
  WARM_PART=$p julia -t 1 --project=/t0/workdir/enrique/year_run/env /t0/workdir/enrique/year_run/warm_glorys_cache.jl > /t0/workdir/enrique/year_run/out/logs/warm_$p.log 2>&1 &
done
wait
echo ALLWARMDONE
