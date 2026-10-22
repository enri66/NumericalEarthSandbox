#!/bin/bash
# Resubmit the run if it is not queued or running, unless it is complete or stopped. Safe to call from cron and from a dependency job.
. /t0/workdir/enrique/year_run/config.env
. $RUN_ROOT/bin/lib.sh
exec 9> $OUT/.supervise.lock; flock -n 9 || exit 0
cd $OUT || exit 1
[ -f COMPLETE ] || [ -f STOP ] && exit 0
squeue -h -u $USER -n mabyear 2>/dev/null | grep -q . && exit 0
if [ -n "$(complete_iterations | awk -v e=$END_ITER '$1 >= e')" ]; then touch COMPLETE; echo "$(date +%F_%T) run complete" >> supervisor.log; exit 0; fi
$RUN_ROOT/bin/prune_checkpoints.sh
$RUN_ROOT/bin/prepare_segment.sh || { echo "$(date +%F_%T) prepare failed ($?)" >> supervisor.log; exit 1; }
job=$(sbatch --parsable $RUN_ROOT/bin/segment.sbatch) || { echo "$(date +%F_%T) sbatch failed" >> supervisor.log; exit 1; }
echo "$(date +%F_%T) submitted job $job" >> supervisor.log
# when this job ends for any reason, look again at once (the antares cron is the second line)
sbatch -p triton24 -n 1 -c 1 -J mabyear_sup --dependency=afterany:$job -o $OUT/logs/sup_%j.log --wrap "$RUN_ROOT/bin/supervise.sh" >/dev/null 2>&1
