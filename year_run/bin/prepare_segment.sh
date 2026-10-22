#!/bin/bash
# Decide the next segment and write $OUT/segment.env (SEGMENT, MAB_TAG, MAB_PICKUP). Exit 3 on a crash loop.
# After a failure at simulated time t_c the new segment picks up from the newest complete checkpoint at or before t_c - 3 days, so
# the 6-day Lanczos daily output (whose filter state is not checkpointed) is valid again by the time t_c is reached again.
. /t0/workdir/enrique/year_run/config.env
. $RUN_ROOT/bin/lib.sh
cd $OUT; mkdir -p logs
prev=$(ls seg*.job 2>/dev/null | sed -E 's/seg([0-9]+)\.job/\1/' | sort -n | tail -1)
if [ -z "$prev" ]; then
    SEGMENT=1; PICKUP=false
else
    SEGMENT=$((prev + 1))
    job=$(cat seg$prev.job); log=$OUT/logs/seg_$job.log
    pick_prev=$(grep -a '^=== segment' $log | head -1 | sed -E 's/.*pickup: ([^ ]+).*/\1/'); [ "$pick_prev" = false ] && pick_prev=0
    tc=$(last_progress_iteration $log); [ -z "$tc" ] && tc=$pick_prev
    [ "$(tail -1 attempts.log 2>/dev/null | cut -d" " -f1)" = "$prev" ] || echo "$prev $pick_prev $tc $(date +%F_%T)" >> attempts.log
    # crash loop: the last two failures did not get a day beyond the failure before them
    if [ $(wc -l < attempts.log) -ge 3 ]; then
        t3=$(tail -1 attempts.log | awk '{print $3}'); t2=$(tail -2 attempts.log | head -1 | awk '{print $3}'); t1=$(tail -3 attempts.log | head -1 | awk '{print $3}')
        if [ $t3 -le $((t2 + DAY_ITERS)) ] && [ $t2 -le $((t1 + DAY_ITERS)) ]; then
            echo "crash loop: failures at iterations $t1 $t2 $t3 made no progress" > STOP; exit 3
        fi
    fi
    target=$((tc - 3 * DAY_ITERS))
    PICKUP=false
    for it in $(complete_iterations); do [ $it -le $target ] && [ $it -gt 0 ] && PICKUP=$it; done
fi
TAG=$OUT/mabyear_s$(printf %02d $SEGMENT)
if [ "$PICKUP" != false ]; then
    # the new segment looks for checkpoints under its own prefix: link the chosen one
    for r in $(seq 0 $((NRANKS - 1))); do
        src=$(ls $OUT/mabyear_s*_checkpoint_rank${r}_iteration${PICKUP}.jld2 | head -1)
        ln -sf $(readlink -f $src) ${TAG}_checkpoint_rank${r}_iteration${PICKUP}.jld2
    done
fi
printf 'SEGMENT=%s\nexport MAB_TAG=%s\nexport MAB_PICKUP=%s\n' $SEGMENT $TAG $PICKUP > segment.env
echo "$(date +%F_%T) prepared segment $SEGMENT pickup $PICKUP" >> supervisor.log
