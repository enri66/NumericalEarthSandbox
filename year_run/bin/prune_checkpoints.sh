#!/bin/bash
# Keep the newest RETAIN complete checkpoints and every 30th simulated day; remove older ones (real files and symlinks of every
# segment). Nothing at or above the newest complete iteration is touched (it may be being written).
. /t0/workdir/enrique/year_run/config.env
. $RUN_ROOT/bin/lib.sh
cd $OUT || exit 1
all=$( ls mabyear_s*_checkpoint_rank*_iteration*.jld2 2>/dev/null | sed -E 's/.*_iteration([0-9]+)\.jld2/\1/' | sort -un )
newest=$(complete_iterations | tail -1)
[ -z "$newest" ] && exit 0
keep=$( { complete_iterations | tail -n $RETAIN; echo "$all" | awk -v m=$MONTH_ITERS '$1 % m == 0'; } | sort -un )
for it in $all; do
    [ $it -ge $newest ] && continue
    echo "$keep" | grep -qx "$it" && continue
    rm -f mabyear_s*_checkpoint_rank*_iteration${it}.jld2
done
