#!/bin/bash
. /t0/workdir/enrique/year_run/config.env
. $RUN_ROOT/bin/lib.sh
cd $OUT || exit 1
echo "flags: $(ls COMPLETE STOP ALERT 2>/dev/null | tr '\n' ' ')"; [ -f STOP ] && cat STOP
echo "queue:"; squeue -u $USER -n mabyear,mabyear_sup -o '  %.8i %.10P %.12j %.2t %.10M %R' 2>/dev/null
echo "last progress: $(tail -1 health.txt 2>/dev/null)"
echo "checkpoints (complete iterations): $(complete_iterations | tr '\n' ' ' | tail -c 200)"
echo "attempts:"; cat attempts.log 2>/dev/null | sed 's/^/  /'; echo "supervisor log:"; tail -5 supervisor.log 2>/dev/null | sed 's/^/  /'
[ -f ALERT ] && { echo "ALERTS:"; tail -3 ALERT; }
du -sh $OUT 2>/dev/null | sed 's/^/size: /'
