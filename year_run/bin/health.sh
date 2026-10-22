#!/bin/bash
# Last progress line of the newest segment log, and an ALERT file if a value leaves its sane range.
. /t0/workdir/enrique/year_run/config.env
cd $OUT || exit 1
log=$(ls -t logs/seg_*.log 2>/dev/null | head -1); [ -z "$log" ] && exit 0
line=$(grep -a -E '^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}  \|u\|' $log | tail -1)
echo "$line" > health.txt
echo "$line" | tr -d '[]|=,∈()s' | awk '{ u=$3; v=$4; emin=$5; emax=$6; smin=$7; smax=$8
  sub(/^u/,"",u); bad=""
  if (u+0 > 6 || v+0 > 6) bad=bad" speed"; if (emin+0 < -5 || emax+0 > 5) bad=bad" eta"; if (smin+0 < -0.1 || smax+0 > 40) bad=bad" salinity"
  if (bad != "") { print strftime("%F_%T")" ALERT:"bad": "$0 >> "ALERT" } }' 2>/dev/null
