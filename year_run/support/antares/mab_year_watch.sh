#!/bin/bash
# Every 10 min from cron: let the supervisor on triton resubmit the year-long MAB run if it is not queued or running, and refresh health.txt.
exec 9> /tmp/mab_year_watch.lock; flock -n 9 || exit 0
[ -f $HOME/mab_year.stop ] && exit 0
echo "$(date +%F_%T)" >> $HOME/mab_year_watch.log
ssh -o ConnectTimeout=30 -o BatchMode=yes triton "/t0/workdir/enrique/year_run/bin/supervise.sh; /t0/workdir/enrique/year_run/bin/health.sh" >> $HOME/mab_year_watch.log 2>&1
tail -n 2000 $HOME/mab_year_watch.log > $HOME/mab_year_watch.log.tmp && mv $HOME/mab_year_watch.log.tmp $HOME/mab_year_watch.log
