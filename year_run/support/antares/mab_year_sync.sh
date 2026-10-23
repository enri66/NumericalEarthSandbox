#!/bin/bash
# Copy the year-long run output (not the checkpoints) from triton to A3. Output files are appended to while the run goes, so
# partial copies are normal until the next sync; the last sync after the run completes is the one that counts.
exec 9> /tmp/mab_year_sync.lock; flock -n 9 || exit 0
DEST=/Volumes/A3/enrique/MAB_year_2019-04
SRC=triton:/t0/workdir/enrique/year_run
mkdir -p $DEST/provenance
rsync -a --partial --exclude="*checkpoint*" --exclude=".supervise.lock" --exclude="*.job" $SRC/out/ $DEST/ >> $HOME/mab_year_sync.log 2>&1
rsync -a $SRC/config.env $SRC/bin $DEST/provenance/ >> $HOME/mab_year_sync.log 2>&1
echo "$(date +%F_%T) sync done (exit $?)" >> $HOME/mab_year_sync.log
tail -n 2000 $HOME/mab_year_sync.log > $HOME/mab_year_sync.log.tmp && mv $HOME/mab_year_sync.log.tmp $HOME/mab_year_sync.log
