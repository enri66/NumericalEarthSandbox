#!/bin/bash
#SBATCH -J mabyear
#SBATCH --exclusive
#SBATCH -o /t0/workdir/enrique/year_run/out/logs/seg_%j.log
# One segment of the year-long run: a fresh start or a pickup (see prepare_segment.sh), to the end of the year or until a failure.
. /t0/workdir/enrique/year_run/config.env
. $ENVDIR/ib_env.sh
. $OUT/segment.env        # SEGMENT, MAB_TAG, MAB_PICKUP
echo $SLURM_JOB_ID > $OUT/seg$SEGMENT.job
cd $SANDBOX
echo "=== segment $SEGMENT job $SLURM_JOB_ID nodes: $SLURM_JOB_NODELIST  pickup: $MAB_PICKUP  tag: $MAB_TAG  commit: $(git log --oneline -1)"
( while sleep 3600; do $RUN_ROOT/bin/prune_checkpoints.sh; done ) &
PRUNER=$!
MPIX=""
for v in $(compgen -e | grep '^MAB_'); do MPIX="$MPIX -x $v"; done
# one rank per triton16 node, three per triton48 node (node53-56), the triton48 nodes last
HF=$OUT/hostfile_$SLURM_JOB_ID
rm -f $HF
for n in $(scontrol show hostnames $SLURM_JOB_NODELIST); do
    case $n in node5[3-6]) ;; *) echo "$n slots=1" >> $HF ;; esac
done
for n in $(scontrol show hostnames $SLURM_JOB_NODELIST); do
    case $n in node5[3-6]) echo "$n slots=3" >> $HF ;; esac
done
echo "ranks: $NRANKS  hostfile: $(tr '\n' ' ' < $HF)"
mpirun -np $NRANKS --hostfile $HF --map-by slot --bind-to none -x LD_LIBRARY_PATH -x PATH -x JULIA_DEPOT_PATH -x MPITRAMPOLINE_LIB -x MPITRAMPOLINE_MPIEXEC $MPIX \
    julia -t $THREADS --project=$ENVDIR $SCRIPT
status=$?
kill $PRUNER 2>/dev/null
echo "=== segment $SEGMENT finished with exit status $status at $(date)"
exit $status
