#!/bin/bash
#SBATCH -J mabyear
#SBATCH --exclusive
#SBATCH --no-requeue
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
# node6 came back from the 2026-10-26 power loss with less memory (49 GB): give it an interior rank (rank 8 of the 7 x 4 layout),
# which holds no sponge columns, instead of an edge or corner rank
T16=(); for n in $(scontrol show hostnames $SLURM_JOB_NODELIST); do case $n in node5[3-6]) ;; node6) ;; *) T16+=($n) ;; esac; done
if scontrol show hostnames $SLURM_JOB_NODELIST | grep -qx node6; then T16=("${T16[@]:0:8}" node6 "${T16[@]:8}"); fi
for n in "${T16[@]}"; do echo "$n slots=1" >> $HF; done
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
