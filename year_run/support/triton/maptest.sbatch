#!/bin/bash
#SBATCH --exclusive
. /t0/workdir/enrique/year_run/config.env
. $ENVDIR/ib_env.sh
HF=/t0/workdir/enrique/year_run/test_sponge/hf_$SLURM_JOB_ID
rm -f $HF
for n in $(scontrol show hostnames $SLURM_JOB_NODELIST); do case $n in node5[3-6]) ;; *) echo "$n slots=1" >> $HF ;; esac; done
for n in $(scontrol show hostnames $SLURM_JOB_NODELIST); do case $n in node5[3-6]) echo "$n slots=3" >> $HF ;; esac; done
mpirun -np $NRANKS --hostfile $HF --map-by slot --bind-to none bash -c "echo \$OMPI_COMM_WORLD_RANK \$(hostname -s)" | sort -n | awk "{printf \"%s:%s \", \$1, \$2} END {print \"\"}"
mpirun -np $NRANKS --hostfile $HF --map-by slot --bind-to none bash -c "hostname -s" | sort | uniq -c | awk "{print \$1}" | sort | uniq -c
