#!/bin/bash
#SBATCH --exclusive
. /t0/workdir/enrique/year_run/config.env
. $ENVDIR/ib_env.sh
HF=/t0/workdir/enrique/year_run/test_sponge/hfc_$SLURM_JOB_ID; rm -f $HF
for n in $(scontrol show hostnames $SLURM_JOB_NODELIST); do case $n in node5[3-6]) ;; *) echo "$n slots=1" >> $HF ;; esac; done
for n in $(scontrol show hostnames $SLURM_JOB_NODELIST); do case $n in node5[3-6]) echo "$n slots=3" >> $HF ;; esac; done
for tls in "ud,sm,self" "tcp,sm,self"; do
  echo "=== UCX_TLS='$tls'"
  if [ -z "$tls" ]; then X=""; else X="-x UCX_TLS=$tls"; fi
  timeout 400 mpirun -np 40 --hostfile $HF --map-by slot --bind-to none -x LD_LIBRARY_PATH -x PATH -x JULIA_DEPOT_PATH -x MPITRAMPOLINE_LIB -x MPITRAMPOLINE_MPIEXEC $X julia -t 1 --project=$ENVDIR /t0/workdir/enrique/year_run/test_sponge/mpicheck.jl 2>&1 > out/logs/mpicheck_$(echo $tls | tr , _)_$SLURM_JOB_ID.out 2>&1; grep -a -E "RESULT|rror|Abort|signal" out/logs/mpicheck_$(echo $tls | tr , _)_$SLURM_JOB_ID.out | head -4
done
