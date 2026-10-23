#!/bin/bash
#SBATCH -J mabsetup
#SBATCH -p triton16
#SBATCH -N 25
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH -o /t0/workdir/enrique/year_run/test_sponge/setup_only.log
. /t0/workdir/enrique/year_run/config.env
. $ENVDIR/ib_env.sh
export MAB_STAGE=model MAB_PICKUP=false MAB_TAG=/t0/workdir/enrique/year_run/test_sponge/setup_only
cd $SANDBOX
echo "=== setup-only test nodes: $SLURM_JOB_NODELIST start $(date)"
MPIX=""; for v in $(compgen -e | grep "^MAB_"); do MPIX="$MPIX -x $v"; done
mpirun -np $SLURM_NNODES --map-by ppr:1:node --bind-to none -x LD_LIBRARY_PATH -x PATH -x JULIA_DEPOT_PATH -x MPITRAMPOLINE_LIB -x MPITRAMPOLINE_MPIEXEC $MPIX julia -t $THREADS --project=$ENVDIR $SCRIPT
echo "=== end $(date) status $?"
