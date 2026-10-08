#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH -N 1
#SBATCH -t 12:00:00
# analysis_package.sh on one node. From sbatch --export: MAB_TAG, MAB_START_DATE, MAB_DAYS, STEPS (see the script), e.g.
#   MAB_TAG=/t0/workdir/enrique/runs/res_test/hr100_fx3 MAB_DAYS=5,10 sbatch -p triton48 -J ana \
#       -o /t0/workdir/enrique/runs/analysis/ana_%j.log --export=ALL scripts/analysis_package.sbatch
cd /t0/workdir/enrique/Models/NumericalEarthSandbox
THREADS=${THREADS:-16} bash scripts/analysis_package.sh
