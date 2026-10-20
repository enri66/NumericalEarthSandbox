#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
# sea_ice_obs_vs_model.jl for a GSL run on triton. From sbatch --export: RUN_TAG (default ice12_winter), ICE_OBS_END
# (optional), SANDBOX (a checkout with this script; default the triton sandbox), e.g.
#   RUN_TAG=ice12_winter sbatch -p triton24 --exclude=node[38,41] -J iceobs \
#       -o /t0/workdir/enrique/runs/seaice/obs/iceobs_%j.log --export=ALL scripts/sea_ice_obs.sbatch
. /t0/workdir/enrique/mpi05_ib/ib_env.sh
export JULIA_PKG_GC_AUTO=false
cd ${SANDBOX:-/t0/workdir/enrique/Models/NumericalEarthSandbox}
julia -t 4 --project=/t0/workdir/enrique/mpi05_ib scripts/sea_ice_obs_vs_model.jl /t0/workdir/enrique/runs/seaice/gsl \
    ${RUN_TAG:-ice12_winter} /t0/workdir/enrique/Data/sea_ice_obs /t0/workdir/enrique/Data/NumericalEarth 2>&1 |
    grep -v -e HostCPUFeatures -e "CPU info" -e "incorrect CPU"
