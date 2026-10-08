#!/bin/bash
# The standard evaluation of a script-05 run, on triton: surface maps against GLORYS (SST, SSH, surface currents),
# sections and maps of the mixed-layer and thermocline depths against GLORYS, the MLD bias against fronts, MLD and
# stratification against Argo, and the near-inertial currents against the OOI Pioneer moorings. Each step runs on its own,
# so one failure does not stop the others; the figures and the log are gathered in ANALYSIS_DIR/<tag>, which antares
# pulls into ~/Data/mab_analysis (antares ~/pull_mab_analysis.sh, from cron). Usage (or through analysis_package.sbatch):
#   MAB_TAG=/t0/workdir/enrique/runs/res_test/hr100_fx3 MAB_START_DATE=2019-08-29 MAB_DAYS=5,10,20 \
#       bash scripts/analysis_package.sh
# MAB_TAG:        the run's absolute path prefix
# MAB_START_DATE: the run's start date
# MAB_DAYS:       days for the maps and sections (default: the last daily frame)
# ANALYSIS_DIR:   where the figures are gathered (default /t0/workdir/enrique/runs/analysis)
# STEPS:          which steps to run (default all): surface,transects,fronts,argo_mld,argo_strat,moorings
. /t0/workdir/enrique/julia_env.sh
cd "$(dirname "$0")/.."
: "${MAB_TAG:?set MAB_TAG to the run's path prefix}"
export MAB_TAG MAB_START_DATE=${MAB_START_DATE:-2019-08-29} MAB_DATA_DIR=${MAB_DATA_DIR:-/t0/workdir/enrique/Data/NumericalEarth}
[ -n "$MAB_DAYS" ] && export MAB_DAYS
export ARGO_DIR=${ARGO_DIR:-/t0/workdir/enrique/Data/Argo/mab} OOI_DIR=${OOI_DIR:-/t0/workdir/enrique/Data/OOI/pioneer}
PROJECT=${PROJECT:-/t0/workdir/enrique/mpi05_ib}
STEPS=${STEPS:-surface,transects,fronts,argo_mld,argo_strat,moorings}
tag=$(basename "$MAB_TAG")
out=${ANALYSIS_DIR:-/t0/workdir/enrique/runs/analysis}/$tag
mkdir -p "$out"
log=$out/analysis.log
echo "=== $(date) $tag  days=${MAB_DAYS:-last}  steps=$STEPS  commit $(git log --oneline -1)" | tee -a "$log"

step() {
    local name=$1 script=$2
    [[ ",$STEPS," == *",$name,"* ]] || return 0
    echo "--- $(date +%H:%M) $name ($script)" | tee -a "$log"
    if julia -t ${THREADS:-8} --project="$PROJECT" "scripts/$script" >> "$log" 2>&1; then
        echo "    ok" | tee -a "$log"
    else
        echo "    FAILED (see $log)" | tee -a "$log"
    fi
}

step surface    surface_vs_glorys.jl
step transects  transects_vs_glorys.jl
step fronts     mld_bias_vs_fronts.jl
step argo_mld   mld_vs_argo.jl
step argo_strat stratification_vs_argo.jl
step moorings   moorings_vs_ooi.jl

cp -p "$(dirname "$MAB_TAG")"/"$tag"_*.png "$out"/ 2>/dev/null
echo "=== $(date) done: $(ls "$out"/*.png 2>/dev/null | wc -l) figures in $out" | tee -a "$log"
