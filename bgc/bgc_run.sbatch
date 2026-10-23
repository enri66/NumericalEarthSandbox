#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
# MAB biogeochemistry run (bgc/mab_bgc.jl) with the production physics of PRODUCTION_CONFIG.md: z*, partial cells,
# bounds-preserving WENO-7 for T and S, GloFAS rivers, SSS restoring with the river mask, T,S,u,v sponge (and the BGC
# sponge with it). One rank per node, PX × PY ranks. From sbatch --export: RUN_TAG, START (YYYY-MM-DD), DAYS, PX, PY,
# CPD (cells per degree), DT (s), THREADS, SPONGE_W and NMATCH (cells), CKPT (days), and optionally STAGE (model = build
# and stop), PROGRESS (minutes), RIVERS, PCO2, BGC_DATA (woa | copernicus), JENV (Julia environment, default
# mpi05_bgc; mpi05_bgc_cmems has the NumericalEarth with CopernicusBGCMonthly) and SANDBOX (the sandbox copy to run
# from, default /t0/workdir/enrique/Models/NumericalEarthSandbox), BGC_BOUNDS (MAB_BGC_BOUNDS, e.g. "NO₃:0:30,O₂:50:450";
# it has commas, so export it in the shell before sbatch rather than in --export), PICKUP (MAB_PICKUP), RIVER_LOADS
# (MAB_RIVER_LOADS, default true: the river water carries tracers), SCALE_NEGATIVES (MAB_BGC_SCALE_NEGATIVES, default false), NO3_SOURCE (MAB_BGC_NO3_SOURCE: woa takes NO₃ from WOA when BGC_DATA=copernicus), RIVER_CHEMISTRY (MAB_RIVER_CHEMISTRY: rc4uscoast, the
# default, = monthly NO₃, NH₄, DOM, Alk, DIC and O₂ from RC4USCoast in Data/RC4USCoast; or assumed = constant NO₃ from
# RIVER_NO3 (MAB_RIVER_NO3, default "Hudson:50,Delaware:50,other:20" mmol/m³; export it in the shell like BGC_BOUNDS) and
# fresh-water O₂). Bounds that include the river values: BGC_BOUNDS="NO₃:0:150,O₂:50:500,DIC:0:2300".
# Example (1/12°, 2 ranks on triton24):
#   sbatch -p triton24 -N 2 --exclude=node[38-41] -J bgc2 -o runs/bgc/bgc2_%j.log \
#     --export=ALL,RUN_TAG=bgc2,START=2019-09-01,DAYS=3,PX=2,PY=1,CPD=12,DT=300,THREADS=24,SPONGE_W=8,NMATCH=4,CKPT=5 \
#     bgc/bgc_run.sbatch
. /t0/workdir/enrique/mpi05_bgc/ib_env.sh
export JULIA_PKG_GC_AUTO=false
cd ${SANDBOX:-/t0/workdir/enrique/Models/NumericalEarthSandbox}
mkdir -p /t0/workdir/enrique/runs/bgc
JENV=${JENV:-/t0/workdir/enrique/mpi05_bgc}
echo "=== nodes: $SLURM_JOB_NODELIST  tag: $RUN_TAG  env: $JENV  BGC data: ${BGC_DATA:-woa}  1/$CPD deg  dt $DT s  start $START  $DAYS days  commit: $(git log --oneline -1)"
export MAB_START_DATE=$START MAB_DAYS=$DAYS MAB_PARTITION_X=$PX MAB_PARTITION_Y=${PY:-1} MAB_NZ=${NZ:-50} MAB_CELLS_PER_DEGREE=$CPD MAB_DT=$DT \
       MAB_PROGRESS_MINUTES=${PROGRESS:-360} MAB_SPONGE_VARS=T,S,u,v MAB_SPONGE_WIDTH=$SPONGE_W MAB_N_MATCH=$NMATCH MAB_CHECKPOINT_EVERY=$CKPT \
       MAB_MOORINGS=pioneer MAB_MASK_CLEANUP=true MAB_RIVERS=${RIVERS:-true} MAB_SSS_PISTON=0.5 MAB_SSS_RIVER_MASK=30,100 \
       MAB_TAG=/t0/workdir/enrique/runs/bgc/$RUN_TAG MAB_BOTTOM=partial MAB_ZSTAR=true MAB_TRACER_BOUNDS=T:-2:31,S:0:37.2 \
       MAB_STAGE=${STAGE:-run} MAB_PCO2=${PCO2:-411} MAB_BGC_DATA=${BGC_DATA:-woa} MAB_BGC_BOUNDS=${BGC_BOUNDS:-} MAB_PICKUP=${PICKUP:-false} \
       MAB_BGC_NO3_SOURCE=${NO3_SOURCE:-} MAB_BGC_SCALE_NEGATIVES=${SCALE_NEGATIVES:-false} MAB_RIVER_LOADS=${RIVER_LOADS:-true} MAB_RIVER_CHEMISTRY=${RIVER_CHEMISTRY:-rc4uscoast} MAB_RIVER_NO3=${RIVER_NO3:-Hudson:50,Delaware:50,other:20}
mpirun -np $SLURM_NNODES --map-by ppr:1:node --bind-to none -x LD_LIBRARY_PATH -x PATH -x JULIA_DEPOT_PATH -x JULIA_PKG_GC_AUTO \
    -x MPITRAMPOLINE_LIB -x MPITRAMPOLINE_MPIEXEC \
    -x MAB_START_DATE -x MAB_DAYS -x MAB_PARTITION_X -x MAB_PARTITION_Y -x MAB_NZ -x MAB_CELLS_PER_DEGREE -x MAB_DT -x MAB_PROGRESS_MINUTES \
    -x MAB_SPONGE_VARS -x MAB_SPONGE_WIDTH -x MAB_N_MATCH -x MAB_CHECKPOINT_EVERY -x MAB_MOORINGS -x MAB_MASK_CLEANUP -x MAB_RIVERS \
    -x MAB_SSS_PISTON -x MAB_SSS_RIVER_MASK -x MAB_TAG -x MAB_BOTTOM -x MAB_ZSTAR -x MAB_TRACER_BOUNDS -x MAB_STAGE -x MAB_PCO2 -x MAB_BGC_DATA -x MAB_BGC_BOUNDS -x MAB_PICKUP -x MAB_BGC_NO3_SOURCE -x MAB_BGC_SCALE_NEGATIVES -x MAB_RIVER_LOADS -x MAB_RIVER_CHEMISTRY -x MAB_RIVER_NO3 \
    julia -t $THREADS --project=$JENV bgc/mab_bgc.jl
