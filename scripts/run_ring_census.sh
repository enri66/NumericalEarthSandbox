#!/bin/bash
# Gulf Stream ring census of the year run (ring_census.jl) on antares, from the synced copy of its output. Any of the settings
# below can be overridden from the environment, e.g. for a 1/12 degree run:
#   MAB_TAGS=<prefix> MAB_CELLS_PER_DEGREE=12 MAB_PARTITION_X=1 MAB_PARTITION_Y=1 RING_MODEL_NAME=<case> RING_OUT=<dir> run_ring_census.sh
D=/Volumes/A3/enrique/MAB_year_2019-04
export MAB_TAGS=${MAB_TAGS:-$D/mabyear_s01,$D/mabyear_s02,$D/mabyear_s03}
export MAB_START_DATE=${MAB_START_DATE:-2019-04-01}
export SAT_SSH=${SAT_SSH:-/Volumes/A3/enrique/Satellite/duacs_l4_0125_ssh_2019-03-28_2020-04-04.nc}
export GLORYS_ROOT=${GLORYS_ROOT:-/Volumes/A3/enrique/GLORYS12}
export RING_OUT=${RING_OUT:-$HOME/Data/mab_analysis/rings_2019}
export RING_VIDEO=${RING_VIDEO:-false}
cd $(dirname $(readlink -f $0))/..
~/.juliaup/bin/julia -t 8 --project=. scripts/ring_census.jl
