#!/bin/bash
# Satellite SSH (DUACS L4 DT, 1/8 deg: adt, sla, ugos, vgos) and SST (OSTIA reprocessed, 1/20 deg) from Copernicus Marine for
# the 2019-04 year run, over a box wider than the model (80-50W, 30-45N) so rings east of the model can be counted too.
# The toolbox needs h5py to write NetCDF: ~/cmenv_subset has it, ~/cmenv does not (and `subset` then exits 0 without a file).
#   download_satellite_ssh_sst.sh [output directory, default /Volumes/A3/enrique/Satellite]
CM=$HOME/cmenv_subset/bin/copernicusmarine
OUT=${1:-/Volumes/A3/enrique/Satellite}
mkdir -p $OUT
B="--minimum-longitude -80 --maximum-longitude -50 --minimum-latitude 30 --maximum-latitude 45 --start-datetime 2019-03-28 --end-datetime 2020-04-04"
echo "=== $(date) start"
$CM subset --dataset-id cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.125deg_P1D -v adt -v sla -v ugos -v vgos $B \
    --service arco-time-series --output-directory $OUT --output-filename duacs_l4_0125_ssh_2019-03-28_2020-04-04.nc \
    --disable-progress-bar --log-level WARN
$CM subset --dataset-id METOFFICE-GLO-SST-L4-REP-OBS-SST -v analysed_sst -v mask $B \
    --service arco-time-series --output-directory $OUT --output-filename ostia_rep_sst_2019-03-28_2020-04-04.nc \
    --disable-progress-bar --log-level WARN
echo "=== $(date) done"
ls -l $OUT
