#!/bin/bash
# The observations bgc/validate_bgc.jl compares a BGC run with, for 2019-01-25 → 2019-04-02 (T0=2019-04-01 T1=2019-06-02 for
# another period; the satellite command at the end has its dates to change too). Run where there is internet (the Mac, or antares), then copy the files to triton:
#   OOI Pioneer BGC sensors  → /t0/workdir/enrique/Data/OOI/pioneer_bgc_2019/   (CSV from OOI's ERDDAP)
#   Copernicus ocean colour  → /t0/workdir/enrique/Data/ocean_colour/            (L3 daily multi-sensor chlorophyll, 4 km)
# Notes from using them in 2019:
#   - the nitrate files must hold the `_suna` column (the main and `_full` nitrate columns are empty in 2019);
#   - the offshore mooring (cp04ossm) has data only until 2019-02-01 in this window;
#   - the near-surface oxygen optodes (dosta_nsif) read 123-131 % saturation at CNSM and ISSM in Feb-Mar 2019, with QC flag 1.
set -e
B="https://erddap.dataexplorer.oceanobservatories.org/erddap/tabledap"
T0=${T0:-2019-01-25}; T1=${T1:-2019-04-02}      # the period (environment T0, T1)
T="&time>=${T0}T00:00:00Z&time<=${T1}T00:00:00Z"
OUT=${1:-ooi_bgc_2019}; mkdir -p "$OUT"; cd "$OUT"
get() { curl -s -m 300 -o "$1.csv" "$B/$2.csv?$3$T"; }
for m in cp01cnsm cp03issm; do      # central and inshore surface moorings (the offshore one is off station)
  get ${m}_flort ooi-$m-rid27-02-flortd000 "time,latitude,longitude,z,mass_concentration_of_chlorophyll_a_in_sea_water,mass_concentration_of_chlorophyll_a_in_sea_water_qc_agg"
  for node in nsif:rid27-04 mfn:mfd37-04; do
    get ${m}_dosta_${node%%:*} ooi-$m-${node##*:}-dostad000 "time,latitude,longitude,z,mole_concentration_of_dissolved_molecular_oxygen_in_sea_water,mole_concentration_of_dissolved_molecular_oxygen_in_sea_water_qc_agg,sea_water_temperature,sea_water_practical_salinity"
  done
  get ${m}_nutnr_all ooi-$m-rid26-07-nutnrb000 "time,latitude,longitude,z,mole_concentration_of_nitrate_in_sea_water_suna,mole_concentration_of_nitrate_in_sea_water_suna_qc_agg"
done

# Satellite chlorophyll. The toolbox needs h5py to write NetCDF: a venv with `pip install copernicusmarine==2.3.0 h5py`
# (antares: ~/cmenv_subset; ~/cmenv is the one the global download watchdogs use and has no h5py).
#   source ~/cmenv_subset/bin/activate
#   copernicusmarine subset --dataset-id cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D --variable CHL --variable flags \
#       --minimum-longitude -76 --maximum-longitude -64 --minimum-latitude 34 --maximum-latitude 42 \
#       --start-datetime 2019-01-25T00:00:00 --end-datetime 2019-04-02T00:00:00 \
#       --output-filename chl_l3_multi_4km_mab_2019.nc --output-directory . --disable-progress-bar
