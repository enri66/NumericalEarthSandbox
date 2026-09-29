#!/bin/bash
# Wait for download_era5_year.jl to finish, sync ERA5 to triton, verify every variable reaches 2020-04-01T01,
# then release the held phase-2 job of the sponge-year run (triton job 11236).
LOG=~/Models/NumericalEarthSandbox/era5_year_download.log
while pgrep -f download_era5_year >/dev/null; do sleep 120; done
echo "$(date) download process finished"
grep -q "ERA5 complete" $LOG || { echo "$(date) download did not complete - not syncing, phase 2 stays held"; exit 1; }
rsync -a ~/Data/NumericalEarth/era5/ triton:Data/NumericalEarth/era5/ || { echo "$(date) rsync failed"; exit 1; }
echo "$(date) rsync done"
ssh -o BatchMode=yes triton bash -s <<"REMOTE"
d=~/Data/NumericalEarth/era5
ok=1
for v in 10m_u_component_of_wind 10m_v_component_of_wind 2m_dewpoint_temperature 2m_temperature surface_pressure \
         surface_solar_radiation_downwards surface_thermal_radiation_downwards total_precipitation; do
    n=$(ls $d | grep -c "^${v}_ERA5HourlySingleLevel_")
    last=$(ls $d | grep "^${v}_ERA5HourlySingleLevel_" | grep -oE "20[0-9]{2}-[0-9]{2}-[0-9]{2}T[0-9]{2}" | sort | tail -1)
    echo "$v: $n files, last $last"
    { [ "$n" -ge 8786 ] && [ "$last" = "2020-04-01T01" ]; } || ok=0
done
if [ $ok = 1 ]; then scontrol release 11236 && echo "RELEASED phase 2 (11236)"; else echo "INCOMPLETE - phase 2 stays held"; fi
squeue -h
REMOTE
