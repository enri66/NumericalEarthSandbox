#!/bin/bash
# Downloads the global GLORYS12 daily means (Copernicus Marine GLOBAL_MULTIYEAR_PHY_001_030,
# cmems_mod_glo_phy_my_0.083deg_P1D-m, 1993 to the present) as Mercator's original daily files, about 1.3 GB each
# with every variable (thetao, so, uo, vo, zos, mlotst, bottomT and the sea-ice fields), into <root>/<YYYY>/<MM>/.
# One month per `copernicusmarine get`; files already present are skipped, so it can be rerun at any time.
# Usage (on antares):  bash scripts/download_glorys_global.sh [root] [first YYYY-MM] [last YYYY-MM]
# The toolbox lives in ~/cmenv (pip install copernicusmarine==2.3.0; antares' Python 3.9 cannot take 2.4).
ROOT=${1:-/Volumes/A3/enrique/GLORYS12}
FIRST=${2:-1993-01}
LAST=${3:-$(date +%Y-%m)}
CM=$HOME/cmenv/bin/copernicusmarine
DATASET=cmems_mod_glo_phy_my_0.083deg_P1D-m

month=$FIRST
failed=""
while [[ ! "$month" > "$LAST" ]]; do
    y=${month%-*}; m=${month#*-}
    mkdir -p "$ROOT/$y/$m"
    t0=$(date +%s)
    if $CM get --dataset-id $DATASET --filter "*_${y}${m}??_*" --output-directory "$ROOT/$y/$m" --no-directories \
               --skip-existing --disable-progress-bar --log-level WARN > "$ROOT/$y/$m/.get.log" 2>&1; then
        n=$(ls "$ROOT/$y/$m"/*.nc 2>/dev/null | wc -l)
        echo "$(date '+%Y-%m-%d %H:%M')  $month  $n files  $(du -sh "$ROOT/$y/$m" | cut -f1)  $(( ($(date +%s) - t0) / 60 )) min"
    elif grep -q "No data to download" "$ROOT/$y/$m/.get.log"; then
        echo "$(date '+%Y-%m-%d %H:%M')  $month  not published yet"
        rm -f "$ROOT/$y/$m/.get.log"; rmdir "$ROOT/$y/$m" 2>/dev/null
    else
        echo "$(date '+%Y-%m-%d %H:%M')  $month  FAILED (see $ROOT/$y/$m/.get.log)"
        failed="$failed $month"
    fi
    month=$(date -d "$month-01 +1 month" +%Y-%m)
done
[ -n "$failed" ] && echo "failed, rerun to retry:$failed"
