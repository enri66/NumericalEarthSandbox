# Shared helpers (sourced). Needs config.env sourced first.

# Iterations that have a checkpoint file for every rank, none of them cut short, ascending. A file is cut short when it is
# smaller than 90% of the largest checkpoint file of the same rank (a power loss during the write left 50% on 2026-10-26).
complete_iterations() {
    ( cd $OUT && for f in mabyear_s*_checkpoint_rank*_iteration*.jld2; do [ -e "$f" ] && stat -L -c '%s %n' "$f"; done ) 2>/dev/null |
    sed -E 's/^([0-9]+) .*_rank([0-9]+)_iteration([0-9]+)\.jld2$/\2 \3 \1/' | sort -u |
    awk -v n=$NRANKS '{ rank=$1; it=$2; size=$3; sz[rank" "it]=size; if (size > mx[rank]) mx[rank]=size; its[it]=1 }
        END { for (it in its) { ok=0; for (r in mx) { k=r" "it; if (k in sz && sz[k] >= 0.9*mx[r]) ok++ } if (ok == n) print it } }' | sort -n
}

# Simulated iteration of the last progress line in a log (empty if none).
last_progress_iteration() {
    local line
    line=$(grep -a -E '^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}  \|u\|' "$1" 2>/dev/null | tail -1 | cut -c1-16)
    [ -z "$line" ] && return 0
    local t0 t1
    t0=$(date -u -d "$START_DATE" +%s); t1=$(date -u -d "$line:00" +%s)
    echo $(( (t1 - t0) / DT ))
}
