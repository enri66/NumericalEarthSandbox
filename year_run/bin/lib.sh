# Shared helpers (sourced). Needs config.env sourced first.

# Iterations that have a checkpoint file for every rank, ascending.
complete_iterations() {
    ( cd $OUT && ls mabyear_s*_checkpoint_rank*_iteration*.jld2 2>/dev/null ) | sed -E 's/.*_rank([0-9]+)_iteration([0-9]+)\.jld2/\1 \2/' | sort -u |
    awk -v n=$NRANKS '{c[$2]++} END {for (i in c) if (c[i] == n) print i}' | sort -n
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
