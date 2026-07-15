#!/bin/bash
# Shared helpers for scripts/{download_qc,assembly,msa,biological_filtering,subtyping,motif_mapping}/compare_*.sh orchestrators.
# Source this from each stage's compare script:
#   source "$(dirname "$0")/../common/lib_compare.sh"

# Wraps a command with GNU time -v (wall-clock + peak RSS) when available,
# falling back to the bash builtin `time` (coarser, no RSS) otherwise.
# Usage: measure_and_run <timelog_path> -- <command...>
measure_and_run() {
    local timelog="$1"
    shift
    if [ "$1" = "--" ]; then shift; fi
    if command -v /usr/bin/time >/dev/null 2>&1; then
        /usr/bin/time -v -o "${timelog}" "$@"
    else
        { time "$@"; } 2>"${timelog}"
    fi
    return $?
}

# Extracts wall-clock seconds and peak RSS (MB) from a timing file into
# WALLCLOCK_SEC and PEAK_RSS_MB. Handles GNU `time -v` output (both metrics)
# and falls back to parsing the bash-builtin `time` format (wall-clock only
# -- RSS isn't available from the builtin) when /usr/bin/time isn't present.
parse_time_metrics() {
    local f="$1"
    WALLCLOCK_SEC=""
    PEAK_RSS_MB=""
    if grep -q "Elapsed (wall clock)" "${f}" 2>/dev/null; then
        local raw
        raw=$(grep "Elapsed (wall clock)" "${f}" | awk -F': ' '{print $2}' | tail -1)
        # format is [h:]mm:ss[.ss] -- convert to seconds
        WALLCLOCK_SEC=$(echo "${raw}" | awk -F: '{
            if (NF==3) print $1*3600+$2*60+$3;
            else if (NF==2) print $1*60+$2;
            else print $1
        }')
        local rss_kb
        rss_kb=$(grep "Maximum resident set size" "${f}" | awk -F': ' '{print $2}' | tail -1)
        if [ -n "${rss_kb}" ]; then
            PEAK_RSS_MB=$(awk -v kb="${rss_kb}" 'BEGIN{printf "%.1f", kb/1024}')
        fi
    elif grep -qE "^real[[:space:]]" "${f}" 2>/dev/null; then
        # bash builtin `time` fallback, e.g. "real  0m1.897s" (no RSS available)
        local raw
        raw=$(grep -E "^real[[:space:]]" "${f}" | tail -1 | awk '{print $2}')
        WALLCLOCK_SEC=$(echo "${raw}" | sed -E 's/([0-9]+)m([0-9.]+)s/\1 \2/' | awk '{printf "%.2f", $1*60+$2}')
    fi
}

# Appends one row to the stage's summary.tsv. Columns:
# stage  tool  sample  wallclock_sec  peak_rss_mb  exit_code  output_valid  key_metric  notes_placeholder
append_summary_row() {
    local stage="$1" tool="$2" sample="$3" wallclock="$4" rss="$5" exit_code="$6" valid="$7" metric="$8"
    if [ ! -f "${SUMMARY_TSV}" ]; then
        printf "stage\ttool\tsample\twallclock_sec\tpeak_rss_mb\texit_code\toutput_valid\tkey_metric\tnotes_placeholder\n" > "${SUMMARY_TSV}"
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t(see ease_of_use_notes.md)\n" \
        "${stage}" "${tool}" "${sample}" "${wallclock}" "${rss}" "${exit_code}" "${valid}" "${metric}" \
        >> "${SUMMARY_TSV}"
}

# Reads scripts/common/subset_samples.tsv and prints accessions for a given
# platform ("illumina" or "nanopore"), one per line.
subset_accessions() {
    local platform="$1" tsv="$2"
    awk -F'\t' -v p="${platform}" '$1==p {print $2}' "${tsv}"
}
