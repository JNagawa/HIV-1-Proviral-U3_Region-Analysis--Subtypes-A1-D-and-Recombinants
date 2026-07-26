#!/bin/bash
# Shared helpers for scripts/{download_qc,assembly,msa,biological_filtering,subtyping,motif_mapping}/compare_*.sh orchestrators.
# Source this from each stage's compare script:
#   source "$(dirname "$0")/../common/lib_compare.sh"

# Wraps a command with GNU time -v (wall-clock + peak RSS) when available,
# falling back to a wall-clock-only measurement (coarser, no RSS) otherwise.
#
# The fallback does NOT use `{ time "$@"; } 2>"${timelog}"`: that idiom
# redirects the wrapped command's own stderr into the timelog too (it shares
# fd 2 with `time`'s own report inside the group), silently swallowing every
# tool's real error output into the .time file instead of the caller's
# .log file. Timing via `date` instead leaves the command's stdout/stderr
# completely untouched.
# Usage: measure_and_run <timelog_path> -- <command...>
measure_and_run() {
    local timelog="$1"                               # arg 1 = file to write timing metrics to
    shift                                             # drop it so "$@" is now the -- and command
    if [ "$1" = "--" ]; then shift; fi               # allow an optional -- separator before the command
    if command -v /usr/bin/time >/dev/null 2>&1; then  # prefer GNU time (gives wall-clock AND peak RSS)
        /usr/bin/time -v -o "${timelog}" "$@"        # run the command under GNU time, metrics -> timelog
        return $?                                    # propagate the command's own exit status
    else
        local start end rc                           # fallback: only wall-clock via date (no RSS)
        start=$(date +%s.%N)                         # record start time with nanosecond precision
        "$@"                                         # run the command untouched (stdout/stderr untouched)
        rc=$?                                         # capture its exit status before overwriting $?
        end=$(date +%s.%N)                           # record end time
        awk -v s="${start}" -v e="${end}" 'BEGIN{d=e-s; printf "real\t%dm%.3fs\n", int(d/60), d-int(d/60)*60}' > "${timelog}"  # write a "real 0m0.000s" line parse_time_metrics can read
        return ${rc}                                 # propagate the command's exit status, not awk's
    fi
}

# Extracts wall-clock seconds and peak RSS (MB) from a timing file into
# WALLCLOCK_SEC and PEAK_RSS_MB. Handles GNU `time -v` output (both metrics)
# and falls back to parsing the bash-builtin `time` format (wall-clock only
# -- RSS isn't available from the builtin) when /usr/bin/time isn't present.
parse_time_metrics() {
    local f="$1"                                     # arg 1 = the .time file to parse
    WALLCLOCK_SEC=""                                 # reset outputs so a failed parse yields empty, not stale
    PEAK_RSS_MB=""                                   # reset the RSS global too
    if grep -q "Elapsed (wall clock)" "${f}" 2>/dev/null; then  # this line only appears in GNU time -v output
        local raw                                    # will hold the raw "[h:]mm:ss[.ss]" string
        raw=$(grep "Elapsed (wall clock)" "${f}" | awk -F': ' '{print $2}' | tail -1)  # pull the time value after the colon
        # format is [h:]mm:ss[.ss] -- convert to seconds
        WALLCLOCK_SEC=$(echo "${raw}" | awk -F: '{
            if (NF==3) print $1*3600+$2*60+$3;
            else if (NF==2) print $1*60+$2;
            else print $1
        }')                                          # handle h:m:s, m:s, or bare-seconds forms -> total seconds
        local rss_kb                                 # will hold peak RSS in kilobytes
        rss_kb=$(grep "Maximum resident set size" "${f}" | awk -F': ' '{print $2}' | tail -1)  # GNU time reports peak RSS in KB
        if [ -n "${rss_kb}" ]; then                  # only convert if the field was present
            PEAK_RSS_MB=$(awk -v kb="${rss_kb}" 'BEGIN{printf "%.1f", kb/1024}')  # KB -> MB for the summary table
        fi
    elif grep -qE "^real[[:space:]]" "${f}" 2>/dev/null; then  # else fall back to the bash-builtin "real  0m1.897s" line
        # bash builtin `time` fallback, e.g. "real  0m1.897s" (no RSS available)
        local raw                                    # will hold the "0m1.897s" string
        raw=$(grep -E "^real[[:space:]]" "${f}" | tail -1 | awk '{print $2}')  # grab the value after "real"
        WALLCLOCK_SEC=$(echo "${raw}" | sed -E 's/([0-9]+)m([0-9.]+)s/\1 \2/' | awk '{printf "%.2f", $1*60+$2}')  # split "Nm S.SSSs" and convert minutes+seconds -> total seconds
    fi
}

# Appends one row to the stage's summary.tsv. Columns:
# stage  tool  sample  wallclock_sec  peak_rss_mb  exit_code  output_valid  key_metric  notes_placeholder
append_summary_row() {
    local stage="$1" tool="$2" sample="$3" wallclock="$4" rss="$5" exit_code="$6" valid="$7" metric="$8"  # positional args -> named locals for readability
    if [ ! -f "${SUMMARY_TSV}" ]; then               # first call creates the file and its header row
        printf "stage\ttool\tsample\twallclock_sec\tpeak_rss_mb\texit_code\toutput_valid\tkey_metric\tnotes_placeholder\n" > "${SUMMARY_TSV}"  # write the tab-separated header once
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t(see ease_of_use_notes.md)\n" \
        "${stage}" "${tool}" "${sample}" "${wallclock}" "${rss}" "${exit_code}" "${valid}" "${metric}" \
        >> "${SUMMARY_TSV}"                           # append this run's row; notes column points to the manual notes file
}

# Reads scripts/common/subset_samples.tsv and prints accessions for a given
# platform ("illumina" or "nanopore"), one per line.
subset_accessions() {
    local platform="$1" tsv="$2"                     # arg 1 = platform filter, arg 2 = the subset TSV path
    awk -F'\t' -v p="${platform}" '$1==p {print $2}' "${tsv}"  # print col2 (accession) of rows whose col1 matches the platform
}
