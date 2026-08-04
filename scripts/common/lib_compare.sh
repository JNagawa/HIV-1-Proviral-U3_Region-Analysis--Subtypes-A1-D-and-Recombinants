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
    # allow an optional -- separator before the command
    if [ "$1" = "--" ]; then shift; fi
    # prefer GNU time (gives wall-clock AND peak RSS)
    if command -v /usr/bin/time >/dev/null 2>&1; then
        # run the command under GNU time, metrics -> timelog
        /usr/bin/time -v -o "${timelog}" "$@"
        return $?                                    # propagate the command's own exit status
    else
        local start end rc                           # fallback: only wall-clock via date (no RSS)
        start=$(date +%s.%N)                         # record start time with nanosecond precision
        # run the command untouched (stdout/stderr untouched)
        "$@"
        # capture its exit status before overwriting $?
        rc=$?
        end=$(date +%s.%N)                           # record end time
        # write a "real 0m0.000s" line parse_time_metrics can read
        awk -v s="${start}" -v e="${end}" 'BEGIN{d=e-s; printf "real\t%dm%.3fs\n", int(d/60), d-int(d/60)*60}' > "${timelog}"
        # propagate the command's exit status, not awk's
        return ${rc}
    fi
}

# Extracts wall-clock seconds and peak RSS (MB) from a timing file into
# WALLCLOCK_SEC and PEAK_RSS_MB. Handles GNU `time -v` output (both metrics)
# and falls back to parsing the bash-builtin `time` format (wall-clock only
# -- RSS isn't available from the builtin) when /usr/bin/time isn't present.
parse_time_metrics() {
    local f="$1"                                     # arg 1 = the .time file to parse
    # reset outputs so a failed parse yields empty, not stale
    WALLCLOCK_SEC=""
    PEAK_RSS_MB=""                                   # reset the RSS global too
    # this line only appears in GNU time -v output
    if grep -q "Elapsed (wall clock)" "${f}" 2>/dev/null; then
        local raw                                    # will hold the raw "[h:]mm:ss[.ss]" string
        # pull the time value after the colon
        raw=$(grep "Elapsed (wall clock)" "${f}" | awk -F': ' '{print $2}' | tail -1)
        # format is [h:]mm:ss[.ss] -- convert to seconds
        WALLCLOCK_SEC=$(echo "${raw}" | awk -F: '{
            if (NF==3) print $1*3600+$2*60+$3;
            else if (NF==2) print $1*60+$2;
            else print $1
        }')                                          # handle h:m:s, m:s, or bare-seconds forms -> total seconds
        local rss_kb                                 # will hold peak RSS in kilobytes
        # GNU time reports peak RSS in KB
        rss_kb=$(grep "Maximum resident set size" "${f}" | awk -F': ' '{print $2}' | tail -1)
        if [ -n "${rss_kb}" ]; then                  # only convert if the field was present
            # KB -> MB for the summary table
            PEAK_RSS_MB=$(awk -v kb="${rss_kb}" 'BEGIN{printf "%.1f", kb/1024}')
        fi
    # else fall back to the bash-builtin "real 0m1.897s" line
    elif grep -qE "^real[[:space:]]" "${f}" 2>/dev/null; then
        # bash builtin `time` fallback, e.g. "real  0m1.897s" (no RSS available)
        local raw                                    # will hold the "0m1.897s" string
        # grab the value after "real"
        raw=$(grep -E "^real[[:space:]]" "${f}" | tail -1 | awk '{print $2}')
        # split "Nm S.SSSs" and convert minutes+seconds -> total seconds
        WALLCLOCK_SEC=$(echo "${raw}" | sed -E 's/([0-9]+)m([0-9.]+)s/\1 \2/' | awk '{printf "%.2f", $1*60+$2}')
    fi
}

# Appends one row to the stage's summary.tsv. Columns:
# stage  tool  sample  wallclock_sec  peak_rss_mb  exit_code  output_valid  key_metric  notes_placeholder
append_summary_row() {
    # positional args -> named locals for readability
    local stage="$1" tool="$2" sample="$3" wallclock="$4" rss="$5" exit_code="$6" valid="$7" metric="$8"
    # first call creates the file and its header row
    if [ ! -f "${SUMMARY_TSV}" ]; then
        # write the tab-separated header once
        printf "stage\ttool\tsample\twallclock_sec\tpeak_rss_mb\texit_code\toutput_valid\tkey_metric\tnotes_placeholder\n" > "${SUMMARY_TSV}"
    fi
    # re-running a stage should refresh its rows, not stack a second copy on top
    # of the stale ones, so drop any existing row for this exact
    # (stage, tool, sample) before appending the new one
    local tmp="${SUMMARY_TSV}.tmp"
    awk -F'\t' -v s="${stage}" -v t="${tool}" -v n="${sample}" \
        '!($1==s && $2==t && $3==n)' "${SUMMARY_TSV}" > "${tmp}" && mv "${tmp}" "${SUMMARY_TSV}"
    # append this run's row; notes column points to the manual notes file
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t(see ease_of_use_notes.md)\n" \
        "${stage}" "${tool}" "${sample}" "${wallclock}" "${rss}" "${exit_code}" "${valid}" "${metric}" \
        >> "${SUMMARY_TSV}"
}

# Reads scripts/common/subset_samples.tsv and prints accessions for a given
# platform ("illumina" or "nanopore"), one per line.
subset_accessions() {
    # arg 1 = platform filter, arg 2 = the subset TSV path
    local platform="$1" tsv="$2"
    # print col2 (accession) of rows whose col1 matches the platform
    awk -F'\t' -v p="${platform}" '$1==p {print $2}' "${tsv}"
}
