#!/bin/bash
# Compares Poplars (Hypermut 3) vs HIVSeqinR vs HIVIntact on the assembled
# subset (HXB2 + per-sample consensus sequences from the assembly stage). No gold
# standard exists for this cohort, so "validity" here is just "did it run
# and produce a classification" -- cross-tool agreement is recorded as a
# qualitative note, not an accuracy score.
# Usage: ./compare_biological_filtering_oxfordnano.sh
# -u errors on unset vars, pipefail fails a pipe if any stage fails; no -e so one tool failing
# doesn't kill the comparison
set -uo pipefail

# absolute path of this script's own dir, so paths work regardless of launch location
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
# repo top-level (three dirs up), the base for every other path below
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
# load shared helpers: measure_and_run, parse_time_metrics, append_summary_row
source "${REPO_ROOT}/scripts/common/lib_compare.sh"

# HXB2 reference genome, used as the first record (reference) for the tools
REF_FASTA="${REPO_ROOT}/data/reference/K03455.1.fasta"
# input: per-sample consensus FASTAs from the ONT assembly stage
ASSEMBLY_DIR="${REPO_ROOT}/results/assembly/oxfordnano/minimap2_out"
# output: all classification results + summary go here
RESULTS_DIR="${REPO_ROOT}/results/biological_filtering/oxfordnano"
# the single TSV every tool appends a timing/validity row to
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
# create the results dir (and parents) if it doesn't exist
mkdir -p "${RESULTS_DIR}"
# seed the manual notes file from the template only on first run
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

# single FASTA fed to all three tools (HXB2 first, then consensuses)
INPUT_FASTA="${RESULTS_DIR}/combined_input.fasta"
# concatenate reference + all consensus seqs; HXB2 first so Poplars uses it as reference
cat "${REF_FASTA}" "${ASSEMBLY_DIR}"/*_consensus.fasta > "${INPUT_FASTA}" 2>/dev/null
# count FASTA records to check we have enough input
N_SEQS=$(grep -c "^>" "${INPUT_FASTA}" 2>/dev/null || echo 0)
# need at least the reference plus one consensus to compare
if [ "${N_SEQS}" -lt 2 ]; then
    # tell the user what to run first
    echo "ERROR: fewer than 2 sequences available (run scripts/assembly/oxfordnano/assembly_oxfordnano.sh (minimap2) first). Found ${N_SEQS}." >&2
    exit 1                                           # bail; nothing to classify
fi

if [ ! -d "${REPO_ROOT}/scripts/tools/Poplars" ]; then  # the tool clones are a prerequisite
    # point the user at the setup script
    echo "NOTE: tools not set up yet -- run ./setup_tools.sh first." >&2
    exit 1                                           # bail; tools not installed
fi

echo "=== Poplars (Hypermut 3) ==="                  # progress marker in the log
OUT="${RESULTS_DIR}/poplars_out.tsv"                 # Poplars hypermutation results table
# file where measure_and_run records wallclock/RSS
TIMELOG="${RESULTS_DIR}/poplars.time"
# run Poplars wrapper under timing, capturing stdout+stderr
measure_and_run "${TIMELOG}" -- "${STAGE_DIR}/run_poplars.sh" "${INPUT_FASTA}" "${OUT}" > "${RESULTS_DIR}/poplars.log" 2>&1
# capture the wrapper's exit status before $? is overwritten
EXIT_CODE=$?
# set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file
parse_time_metrics "${TIMELOG}"
VALID=0; METRIC="n/a"                                 # assume invalid until proven otherwise
# valid if output non-empty; record row count as the key metric
if [ -s "${OUT}" ]; then VALID=1; METRIC="$(wc -l < "${OUT}") result rows"; fi
# write Poplars' row to summary.tsv
append_summary_row "biological_filtering_oxfordnano" "poplars" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"

echo "=== HIVSeqinR ==="                              # progress marker in the log
OUTDIR="${RESULTS_DIR}/hivseqinr_out"                 # per-tool output dir for HIVSeqinR
TIMELOG="${RESULTS_DIR}/hivseqinr.time"               # timing file for this run
# run HIVSeqinR wrapper under timing, capturing output
measure_and_run "${TIMELOG}" -- "${STAGE_DIR}/run_hivseqinr.sh" "${INPUT_FASTA}" "${OUTDIR}" > "${RESULTS_DIR}/hivseqinr_wrapper.log" 2>&1
EXIT_CODE=$?                                          # capture the wrapper's exit status
# refresh WALLCLOCK_SEC / PEAK_RSS_MB from the timing file
parse_time_metrics "${TIMELOG}"
VALID=0; METRIC="n/a"                                 # default to invalid until checked
CSV="${OUTDIR}/Output_MyBigSummary_DF_FINAL.csv"      # HIVSeqinR's final summary table
# valid if CSV exists; subtract 1 for the header row to get sequences classified
if [ -s "${CSV}" ]; then VALID=1; METRIC="$(($(wc -l < "${CSV}") - 1)) classified"; fi
# write HIVSeqinR's row to summary.tsv
append_summary_row "biological_filtering_oxfordnano" "hivseqinr" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"

echo "=== HIVIntact ==="                              # progress marker in the log
OUTDIR="${RESULTS_DIR}/hivintact_out"                 # per-tool output dir for HIVIntact
TIMELOG="${RESULTS_DIR}/hivintact.time"               # timing file for this run
# run HIVIntact wrapper with subtype B
measure_and_run "${TIMELOG}" -- "${STAGE_DIR}/run_hivintact.sh" "${INPUT_FASTA}" "${OUTDIR}" B > "${RESULTS_DIR}/hivintact_wrapper.log" 2>&1
EXIT_CODE=$?                                          # capture the wrapper's exit status
# refresh WALLCLOCK_SEC / PEAK_RSS_MB from the timing file
parse_time_metrics "${TIMELOG}"
VALID=0; METRIC="n/a"                                 # default to invalid until checked
# valid if either the intact or non-intact output was produced
if [ -s "${OUTDIR}/intact.fasta" ] || [ -s "${OUTDIR}/nonintact.fasta" ]; then
    VALID=1                                           # mark this run as valid
    # count intact proviruses (0 if file missing)
    N_INTACT=$(grep -c "^>" "${OUTDIR}/intact.fasta" 2>/dev/null || echo 0)
    # count non-intact proviruses (0 if file missing)
    N_NONINTACT=$(grep -c "^>" "${OUTDIR}/nonintact.fasta" 2>/dev/null || echo 0)
    # record the intact/non-intact split as the key metric
    METRIC="${N_INTACT} intact, ${N_NONINTACT} non-intact"
fi
# write HIVIntact's row to summary.tsv
append_summary_row "biological_filtering_oxfordnano" "hivintact" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"

# final confirmation pointing the user at the results table
echo "Done. See ${SUMMARY_TSV}"
