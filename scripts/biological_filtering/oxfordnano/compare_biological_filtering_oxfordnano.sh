#!/bin/bash
# Compares Poplars (Hypermut 3) vs HIVSeqinR vs HIVIntact on the assembled
# subset (HXB2 + per-sample consensus sequences from the assembly stage). No gold
# standard exists for this cohort, so "validity" here is just "did it run
# and produce a classification" -- cross-tool agreement is recorded as a
# qualitative note, not an accuracy score.
# Usage: ./compare_biological_filtering_oxfordnano.sh
set -uo pipefail                                     # -u errors on unset vars, pipefail fails a pipe if any stage fails; no -e so one tool failing doesn't kill the comparison

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"           # absolute path of this script's own dir, so paths work regardless of launch location
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"     # repo top-level (three dirs up), the base for every other path below
source "${REPO_ROOT}/scripts/common/lib_compare.sh"  # load shared helpers: measure_and_run, parse_time_metrics, append_summary_row

REF_FASTA="${REPO_ROOT}/data/reference/K03455.1.fasta"                 # HXB2 reference genome, used as the first record (reference) for the tools
ASSEMBLY_DIR="${REPO_ROOT}/results/assembly/oxfordnano/minimap2_out"   # input: per-sample consensus FASTAs from the ONT assembly stage
RESULTS_DIR="${REPO_ROOT}/results/biological_filtering/oxfordnano"     # output: all classification results + summary go here
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"                               # the single TSV every tool appends a timing/validity row to
mkdir -p "${RESULTS_DIR}"                            # create the results dir (and parents) if it doesn't exist
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"  # seed the manual notes file from the template only on first run

INPUT_FASTA="${RESULTS_DIR}/combined_input.fasta"    # single FASTA fed to all three tools (HXB2 first, then consensuses)
cat "${REF_FASTA}" "${ASSEMBLY_DIR}"/*_consensus.fasta > "${INPUT_FASTA}" 2>/dev/null  # concatenate reference + all consensus seqs; HXB2 first so Poplars uses it as reference
N_SEQS=$(grep -c "^>" "${INPUT_FASTA}" 2>/dev/null || echo 0)  # count FASTA records to check we have enough input
if [ "${N_SEQS}" -lt 2 ]; then                       # need at least the reference plus one consensus to compare
    echo "ERROR: fewer than 2 sequences available (run scripts/assembly/oxfordnano/assembly_oxfordnano.sh (minimap2) first). Found ${N_SEQS}." >&2  # tell the user what to run first
    exit 1                                           # bail; nothing to classify
fi

if [ ! -d "${REPO_ROOT}/scripts/tools/Poplars" ]; then  # the tool clones are a prerequisite
    echo "NOTE: tools not set up yet -- run ./setup_tools.sh first." >&2  # point the user at the setup script
    exit 1                                           # bail; tools not installed
fi

echo "=== Poplars (Hypermut 3) ==="                  # progress marker in the log
OUT="${RESULTS_DIR}/poplars_out.tsv"                 # Poplars hypermutation results table
TIMELOG="${RESULTS_DIR}/poplars.time"                # file where measure_and_run records wallclock/RSS
measure_and_run "${TIMELOG}" -- "${STAGE_DIR}/run_poplars.sh" "${INPUT_FASTA}" "${OUT}" > "${RESULTS_DIR}/poplars.log" 2>&1  # run Poplars wrapper under timing, capturing stdout+stderr
EXIT_CODE=$?                                          # capture the wrapper's exit status before $? is overwritten
parse_time_metrics "${TIMELOG}"                       # set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file
VALID=0; METRIC="n/a"                                 # assume invalid until proven otherwise
if [ -s "${OUT}" ]; then VALID=1; METRIC="$(wc -l < "${OUT}") result rows"; fi  # valid if output non-empty; record row count as the key metric
append_summary_row "biological_filtering_oxfordnano" "poplars" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"  # write Poplars' row to summary.tsv

echo "=== HIVSeqinR ==="                              # progress marker in the log
OUTDIR="${RESULTS_DIR}/hivseqinr_out"                 # per-tool output dir for HIVSeqinR
TIMELOG="${RESULTS_DIR}/hivseqinr.time"               # timing file for this run
measure_and_run "${TIMELOG}" -- "${STAGE_DIR}/run_hivseqinr.sh" "${INPUT_FASTA}" "${OUTDIR}" > "${RESULTS_DIR}/hivseqinr_wrapper.log" 2>&1  # run HIVSeqinR wrapper under timing, capturing output
EXIT_CODE=$?                                          # capture the wrapper's exit status
parse_time_metrics "${TIMELOG}"                       # refresh WALLCLOCK_SEC / PEAK_RSS_MB from the timing file
VALID=0; METRIC="n/a"                                 # default to invalid until checked
CSV="${OUTDIR}/Output_MyBigSummary_DF_FINAL.csv"      # HIVSeqinR's final summary table
if [ -s "${CSV}" ]; then VALID=1; METRIC="$(($(wc -l < "${CSV}") - 1)) classified"; fi  # valid if CSV exists; subtract 1 for the header row to get sequences classified
append_summary_row "biological_filtering_oxfordnano" "hivseqinr" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"  # write HIVSeqinR's row to summary.tsv

echo "=== HIVIntact ==="                              # progress marker in the log
OUTDIR="${RESULTS_DIR}/hivintact_out"                 # per-tool output dir for HIVIntact
TIMELOG="${RESULTS_DIR}/hivintact.time"               # timing file for this run
measure_and_run "${TIMELOG}" -- "${STAGE_DIR}/run_hivintact.sh" "${INPUT_FASTA}" "${OUTDIR}" B > "${RESULTS_DIR}/hivintact_wrapper.log" 2>&1  # run HIVIntact wrapper with subtype B
EXIT_CODE=$?                                          # capture the wrapper's exit status
parse_time_metrics "${TIMELOG}"                       # refresh WALLCLOCK_SEC / PEAK_RSS_MB from the timing file
VALID=0; METRIC="n/a"                                 # default to invalid until checked
if [ -s "${OUTDIR}/intact.fasta" ] || [ -s "${OUTDIR}/nonintact.fasta" ]; then  # valid if either the intact or non-intact output was produced
    VALID=1                                           # mark this run as valid
    N_INTACT=$(grep -c "^>" "${OUTDIR}/intact.fasta" 2>/dev/null || echo 0)        # count intact proviruses (0 if file missing)
    N_NONINTACT=$(grep -c "^>" "${OUTDIR}/nonintact.fasta" 2>/dev/null || echo 0)  # count non-intact proviruses (0 if file missing)
    METRIC="${N_INTACT} intact, ${N_NONINTACT} non-intact"  # record the intact/non-intact split as the key metric
fi
append_summary_row "biological_filtering_oxfordnano" "hivintact" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"  # write HIVIntact's row to summary.tsv

echo "Done. See ${SUMMARY_TSV}"                       # final confirmation pointing the user at the results table
