#!/bin/bash
# Compares Poplars (Hypermut 3) vs HIVSeqinR vs HIVIntact on the assembled
# subset (HXB2 + per-sample consensus sequences from stage 02). No gold
# standard exists for this cohort, so "validity" here is just "did it run
# and produce a classification" -- cross-tool agreement is recorded as a
# qualitative note, not an accuracy score.
# Usage: ./compare_04.sh
set -uo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
source "${REPO_ROOT}/scripts/tool_comparison/common/lib_compare.sh"

REF_FASTA="${REPO_ROOT}/data/reference/K03455.1.fasta"
ASSEMBLY_DIR="${REPO_ROOT}/results/tool_comparison/02_assembly_illumina/bwa_consensus_out"
RESULTS_DIR="${REPO_ROOT}/results/tool_comparison/04_biological_filtering"
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
mkdir -p "${RESULTS_DIR}"
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/tool_comparison/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

INPUT_FASTA="${RESULTS_DIR}/combined_input.fasta"
cat "${REF_FASTA}" "${ASSEMBLY_DIR}"/*_consensus.fasta > "${INPUT_FASTA}" 2>/dev/null
N_SEQS=$(grep -c "^>" "${INPUT_FASTA}" 2>/dev/null || echo 0)
if [ "${N_SEQS}" -lt 2 ]; then
    echo "ERROR: fewer than 2 sequences available (run stage 02's bwa_consensus first). Found ${N_SEQS}." >&2
    exit 1
fi

if [ ! -d "${STAGE_DIR}/tools/Poplars" ]; then
    echo "NOTE: tools not set up yet -- run ./setup_tools.sh first." >&2
    exit 1
fi

echo "=== Poplars (Hypermut 3) ==="
OUT="${RESULTS_DIR}/poplars_out.tsv"
TIMELOG="${RESULTS_DIR}/poplars.time"
measure_and_run "${TIMELOG}" -- "${STAGE_DIR}/run_poplars.sh" "${INPUT_FASTA}" "${OUT}" > "${RESULTS_DIR}/poplars.log" 2>&1
EXIT_CODE=$?
parse_time_metrics "${TIMELOG}"
VALID=0; METRIC="n/a"
if [ -s "${OUT}" ]; then VALID=1; METRIC="$(wc -l < "${OUT}") result rows"; fi
append_summary_row "04_biological_filtering" "poplars" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"

echo "=== HIVSeqinR ==="
OUTDIR="${RESULTS_DIR}/hivseqinr_out"
TIMELOG="${RESULTS_DIR}/hivseqinr.time"
measure_and_run "${TIMELOG}" -- "${STAGE_DIR}/run_hivseqinr.sh" "${INPUT_FASTA}" "${OUTDIR}" > "${RESULTS_DIR}/hivseqinr_wrapper.log" 2>&1
EXIT_CODE=$?
parse_time_metrics "${TIMELOG}"
VALID=0; METRIC="n/a"
CSV="${OUTDIR}/Output_MyBigSummary_DF_FINAL.csv"
if [ -s "${CSV}" ]; then VALID=1; METRIC="$(($(wc -l < "${CSV}") - 1)) classified"; fi
append_summary_row "04_biological_filtering" "hivseqinr" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"

echo "=== HIVIntact ==="
OUTDIR="${RESULTS_DIR}/hivintact_out"
TIMELOG="${RESULTS_DIR}/hivintact.time"
measure_and_run "${TIMELOG}" -- "${STAGE_DIR}/run_hivintact.sh" "${INPUT_FASTA}" "${OUTDIR}" B > "${RESULTS_DIR}/hivintact_wrapper.log" 2>&1
EXIT_CODE=$?
parse_time_metrics "${TIMELOG}"
VALID=0; METRIC="n/a"
if [ -s "${OUTDIR}/intact.fasta" ] || [ -s "${OUTDIR}/nonintact.fasta" ]; then
    VALID=1
    N_INTACT=$(grep -c "^>" "${OUTDIR}/intact.fasta" 2>/dev/null || echo 0)
    N_NONINTACT=$(grep -c "^>" "${OUTDIR}/nonintact.fasta" 2>/dev/null || echo 0)
    METRIC="${N_INTACT} intact, ${N_NONINTACT} non-intact"
fi
append_summary_row "04_biological_filtering" "hivintact" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"

echo "Done. See ${SUMMARY_TSV}"
