#!/bin/bash
# PacBio arm: compares Poplars (Hypermut 3) vs HIVSeqinR vs HIVIntact on the
# PacBio assembled subset (HXB2 + per-sample proviral consensus from
# minimap2_consensus_out). Reuses the Illumina arm's run_<tool>.sh scripts and
# the shared scripts/tools/{Poplars,HIVSeqinR,HIVIntact} clones -- intactness/
# hypermutation classification is platform-agnostic (operates on FASTA
# assemblies), only the input consensus sequences differ. As on the Illumina
# arm, no gold standard exists for this cohort, so "validity" is just "ran and
# produced a classification"; the review's recommended order is Poplars then
# HIVSeqinR, with HIVIntact as an optional cross-check (subtype-B-biased, so
# read with caution on A1/D).
# Usage: ./compare_biological_filtering_pacbio.sh   (via sbatch scripts/utils/run_comparison_step.slurm.sh)
set -uo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
source "${REPO_ROOT}/scripts/common/lib_compare.sh"
RUN_DIR="${REPO_ROOT}/scripts/biological_filtering/illumina"   # reuse shared run_<tool>.sh

REF_FASTA="${REPO_ROOT}/data/reference/K03455.1.fasta"
ASSEMBLY_DIR="${REPO_ROOT}/results/assembly/pacbio/minimap2_consensus_out"
RESULTS_DIR="${REPO_ROOT}/results/biological_filtering/pacbio"
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
mkdir -p "${RESULTS_DIR}"
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

INPUT_FASTA="${RESULTS_DIR}/combined_input.fasta"
cat "${REF_FASTA}" "${ASSEMBLY_DIR}"/*_consensus.fasta > "${INPUT_FASTA}" 2>/dev/null
N_SEQS=$(grep -c "^>" "${INPUT_FASTA}" 2>/dev/null || echo 0)
if [ "${N_SEQS}" -lt 2 ]; then
    echo "ERROR: fewer than 2 sequences available (run scripts/assembly/pacbio/compare_assembly_pacbio.sh first). Found ${N_SEQS}." >&2
    exit 1
fi
if [ ! -d "${REPO_ROOT}/scripts/tools/Poplars" ]; then
    echo "NOTE: tools not set up yet -- run scripts/biological_filtering/illumina/setup_tools.sh first." >&2
    exit 1
fi

echo "=== Poplars (Hypermut 3) ==="
OUT="${RESULTS_DIR}/poplars_out.tsv"
TIMELOG="${RESULTS_DIR}/poplars.time"
measure_and_run "${TIMELOG}" -- "${RUN_DIR}/run_poplars.sh" "${INPUT_FASTA}" "${OUT}" > "${RESULTS_DIR}/poplars.log" 2>&1
EXIT_CODE=$?
parse_time_metrics "${TIMELOG}"
VALID=0; METRIC="n/a"
if [ -s "${OUT}" ]; then VALID=1; METRIC="$(wc -l < "${OUT}") result rows"; fi
append_summary_row "biological_filtering_pacbio" "poplars" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"

echo "=== HIVSeqinR ==="
OUTDIR="${RESULTS_DIR}/hivseqinr_out"
TIMELOG="${RESULTS_DIR}/hivseqinr.time"
measure_and_run "${TIMELOG}" -- "${RUN_DIR}/run_hivseqinr.sh" "${INPUT_FASTA}" "${OUTDIR}" > "${RESULTS_DIR}/hivseqinr_wrapper.log" 2>&1
EXIT_CODE=$?
parse_time_metrics "${TIMELOG}"
VALID=0; METRIC="n/a"
CSV="${OUTDIR}/Output_MyBigSummary_DF_FINAL.csv"
if [ -s "${CSV}" ]; then VALID=1; METRIC="$(($(wc -l < "${CSV}") - 1)) classified"; fi
append_summary_row "biological_filtering_pacbio" "hivseqinr" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"

echo "=== HIVIntact ==="
OUTDIR="${RESULTS_DIR}/hivintact_out"
TIMELOG="${RESULTS_DIR}/hivintact.time"
measure_and_run "${TIMELOG}" -- "${RUN_DIR}/run_hivintact.sh" "${INPUT_FASTA}" "${OUTDIR}" A1 > "${RESULTS_DIR}/hivintact_wrapper.log" 2>&1
EXIT_CODE=$?
parse_time_metrics "${TIMELOG}"
VALID=0; METRIC="n/a"
if [ -s "${OUTDIR}/intact.fasta" ] || [ -s "${OUTDIR}/nonintact.fasta" ]; then
    VALID=1
    N_INTACT=$(grep -c "^>" "${OUTDIR}/intact.fasta" 2>/dev/null); N_INTACT="${N_INTACT:-0}"
    N_NONINTACT=$(grep -c "^>" "${OUTDIR}/nonintact.fasta" 2>/dev/null); N_NONINTACT="${N_NONINTACT:-0}"
    METRIC="${N_INTACT} intact, ${N_NONINTACT} non-intact"
fi
append_summary_row "biological_filtering_pacbio" "hivintact" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"

echo "Done. See ${SUMMARY_TSV}"
