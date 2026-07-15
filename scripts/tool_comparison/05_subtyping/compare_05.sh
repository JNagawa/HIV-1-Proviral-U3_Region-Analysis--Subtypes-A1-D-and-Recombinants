#!/bin/bash
# Compares jpHMM (per-sample recombination/subtype calls) vs IQ-TREE2
# (ML phylogeny to confirm subtype clustering) on the Illumina subset.
# COMET and REGA v3 are excluded (web-only, run those separately).
# Usage: ./compare_05.sh
set -uo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
source "${REPO_ROOT}/scripts/tool_comparison/common/lib_compare.sh"

ASSEMBLY_DIR="${REPO_ROOT}/results/tool_comparison/02_assembly_illumina/bwa_consensus_out"
ALIGNMENT="${REPO_ROOT}/results/tool_comparison/03_msa/mafft_aligned.fasta"
RESULTS_DIR="${REPO_ROOT}/results/tool_comparison/05_subtyping"
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
mkdir -p "${RESULTS_DIR}"
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/tool_comparison/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

THREADS="${THREADS:-4}"
export THREADS

echo "=== jpHMM (per sample) ==="
for SRR in $(subset_accessions illumina "${REPO_ROOT}/scripts/tool_comparison/subset_samples.tsv"); do
    QUERY="${ASSEMBLY_DIR}/${SRR}_consensus.fasta"
    [ -s "${QUERY}" ] || { echo "NOTE: no consensus for ${SRR}, skipping (run stage 02 first)." >&2; continue; }

    OUTDIR="${RESULTS_DIR}/jphmm_out/${SRR}"
    TIMELOG="${RESULTS_DIR}/jphmm_${SRR}.time"
    LOG="${RESULTS_DIR}/jphmm_${SRR}.log"

    measure_and_run "${TIMELOG}" -- "${STAGE_DIR}/run_jphmm.sh" "${QUERY}" "${OUTDIR}" > "${LOG}" 2>&1
    EXIT_CODE=$?
    parse_time_metrics "${TIMELOG}"

    VALID=0
    METRIC="n/a"
    RECOMB_FILE=$(compgen -G "${OUTDIR}/*recombination*" | head -1)
    if [ -n "${RECOMB_FILE}" ] && [ -s "${RECOMB_FILE}" ]; then
        VALID=1
        METRIC="see ${RECOMB_FILE}"
    fi
    append_summary_row "05_subtyping" "jphmm" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"
done

echo "=== IQ-TREE2 (whole-subset ML tree, from stage 03's MAFFT alignment) ==="
if [ -s "${ALIGNMENT}" ]; then
    TIMELOG="${RESULTS_DIR}/iqtree.time"
    LOG="${RESULTS_DIR}/iqtree.log"
    measure_and_run "${TIMELOG}" -- "${STAGE_DIR}/run_iqtree.sh" "${ALIGNMENT}" "${RESULTS_DIR}/iqtree_out" subset > "${LOG}" 2>&1
    EXIT_CODE=$?
    parse_time_metrics "${TIMELOG}"

    VALID=0
    METRIC="n/a"
    TREEFILE="${RESULTS_DIR}/iqtree_out/subset.treefile"
    if [ -s "${TREEFILE}" ]; then
        VALID=1
        METRIC="tree written, $(grep -o "bp\." "${RESULTS_DIR}/iqtree_out/subset.iqtree" 2>/dev/null | wc -l) bootstrap refs"
    fi
    append_summary_row "05_subtyping" "iqtree2" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"
else
    echo "NOTE: no alignment at ${ALIGNMENT}, run stage 03 first." >&2
fi

echo "Done. See ${SUMMARY_TSV}"
