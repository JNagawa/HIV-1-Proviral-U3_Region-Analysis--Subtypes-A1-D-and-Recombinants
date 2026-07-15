#!/bin/bash
# Compares FIMO vs MOODS vs TFBSTools -- PWM motif scanning for the 6 core
# JASPAR TFs (see setup_jaspar.sh) against U3 sequences extracted from
# stage 03's alignment via scripts/utils/extract_u3_by_hxb2_anchor.sh.
#
# Per the README, each tool runs a positive control first: scanning HXB2's
# own U3 alone, where the target motifs are known/expected to be present.
# A tool that finds zero hits on its own positive control is flagged before
# its subset numbers are trusted.
# Usage: ./compare_06.sh
set -uo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
source "${REPO_ROOT}/scripts/tool_comparison/common/lib_compare.sh"

ALIGNMENT="${REPO_ROOT}/results/tool_comparison/03_msa/mafft_aligned.fasta"
RESULTS_DIR="${REPO_ROOT}/results/tool_comparison/06_motif_mapping"
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
mkdir -p "${RESULTS_DIR}"
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/tool_comparison/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

if [ ! -s "${ALIGNMENT}" ]; then
    echo "ERROR: no alignment at ${ALIGNMENT} (run stage 03 first)." >&2
    exit 1
fi

echo "=== Extracting U3 (anchored to HXB2's own annotated LTR/R-region) ==="
U3_FASTA="${RESULTS_DIR}/U3_extracted.fasta"
bash "${REPO_ROOT}/scripts/utils/extract_u3_by_hxb2_anchor.sh" \
    --alignment "${ALIGNMENT}" \
    --hxb2-id K03455.1 \
    --gb-cache "${REPO_ROOT}/data/reference/K03455.1.gb" \
    --out-gapped "${RESULTS_DIR}/U3_aligned.fasta" \
    --out "${U3_FASTA}" \
    --warnings-log "${RESULTS_DIR}/u3_extraction_warnings.log"
if [ ! -s "${U3_FASTA}" ]; then
    echo "ERROR: U3 extraction produced no sequences, see extraction output above." >&2
    exit 1
fi

HXB2_U3_ONLY="${RESULTS_DIR}/hxb2_u3_only.fasta"
seqkit grep -n -r -p "^K03455" "${U3_FASTA}" > "${HXB2_U3_ONLY}"
if [ ! -s "${HXB2_U3_ONLY}" ]; then
    echo "ERROR: could not isolate HXB2's own U3 from ${U3_FASTA} for the positive control." >&2
    exit 1
fi

run_tool() {
    local tool="$1" sample="$2" input="$3" out="$4"
    local timelog="${RESULTS_DIR}/${tool}_${sample}.time"
    local log="${RESULTS_DIR}/${tool}_${sample}.log"
    local exit_code n_hits valid metric

    case "${tool}" in
        fimo)
            measure_and_run "${timelog}" -- "${STAGE_DIR}/run_fimo.sh" "${input}" "${out}" > "${log}" 2>&1
            exit_code=$?
            n_hits=0
            [ -s "${out}/fimo.tsv" ] && n_hits=$(grep -vc "^#\|^motif_id" "${out}/fimo.tsv" 2>/dev/null || echo 0)
            ;;
        moods)
            measure_and_run "${timelog}" -- "${STAGE_DIR}/run_moods.sh" "${input}" "${out}" > "${log}" 2>&1
            exit_code=$?
            n_hits=0
            [ -s "${out}" ] && n_hits=$(wc -l < "${out}")
            ;;
        tfbstools)
            measure_and_run "${timelog}" -- "${STAGE_DIR}/run_tfbstools.sh" "${input}" "${out}" > "${log}" 2>&1
            exit_code=$?
            n_hits=0
            [ -s "${out}" ] && n_hits=$(grep -vc "^#" "${out}" 2>/dev/null || echo 0)
            ;;
    esac

    parse_time_metrics "${timelog}"
    valid=0
    metric="0 hits"
    if [ "${n_hits}" -gt 0 ] 2>/dev/null; then
        valid=1
        metric="${n_hits} hits"
    fi
    append_summary_row "06_motif_mapping" "${tool}" "${sample}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${exit_code}" "${valid}" "${metric}"

    if [ "${sample}" = "positive_control" ] && [ "${valid}" -eq 0 ]; then
        echo "WARNING: ${tool} found zero motif hits on HXB2's own U3 (positive control failed) -- treat its subset numbers as suspect until this is investigated." >&2
    fi
}

for TOOL in fimo moods tfbstools; do
    echo "=== ${TOOL}: positive control (HXB2 own U3) ==="
    case "${TOOL}" in
        fimo) PC_OUT="${RESULTS_DIR}/fimo_positive_control_out" ;;
        moods) PC_OUT="${RESULTS_DIR}/moods_positive_control.tsv" ;;
        tfbstools) PC_OUT="${RESULTS_DIR}/tfbstools_positive_control.gff3" ;;
    esac
    run_tool "${TOOL}" "positive_control" "${HXB2_U3_ONLY}" "${PC_OUT}"

    echo "=== ${TOOL}: full subset ==="
    case "${TOOL}" in
        fimo) OUT="${RESULTS_DIR}/fimo_out" ;;
        moods) OUT="${RESULTS_DIR}/moods_out.tsv" ;;
        tfbstools) OUT="${RESULTS_DIR}/tfbstools_out.gff3" ;;
    esac
    run_tool "${TOOL}" "all_subset" "${U3_FASTA}" "${OUT}"
done

echo "Done. See ${SUMMARY_TSV}"
