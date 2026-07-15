#!/bin/bash
# Compares MAFFT vs MUSCLE vs Clustal Omega, aligning HXB2 + whatever
# per-sample consensus sequences are available from the assembly stage (preferring
# bwa_consensus_out since it doesn't depend on SPAdes/SHIVER succeeding).
# Usage: ./compare_msa_illumina.sh
set -uo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
source "${REPO_ROOT}/scripts/common/lib_compare.sh"

REF_FASTA="${REPO_ROOT}/data/reference/K03455.1.fasta"
ASSEMBLY_DIR="${REPO_ROOT}/results/assembly/illumina/bwa_consensus_out"
RESULTS_DIR="${REPO_ROOT}/results/msa/illumina"
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
mkdir -p "${RESULTS_DIR}"
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

COMBINED="${RESULTS_DIR}/combined_input.fasta"
cat "${REF_FASTA}" "${ASSEMBLY_DIR}"/*_consensus.fasta > "${COMBINED}" 2>/dev/null
N_SEQS=$(grep -c "^>" "${COMBINED}" 2>/dev/null || echo 0)
if [ "${N_SEQS}" -lt 2 ]; then
    echo "ERROR: fewer than 2 sequences available (run scripts/assembly/illumina/compare_assembly_illumina.sh (bwa_consensus) first). Found ${N_SEQS}." >&2
    exit 1
fi
echo "Aligning ${N_SEQS} sequences (HXB2 + subset consensus sequences)."

THREADS="${THREADS:-4}"
export THREADS

for TOOL in mafft muscle clustalo; do
    OUT="${RESULTS_DIR}/${TOOL}_aligned.fasta"
    TIMELOG="${RESULTS_DIR}/${TOOL}.time"
    LOG="${RESULTS_DIR}/${TOOL}.log"

    echo "=== ${TOOL} ==="
    measure_and_run "${TIMELOG}" -- "${STAGE_DIR}/run_${TOOL}.sh" "${COMBINED}" "${OUT}" "${LOG}"
    EXIT_CODE=$?
    parse_time_metrics "${TIMELOG}"

    VALID=0
    METRIC="n/a"
    if [ -s "${OUT}" ]; then
        OUT_SEQS=$(grep -c "^>" "${OUT}")
        COLS=$(seqkit stats -T "${OUT}" 2>/dev/null | tail -1 | cut -f7 | tr -d ',')
        GAP_CHARS=$(grep -v "^>" "${OUT}" | tr -cd '-' | wc -c)
        TOTAL_CHARS=$(grep -v "^>" "${OUT}" | tr -d '\n' | wc -c)
        GAP_PCT=$(awk -v g="${GAP_CHARS}" -v t="${TOTAL_CHARS}" 'BEGIN{if(t>0) printf "%.1f", 100*g/t; else print "0"}')
        if [ "${OUT_SEQS}" = "${N_SEQS}" ]; then
            VALID=1
        fi
        METRIC="${OUT_SEQS}/${N_SEQS} seqs retained, ${COLS:-?} columns, ${GAP_PCT}% gap"
    fi

    append_summary_row "msa_illumina" "${TOOL}" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"
done

echo "Done. See ${SUMMARY_TSV}"
