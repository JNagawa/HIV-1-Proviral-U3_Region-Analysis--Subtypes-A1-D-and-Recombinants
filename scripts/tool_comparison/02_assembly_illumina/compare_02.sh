#!/bin/bash
# Compares BWA+bcftools-consensus vs SPAdes vs SHIVER on the Illumina subset.
# Uses each tool's TRIMMED reads from stage 01 if available (fastp output,
# since the review recommends fastp before assembly), else falls back to raw
# reads with a warning.
# Usage: ./compare_02.sh
set -uo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
source "${REPO_ROOT}/scripts/tool_comparison/common/lib_compare.sh"

RAW_DIR="${REPO_ROOT}/data/raw/illumina"
FASTP_DIR="${REPO_ROOT}/results/tool_comparison/01_qc_trim_illumina/fastp_out"
REF_FASTA="${REPO_ROOT}/data/reference/K03455.1.fasta"
RESULTS_DIR="${REPO_ROOT}/results/tool_comparison/02_assembly_illumina"
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
mkdir -p "${RESULTS_DIR}"
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/tool_comparison/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

THREADS="${THREADS:-4}"
export THREADS

for SRR in $(subset_accessions illumina "${REPO_ROOT}/scripts/tool_comparison/subset_samples.tsv"); do
    R1="${FASTP_DIR}/${SRR}_1.trimmed.fastq.gz"
    R2="${FASTP_DIR}/${SRR}_2.trimmed.fastq.gz"
    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ]; then
        echo "NOTE: no fastp-trimmed reads for ${SRR}, run stage 01 first for the intended input. Falling back to raw reads." >&2
        R1="${RAW_DIR}/${SRR}_1.fastq.gz"
        R2="${RAW_DIR}/${SRR}_2.fastq.gz"
    fi
    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ]; then
        echo "WARNING: no reads at all for ${SRR}, skipping." >&2
        continue
    fi

    for TOOL in bwa_consensus spades shiver; do
        OUTDIR="${RESULTS_DIR}/${TOOL}_out"
        TIMELOG="${RESULTS_DIR}/${TOOL}_${SRR}.time"
        LOG="${RESULTS_DIR}/${TOOL}_${SRR}.log"
        CONSENSUS="${OUTDIR}/${SRR}_consensus.fasta"

        echo "=== ${TOOL} on ${SRR} ==="
        case "${TOOL}" in
            bwa_consensus)
                measure_and_run "${TIMELOG}" -- \
                    "${STAGE_DIR}/run_bwa_consensus.sh" "${SRR}" "${R1}" "${R2}" "${REF_FASTA}" "${OUTDIR}" > "${LOG}" 2>&1 ;;
            spades)
                measure_and_run "${TIMELOG}" -- \
                    "${STAGE_DIR}/run_spades.sh" "${SRR}" "${R1}" "${R2}" "${REF_FASTA}" "${OUTDIR}" > "${LOG}" 2>&1 ;;
            shiver)
                measure_and_run "${TIMELOG}" -- \
                    "${STAGE_DIR}/run_shiver.sh" "${SRR}" "${R1}" "${R2}" "${OUTDIR}" > "${LOG}" 2>&1 ;;
        esac
        EXIT_CODE=$?
        parse_time_metrics "${TIMELOG}"

        VALID=0
        METRIC="n/a"
        if [ -s "${CONSENSUS}" ]; then
            LEN=$(seqkit stats -T "${CONSENSUS}" 2>/dev/null | tail -1 | cut -f5)
            N_PCT=$(seqkit fx2tab -n -g -B N "${CONSENSUS}" 2>/dev/null | awk -F'\t' '{print $NF}' | tail -1)
            if [ -n "${LEN}" ] && [ "${LEN}" -ge 8000 ] && [ "${LEN}" -le 10000 ]; then
                VALID=1
            fi
            METRIC="${LEN}bp, ${N_PCT:-?}% N"
        fi

        append_summary_row "02_assembly_illumina" "${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"
    done
done

echo "Done. See ${SUMMARY_TSV}"
