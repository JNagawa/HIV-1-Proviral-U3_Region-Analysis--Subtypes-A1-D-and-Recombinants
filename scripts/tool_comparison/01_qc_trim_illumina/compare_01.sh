#!/bin/bash
# Compares fastp vs Trimmomatic on the Illumina sample subset.
# Usage: ./compare_01.sh   (run from this directory, with HIV_U3analysis activated)
set -uo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
source "${REPO_ROOT}/scripts/tool_comparison/common/lib_compare.sh"

RAW_DIR="${REPO_ROOT}/data/raw/illumina"
RESULTS_DIR="${REPO_ROOT}/results/tool_comparison/01_qc_trim_illumina"
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
mkdir -p "${RESULTS_DIR}"
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/tool_comparison/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

ADAPTER_FILE=$(compgen -G "${CONDA_PREFIX}/share/trimmomatic*/adapters/TruSeq3-PE-2.fa" | head -1)
THREADS="${THREADS:-4}"
export THREADS

for SRR in $(subset_accessions illumina "${REPO_ROOT}/scripts/tool_comparison/subset_samples.tsv"); do
    R1="${RAW_DIR}/${SRR}_1.fastq.gz"
    R2="${RAW_DIR}/${SRR}_2.fastq.gz"
    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ]; then
        echo "WARNING: raw FASTQs missing for ${SRR}, skipping." >&2
        continue
    fi

    for TOOL in fastp trimmomatic; do
        OUTDIR="${RESULTS_DIR}/${TOOL}_out"
        TIMELOG="${RESULTS_DIR}/${TOOL}_${SRR}.time"
        LOG="${RESULTS_DIR}/${TOOL}_${SRR}.log"

        echo "=== ${TOOL} on ${SRR} ==="
        if [ "${TOOL}" = "fastp" ]; then
            measure_and_run "${TIMELOG}" -- \
                "${STAGE_DIR}/run_fastp.sh" "${SRR}" "${R1}" "${R2}" "${OUTDIR}" > "${LOG}" 2>&1
        else
            measure_and_run "${TIMELOG}" -- \
                "${STAGE_DIR}/run_trimmomatic.sh" "${SRR}" "${R1}" "${R2}" "${OUTDIR}" "${ADAPTER_FILE}" > "${LOG}" 2>&1
        fi
        EXIT_CODE=$?
        parse_time_metrics "${TIMELOG}"

        # Validity check: R1/R2 output pair counts match and are non-zero;
        # pull the relevant surviving-reads metric from each tool's own report.
        R1_OUT="${OUTDIR}/${SRR}_1.trimmed.fastq.gz"
        R2_OUT="${OUTDIR}/${SRR}_2.trimmed.fastq.gz"
        VALID=0
        METRIC="n/a"
        if [ -s "${R1_OUT}" ] && [ -s "${R2_OUT}" ] && gzip -t "${R1_OUT}" 2>/dev/null && gzip -t "${R2_OUT}" 2>/dev/null; then
            R1_READS=$(zcat "${R1_OUT}" | wc -l)
            R2_READS=$(zcat "${R2_OUT}" | wc -l)
            if [ "${R1_READS}" = "${R2_READS}" ] && [ "${R1_READS}" -gt 0 ]; then
                VALID=1
            fi
            if [ "${TOOL}" = "fastp" ] && [ -s "${OUTDIR}/${SRR}_fastp.json" ]; then
                PASSED=$(grep -o '"passed_filter_reads"[^,}]*' "${OUTDIR}/${SRR}_fastp.json" | head -1 | grep -o '[0-9]*$')
                [ -n "${PASSED}" ] && METRIC="${PASSED} reads passed"
            elif [ "${TOOL}" = "trimmomatic" ] && [ -f "${OUTDIR}/${SRR}_trimmomatic_summary.txt" ]; then
                METRIC=$(grep "Both Surviving Read Percent" "${OUTDIR}/${SRR}_trimmomatic_summary.txt" | awk '{print $NF"% surviving"}')
            fi
        fi

        append_summary_row "01_qc_trim_illumina" "${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"
    done
done

echo "Done. See ${SUMMARY_TSV}"
