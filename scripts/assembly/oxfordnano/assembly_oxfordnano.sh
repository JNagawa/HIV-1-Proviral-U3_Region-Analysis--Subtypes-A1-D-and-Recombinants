#!/bin/bash
# For the assembly step of  Nanopore sequences, minimap2 and bcft00ls were used.
# reference-mapping via minimap2 + bcftools consensus)
# Usage: ./assembly_oxfordnano.sh
set -uo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
source "${REPO_ROOT}/scripts/common/lib_compare.sh"

FILTERED_DIR="${REPO_ROOT}/results/download_qc/oxfordnano/porechop_nanofilt_out"
RAW_DIR="${REPO_ROOT}/data/raw/oxnano"
REF_FASTA="${REPO_ROOT}/data/reference/K03455.1.fasta"
RESULTS_DIR="${REPO_ROOT}/results/assembly/oxfordnano"
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
mkdir -p "${RESULTS_DIR}"
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

THREADS="${THREADS:-4}"
export THREADS

run_minimap2_consensus() {
    local srr="$1" fastq="$2" ref_fasta="$3" outdir="$4"
    mkdir -p "${outdir}"
    local bam="${outdir}/${srr}.sorted.bam" vcf="${outdir}/${srr}.vcf.gz" consensus="${outdir}/${srr}_consensus.fasta"

    minimap2 -ax map-ont -t "${THREADS:-4}" "${ref_fasta}" "${fastq}" 2>"${outdir}/${srr}_minimap2.log" | \
        samtools view -b - | samtools sort -o "${bam}" || exit 1
    samtools index "${bam}" || exit 1

    bcftools mpileup -Ou -f "${ref_fasta}" "${bam}" 2>"${outdir}/${srr}_bcftools.log" | \
        bcftools call -c -Oz -o "${vcf}" 2>>"${outdir}/${srr}_bcftools.log" || exit 1
    tabix -p vcf "${vcf}" || exit 1
    cat "${ref_fasta}" | bcftools consensus "${vcf}" > "${consensus}" 2>>"${outdir}/${srr}_bcftools.log" || exit 1
    sed -i "1s/.*/>${srr}/" "${consensus}"
}
export -f run_minimap2_consensus

for SRR in $(subset_accessions nanopore "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do
    FASTQ="${FILTERED_DIR}/${SRR}_filtered.fastq.gz"
    if [ ! -s "${FASTQ}" ]; then
        echo "NOTE: no filtered reads for ${SRR}, run scripts/download_qc/oxfordnano/download_qc_oxfordnano.sh first. Falling back to raw reads." >&2
        FASTQ="${RAW_DIR}/${SRR}.fastq.gz"
    fi
    if [ ! -s "${FASTQ}" ]; then
        echo "WARNING: no reads at all for ${SRR}, skipping." >&2
        continue
    fi

    OUTDIR="${RESULTS_DIR}/minimap2_out"
    TIMELOG="${RESULTS_DIR}/minimap2_${SRR}.time"
    LOG="${RESULTS_DIR}/minimap2_${SRR}.log"
    CONSENSUS="${OUTDIR}/${SRR}_consensus.fasta"

    echo "=== minimap2 consensus on ${SRR} ==="
    measure_and_run "${TIMELOG}" -- \
        bash -c 'run_minimap2_consensus "$@"' _ "${SRR}" "${FASTQ}" "${REF_FASTA}" "${OUTDIR}" > "${LOG}" 2>&1
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

    append_summary_row "assembly_oxfordnano" "minimap2" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"
done

echo "Done. See ${SUMMARY_TSV}"
