#!/bin/bash
# Download-QC step of the tool-comparison harness, Nanopore side. Unlike
# download_qc/illumina (fastp vs Trimmomatic), the tools-review document
# offers no competing alternative to Porechop_ABI (adapter trimming) or
# NanoFilt (quality/length filtering) for Nanopore data -- so this is a
# validity check of the single recommended chain, not a "tool A vs tool B"
# comparison. Runs NanoPlot/NanoQC/NanoStat (pre) -> Porechop_ABI -> NanoFilt
# -> NanoPlot/NanoStat (post) -> MultiQC on the Nanopore sample subset.
# Usage: ./download_qc_oxfordnano.sh
set -uo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
source "${REPO_ROOT}/scripts/common/lib_compare.sh"

RAW_DIR="${REPO_ROOT}/data/raw/oxnano"
RESULTS_DIR="${REPO_ROOT}/results/download_qc/oxfordnano"
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
mkdir -p "${RESULTS_DIR}"
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

# Same thresholds as the production pipeline (scripts/pipelines/oxnano_u3analysis.sh) --
# MIN_QUALITY=7 is deliberately more permissive than Illumina's Q20 bar
# (see download_qc/illumina): Nanopore's per-base error rate is inherently
# higher, so a Q20 floor here would discard nearly all reads. MIN_LENGTH=200
# keeps short/degraded fragments out while still retaining partial reads
# (a stricter near-full-length threshold is a legacy/archived option --
# see scripts/archive/README.md's "Open discrepancy" note).
MIN_QUALITY=7
MIN_LENGTH=200
THREADS="${THREADS:-4}"
export THREADS

for SRR in $(subset_accessions nanopore "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do
    FASTQ="${RAW_DIR}/${SRR}.fastq.gz"
    if [ ! -s "${FASTQ}" ]; then
        echo "WARNING: raw FASTQ missing for ${SRR}, skipping." >&2
        continue
    fi

    OUTDIR="${RESULTS_DIR}/porechop_nanofilt_out"
    TIMELOG="${RESULTS_DIR}/porechop_nanofilt_${SRR}.time"
    LOG="${RESULTS_DIR}/porechop_nanofilt_${SRR}.log"
    mkdir -p "${OUTDIR}"
    TRIMMED="${OUTDIR}/${SRR}_trimmed.fastq.gz"
    FILTERED="${OUTDIR}/${SRR}_filtered.fastq.gz"

    echo "=== NanoPlot/NanoQC/NanoStat (pre) on ${SRR} ==="
    NanoPlot --fastq "${FASTQ}" --outdir "${RESULTS_DIR}/nanoplot_pre_out/${SRR}" --prefix "${SRR}_pre_" \
        --threads "${THREADS}" --loglength --plots dot --title "${SRR} - Pre-filtering QC" \
        > "${RESULTS_DIR}/nanoplot_pre_${SRR}.log" 2>&1
    nanoQC -o "${RESULTS_DIR}/nanoqc_out/${SRR}" "${FASTQ}" > "${RESULTS_DIR}/nanoqc_${SRR}.log" 2>&1
    NanoStat --fastq "${FASTQ}" --outdir "${RESULTS_DIR}/nanostat_out" --name "${SRR}_pre_stats.txt" \
        --threads "${THREADS}" > "${RESULTS_DIR}/nanostat_pre_${SRR}.log" 2>&1

    echo "=== Porechop_ABI + NanoFilt on ${SRR} ==="
    measure_and_run "${TIMELOG}" -- bash -c '
        set -uo pipefail
        fastq="$1" trimmed="$2" filtered="$3" min_quality="$4" min_length="$5" threads="$6"
        if porechop_abi --input "${fastq}" --output "${trimmed}" --threads "${threads}"; then
            :
        elif porechop --input "${fastq}" --output "${trimmed}" --threads "${threads}"; then
            :
        else
            echo "WARNING: adapter trimming unavailable, using raw reads for filtering" >&2
            trimmed="${fastq}"
        fi
        gunzip -c "${trimmed}" | NanoFilt --quality "${min_quality}" --length "${min_length}" | gzip > "${filtered}"
    ' _ "${FASTQ}" "${TRIMMED}" "${FILTERED}" "${MIN_QUALITY}" "${MIN_LENGTH}" "${THREADS}" \
        > "${LOG}" 2>&1
    EXIT_CODE=$?
    parse_time_metrics "${TIMELOG}"
    rm -f "${TRIMMED}"

    VALID=0
    METRIC="n/a"
    if [ -s "${FILTERED}" ] && gzip -t "${FILTERED}" 2>/dev/null; then
        RAW_READS=$(zcat "${FASTQ}" | awk 'END{print NR/4}')
        FILT_READS=$(zcat "${FILTERED}" | awk 'END{print NR/4}')
        if [ "${FILT_READS}" -gt 0 ]; then
            VALID=1
            RETAINED=$(awk -v f="${FILT_READS}" -v r="${RAW_READS}" 'BEGIN{printf "%.1f", 100*f/r}')
            METRIC="${FILT_READS}/${RAW_READS} reads retained (${RETAINED}%)"
        fi

        echo "=== NanoPlot/NanoStat (post) on ${SRR} ==="
        NanoPlot --fastq "${FILTERED}" --outdir "${RESULTS_DIR}/nanoplot_post_out/${SRR}" --prefix "${SRR}_post_" \
            --threads "${THREADS}" --loglength --plots dot --title "${SRR} - Post-filtering QC" \
            > "${RESULTS_DIR}/nanoplot_post_${SRR}.log" 2>&1
        NanoStat --fastq "${FILTERED}" --outdir "${RESULTS_DIR}/nanostat_out" --name "${SRR}_post_stats.txt" \
            --threads "${THREADS}" > "${RESULTS_DIR}/nanostat_post_${SRR}.log" 2>&1
    fi

    append_summary_row "download_qc_oxfordnano" "porechop_abi+nanofilt" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"
done

echo "=== MultiQC ==="
multiqc "${RESULTS_DIR}" --outdir "${RESULTS_DIR}/multiqc_out" --filename "nanopore_qc_report" \
    --title "PRJNA765218 - Nanopore QC Summary (subset)" --force > "${RESULTS_DIR}/multiqc.log" 2>&1

echo "Done. See ${SUMMARY_TSV}"
