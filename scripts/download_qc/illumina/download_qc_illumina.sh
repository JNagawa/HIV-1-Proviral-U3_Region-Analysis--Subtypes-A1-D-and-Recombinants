#!/bin/bash
# Download-QC step of the tool-comparison harness: runs FastQC, fastp, and
# Trimmomatic directly, then compares fastp vs Trimmomatic on the Illumina
# sample subset. Both trimmers are held to the same quality/length bar (see
# TRIM_* below) so the comparison reflects the algorithms, not mismatched
# settings -- matches scripts/pipelines/illumina/01_download_qc_trim.sh's
# production thresholds.
# Usage: ./download_qc_illumina.sh   (run from this directory, with HIV_U3analysis activated)
set -uo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
source "${REPO_ROOT}/scripts/common/lib_compare.sh"

RAW_DIR="${REPO_ROOT}/data/raw/illumina"
RESULTS_DIR="${REPO_ROOT}/results/download_qc/illumina"
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
KRAKEN2_DB="${REPO_ROOT}/data/reference/kraken2_standard_16gb_db"
mkdir -p "${RESULTS_DIR}"
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

# Phred quality score Q is defined as Q = -10*log10(P_error) (Ewing & Green
# 1998) -- Q15 corresponds to a 1-in-32 (~97%) per-base call accuracy, Q20
# to 1-in-100 (99%). Q15 is used below (sliding-window trim AND whole-read
# average) as the quality bar, matching Trimmomatic's own textbook
# SLIDINGWINDOW:4:15 example -- both trimmers are held to the same
# threshold so the comparison reflects the algorithms, not mismatched
# settings.

# LEADING/TRAILING:3 -- Trimmomatic's own manual example value: cut a base
# from the very start/end of a read the moment its quality drops below Q3.
# Q3 is barely above "no confidence at all" -- this is deliberately lenient,
# just removing the occasional genuinely-unusable edge base; the real
# quality bar is enforced by SLIDINGWINDOW and AVGQUAL below, not this.
TRIM_LEADING=3
TRIM_TRAILING=3
# SLIDINGWINDOW:4:15 -- slide a 4-base window along the read; once the
# window's mean quality falls below Q15 (~97% accuracy), trim from there to
# the read's end. 4bp is Trimmomatic's own recommended window size; Q15 is
# Trimmomatic's own manual example threshold.
# fastp's --cut_right/--cut_right_window_size/--cut_right_mean_quality below
# reproduce this exact same window+threshold for a fair comparison.
TRIM_SLIDINGWINDOW="4:15"
# MINLEN:50 -- discard a read/pair if trimming leaves it under 50bp. Reads
# are 251bp to start, and the reference (HXB2, K03455.1) is only ~9.7kb, so
# very short leftover fragments map ambiguously (multi-mapping) with
# BWA-MEM; 50bp keeps enough sequence for confident, unique placement while
# still retaining reads that only lost their tail to quality trimming.
TRIM_MINLEN=50
# AVGQUAL:15 -- independent of the sliding window, also require the read's
# OVERALL average quality to be >=Q15 after trimming. A read can pass a 4bp
# sliding window check while still being poor quality on average; this is
# the final whole-read quality gate applied to both trimmers.
TRIM_AVGQUAL=15

ADAPTER_FILE=$(compgen -G "${CONDA_PREFIX}/share/trimmomatic*/adapters/TruSeq3-PE-2.fa" | head -1)
THREADS="${THREADS:-4}"
export THREADS

for SRR in $(subset_accessions illumina "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do
    R1="${RAW_DIR}/${SRR}_1.fastq.gz"
    R2="${RAW_DIR}/${SRR}_2.fastq.gz"
    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ]; then
        echo "WARNING: raw FASTQs missing for ${SRR}, skipping." >&2
        continue
    fi

    # --- Pre-trim FastQC: one shared baseline per sample, not tied to
    # either trimmer (there's only one "before", regardless of which
    # trimmer runs next). ---
    FASTQC_PRE_DIR="${RESULTS_DIR}/fastqc_pre_out"
    mkdir -p "${FASTQC_PRE_DIR}"
    if [ ! -s "${FASTQC_PRE_DIR}/${SRR}_1_fastqc.html" ]; then
        echo "=== FastQC (pre-trim) on ${SRR} ==="
        fastqc "${R1}" "${R2}" --outdir "${FASTQC_PRE_DIR}" --threads "${THREADS}" --quiet \
            > "${RESULTS_DIR}/fastqc_pre_${SRR}.log" 2>&1
    fi

    for TOOL in fastp trimmomatic; do
        OUTDIR="${RESULTS_DIR}/${TOOL}_out"
        TIMELOG="${RESULTS_DIR}/${TOOL}_${SRR}.time"
        LOG="${RESULTS_DIR}/${TOOL}_${SRR}.log"
        mkdir -p "${OUTDIR}"

        echo "=== ${TOOL} on ${SRR} ==="
        if [ "${TOOL}" = "fastp" ]; then
            # fastp's cut_right is off by default -- enable it so it does a
            # real sliding-window trim equivalent to Trimmomatic's
            # SLIDINGWINDOW, rather than being compared on easier defaults.
            measure_and_run "${TIMELOG}" -- \
                fastp \
                    -i "${R1}" -I "${R2}" \
                    -o "${OUTDIR}/${SRR}_1.trimmed.fastq.gz" -O "${OUTDIR}/${SRR}_2.trimmed.fastq.gz" \
                    --cut_right --cut_right_window_size 4 --cut_right_mean_quality "${TRIM_AVGQUAL}" \
                    --average_qual "${TRIM_AVGQUAL}" \
                    --length_required "${TRIM_MINLEN}" \
                    --dedup \
                    --json "${OUTDIR}/${SRR}_fastp.json" --html "${OUTDIR}/${SRR}_fastp.html" \
                    --thread "${THREADS}" \
                > "${LOG}" 2>&1
        else
            TRIM_STEPS=""
            if [ -n "${ADAPTER_FILE}" ] && [ -f "${ADAPTER_FILE}" ]; then
                TRIM_STEPS="ILLUMINACLIP:${ADAPTER_FILE}:2:30:10:2:True "
            fi
            TRIM_STEPS+="LEADING:${TRIM_LEADING} TRAILING:${TRIM_TRAILING} "
            TRIM_STEPS+="SLIDINGWINDOW:${TRIM_SLIDINGWINDOW} "
            TRIM_STEPS+="AVGQUAL:${TRIM_AVGQUAL} "
            TRIM_STEPS+="MINLEN:${TRIM_MINLEN}"

            measure_and_run "${TIMELOG}" -- \
                trimmomatic PE \
                    -Xmx48g \
                    -threads "${THREADS}" \
                    -phred33 \
                    -summary "${OUTDIR}/${SRR}_trimmomatic_summary.txt" \
                    "${R1}" "${R2}" \
                    "${OUTDIR}/${SRR}_1.trimmed.fastq.gz" "${OUTDIR}/${SRR}_1.unpaired.fastq.gz" \
                    "${OUTDIR}/${SRR}_2.trimmed.fastq.gz" "${OUTDIR}/${SRR}_2.unpaired.fastq.gz" \
                    ${TRIM_STEPS} \
                > "${LOG}" 2>&1
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

            # --- Post-trim FastQC on this tool's own output ---
            FASTQC_POST_DIR="${RESULTS_DIR}/fastqc_post_${TOOL}_out"
            mkdir -p "${FASTQC_POST_DIR}"
            fastqc "${R1_OUT}" "${R2_OUT}" --outdir "${FASTQC_POST_DIR}" --threads "${THREADS}" --quiet \
                > "${RESULTS_DIR}/fastqc_post_${TOOL}_${SRR}.log" 2>&1
        fi

        append_summary_row "download_qc_illumina" "${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"

        # --- Deduplication + Kraken2 host/bacterial contamination filtering,
        # applied downstream of whichever trimmer just ran (fastp's own
        # --dedup above already deduplicates; Trimmomatic has no native
        # dedup capability, so a fastp dedup-only pass runs first for that
        # arm). This mirrors the review's fastp -> Kraken2 -> SHIVER chain
        # for both trimmer arms, so the comparison stays fair. ---
        if [ "${VALID}" -eq 1 ]; then
            DEDUP_R1="${R1_OUT}" DEDUP_R2="${R2_OUT}"
            if [ "${TOOL}" = "trimmomatic" ]; then
                DEDUP_R1="${OUTDIR}/${SRR}_1.dedup.fastq.gz"
                DEDUP_R2="${OUTDIR}/${SRR}_2.dedup.fastq.gz"
                bash "${REPO_ROOT}/scripts/utils/fastp_dedup.sh" \
                    "${R1_OUT}" "${R2_OUT}" "${DEDUP_R1}" "${DEDUP_R2}" \
                    "${OUTDIR}/${SRR}_dedup_fastp" \
                    > "${RESULTS_DIR}/dedup_${TOOL}_${SRR}.log" 2>&1
            fi

            KRAKEN_OUTDIR="${RESULTS_DIR}/kraken2_${TOOL}_out"
            KRAKEN_TIMELOG="${RESULTS_DIR}/kraken2_${TOOL}_${SRR}.time"
            measure_and_run "${KRAKEN_TIMELOG}" -- \
                bash "${REPO_ROOT}/scripts/utils/kraken2_filter_reads.sh" \
                    "${DEDUP_R1}" "${DEDUP_R2}" "${KRAKEN_OUTDIR}" "${SRR}" "${KRAKEN2_DB}" \
                > "${RESULTS_DIR}/kraken2_${TOOL}_${SRR}.log" 2>&1
            KRAKEN_EXIT=$?
            parse_time_metrics "${KRAKEN_TIMELOG}"

            KRAKEN_VALID=0
            KRAKEN_METRIC="n/a"
            KRAKEN_R1_OUT="${KRAKEN_OUTDIR}/${SRR}_1.kraken_filtered.fastq.gz"
            KRAKEN_R2_OUT="${KRAKEN_OUTDIR}/${SRR}_2.kraken_filtered.fastq.gz"
            if [ -s "${KRAKEN_R1_OUT}" ] && [ -s "${KRAKEN_R2_OUT}" ]; then
                KRAKEN_VALID=1
                KRAKEN_METRIC=$(grep -o '^Kraken2 filtering.*retained[^.]*' "${RESULTS_DIR}/kraken2_${TOOL}_${SRR}.log" | tail -1)
                [ -z "${KRAKEN_METRIC}" ] && KRAKEN_METRIC="filtered, see ${KRAKEN_OUTDIR}/${SRR}.kreport"
            fi
            append_summary_row "download_qc_illumina" "kraken2_after_${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${KRAKEN_EXIT}" "${KRAKEN_VALID}" "${KRAKEN_METRIC}"
        fi
    done
done

echo "Done. See ${SUMMARY_TSV}"
