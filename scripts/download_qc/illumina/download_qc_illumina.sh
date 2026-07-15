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
mkdir -p "${RESULTS_DIR}"
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

# Phred quality score Q is defined as Q = -10*log10(P_error) (Ewing & Green
# 1998) -- Q20 corresponds to a 1-in-100 (99%) per-base call accuracy, Q30
# to 1-in-1000 (99.9%). Q20 is used below (sliding-window trim AND
# whole-read average) as the "high-confidence" bar: this pipeline feeds a
# single consensus sequence per sample into subtype calls and TFBS
# motif-mapping in the U3 region, where a handful of miscalled bases can
# flip a motif match -- so both trimmers are held to Q20 rather than
# Trimmomatic's own textbook example (SLIDINGWINDOW:4:15, i.e. Q15/~97%).

# LEADING/TRAILING:3 -- Trimmomatic's own manual example value: cut a base
# from the very start/end of a read the moment its quality drops below Q3.
# Q3 is barely above "no confidence at all" -- this is deliberately lenient,
# just removing the occasional genuinely-unusable edge base; the real
# quality bar is enforced by SLIDINGWINDOW and AVGQUAL below, not this.
TRIM_LEADING=3
TRIM_TRAILING=3
# SLIDINGWINDOW:4:20 -- slide a 4-base window along the read; once the
# window's mean quality falls below Q20 (99% accuracy), trim from there to
# the read's end. 4bp is Trimmomatic's own recommended window size; Q20
# (rather than the manual example's Q15) is our deliberately stricter choice.
# fastp's --cut_right/--cut_right_window_size/--cut_right_mean_quality below
# reproduce this exact same window+threshold for a fair comparison.
TRIM_SLIDINGWINDOW="4:20"
# MINLEN:50 -- discard a read/pair if trimming leaves it under 50bp. Reads
# are 251bp to start, and the reference (HXB2, K03455.1) is only ~9.7kb, so
# very short leftover fragments map ambiguously (multi-mapping) with
# BWA-MEM; 50bp keeps enough sequence for confident, unique placement while
# still retaining reads that only lost their tail to quality trimming.
TRIM_MINLEN=50
# AVGQUAL:20 -- independent of the sliding window, also require the read's
# OVERALL average quality to be >=Q20 after trimming. A read can pass a 4bp
# sliding window check while still being poor quality on average; this is
# the final whole-read quality gate applied to both trimmers.
TRIM_AVGQUAL=20

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
    done
done

echo "Done. See ${SUMMARY_TSV}"
