#!/bin/bash
# Deduplicate paired-end reads with fastp, without any additional
# trimming/filtering (that's already been done by whichever trimmer ran
# before this). Used after Trimmomatic, which has no native deduplication
# capability -- when fastp itself is the trimmer, its own --dedup flag is
# used directly instead of this script, since a second fastp pass would be
# redundant.
# Usage: fastp_dedup.sh <R1.fastq.gz> <R2.fastq.gz> <R1_OUT> <R2_OUT> <REPORT_PREFIX>
set -uo pipefail
R1="$1" R2="$2" R1_OUT="$3" R2_OUT="$4" REPORT_PREFIX="$5"

fastp \
    -i "${R1}" -I "${R2}" \
    -o "${R1_OUT}" -O "${R2_OUT}" \
    --dedup \
    --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering \
    --json "${REPORT_PREFIX}.json" --html "${REPORT_PREFIX}.html" \
    --thread "${THREADS:-4}"
