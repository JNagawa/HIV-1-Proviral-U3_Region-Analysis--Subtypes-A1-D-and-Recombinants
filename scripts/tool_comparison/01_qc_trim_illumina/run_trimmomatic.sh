#!/bin/bash
# Usage: run_trimmomatic.sh <SRR> <R1.fastq.gz> <R2.fastq.gz> <OUTDIR> <ADAPTER_FILE>
set -uo pipefail
SRR="$1" R1="$2" R2="$3" OUTDIR="$4" ADAPTER_FILE="${5:-}"
mkdir -p "${OUTDIR}"

TRIM_STEPS=""
if [ -n "${ADAPTER_FILE}" ] && [ -f "${ADAPTER_FILE}" ]; then
    TRIM_STEPS="ILLUMINACLIP:${ADAPTER_FILE}:2:30:10:2:True "
fi
TRIM_STEPS+="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 AVGQUAL:20 MINLEN:50"

# -Xmx48g: see illumina_u3analysis.sh -- the bioconda wrapper defaults to a
# 1GB heap regardless of available system memory, which crashes on datasets
# this size.
trimmomatic PE \
    -Xmx48g \
    -threads "${THREADS:-4}" \
    -phred33 \
    -summary "${OUTDIR}/${SRR}_trimmomatic_summary.txt" \
    "${R1}" "${R2}" \
    "${OUTDIR}/${SRR}_1.trimmed.fastq.gz" "${OUTDIR}/${SRR}_1.unpaired.fastq.gz" \
    "${OUTDIR}/${SRR}_2.trimmed.fastq.gz" "${OUTDIR}/${SRR}_2.unpaired.fastq.gz" \
    ${TRIM_STEPS}
