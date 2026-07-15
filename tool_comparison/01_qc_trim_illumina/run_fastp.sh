#!/bin/bash
# Usage: run_fastp.sh <SRR> <R1.fastq.gz> <R2.fastq.gz> <OUTDIR>
set -uo pipefail
SRR="$1" R1="$2" R2="$3" OUTDIR="$4"
mkdir -p "${OUTDIR}"

fastp \
    -i "${R1}" -I "${R2}" \
    -o "${OUTDIR}/${SRR}_1.trimmed.fastq.gz" -O "${OUTDIR}/${SRR}_2.trimmed.fastq.gz" \
    --json "${OUTDIR}/${SRR}_fastp.json" --html "${OUTDIR}/${SRR}_fastp.html" \
    --thread "${THREADS:-4}"
