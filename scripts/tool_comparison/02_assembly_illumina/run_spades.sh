#!/bin/bash
# De novo assembly with SPAdes, then pick the single best-matching contig
# against HXB2 (by blast) as the "assembly" output for comparison.
# Usage: run_spades.sh <SRR> <R1.fastq.gz> <R2.fastq.gz> <REF_FASTA> <OUTDIR>
set -uo pipefail
SRR="$1" R1="$2" R2="$3" REF_FASTA="$4" OUTDIR="$5"
SPADES_DIR="${OUTDIR}/${SRR}_spades"
mkdir -p "${OUTDIR}"

spades.py --careful -1 "${R1}" -2 "${R2}" -o "${SPADES_DIR}" -t "${THREADS:-4}" \
    > "${OUTDIR}/${SRR}_spades.log" 2>&1
CONTIGS="${SPADES_DIR}/contigs.fasta"
[ -s "${CONTIGS}" ] || { echo "ERROR: SPAdes produced no contigs for ${SRR}" >&2; exit 1; }

# Pick the contig with the best blast hit against HXB2 as the near-full-genome
# candidate (SPAdes commonly returns many short/host-contaminant contigs).
if command -v makeblastdb >/dev/null 2>&1 && command -v blastn >/dev/null 2>&1; then
    makeblastdb -in "${REF_FASTA}" -dbtype nucl -out "${OUTDIR}/${SRR}_hxb2db" >/dev/null 2>&1
    BEST_ID=$(blastn -query "${CONTIGS}" -db "${OUTDIR}/${SRR}_hxb2db" -outfmt "6 qseqid length bitscore" 2>/dev/null \
        | sort -k3,3 -rn | head -1 | cut -f1)
    if [ -n "${BEST_ID}" ]; then
        seqkit grep -n -p "${BEST_ID}" "${CONTIGS}" > "${OUTDIR}/${SRR}_consensus.fasta"
        sed -i "1s/.*/>${SRR}/" "${OUTDIR}/${SRR}_consensus.fasta"
        exit 0
    fi
fi

# Fallback (no blast available): just take the longest contig.
seqkit sort -l -r "${CONTIGS}" 2>/dev/null | seqkit head -n 1 > "${OUTDIR}/${SRR}_consensus.fasta"
sed -i "1s/.*/>${SRR}/" "${OUTDIR}/${SRR}_consensus.fasta"
