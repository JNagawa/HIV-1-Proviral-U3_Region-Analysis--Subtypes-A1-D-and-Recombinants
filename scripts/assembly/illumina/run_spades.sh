#!/bin/bash
# De novo assembly with SPAdes, then pick the single best-matching contig
# against HXB2 (by blast) as the "assembly" output for comparison.
# Usage: run_spades.sh <SRR> <R1.fastq.gz> <R2.fastq.gz> <REF_FASTA> <OUTDIR>
# -u errors on unset vars, pipefail fails a pipe if any stage fails
set -uo pipefail
# positional args: sample id, paired reads, reference, output dir
SRR="$1" R1="$2" R2="$3" REF_FASTA="$4" OUTDIR="$5"
SPADES_DIR="${OUTDIR}/${SRR}_spades"                 # per-sample SPAdes working dir
mkdir -p "${OUTDIR}"                                 # ensure the output dir exists

# assemble reads into contigs; --careful reduces mismatches/indels
spades.py --careful -1 "${R1}" -2 "${R2}" -o "${SPADES_DIR}" -t "${THREADS:-4}" \
    > "${OUTDIR}/${SRR}_spades.log" 2>&1
CONTIGS="${SPADES_DIR}/contigs.fasta"                # SPAdes' contig output file
# fail if no contigs were produced
[ -s "${CONTIGS}" ] || { echo "ERROR: SPAdes produced no contigs for ${SRR}" >&2; exit 1; }

# Pick the contig with the best blast hit against HXB2 as the near-full-genome
# candidate (SPAdes commonly returns many short/host-contaminant contigs).
# only do blast selection if the blast tools exist
if command -v makeblastdb >/dev/null 2>&1 && command -v blastn >/dev/null 2>&1; then
    # build a blast DB from HXB2
    makeblastdb -in "${REF_FASTA}" -dbtype nucl -out "${OUTDIR}/${SRR}_hxb2db" >/dev/null 2>&1
    # blast contigs vs HXB2, sort by bitscore, take the top hit's contig id
    BEST_ID=$(blastn -query "${CONTIGS}" -db "${OUTDIR}/${SRR}_hxb2db" -outfmt "6 qseqid length bitscore" 2>/dev/null \
        | sort -k3,3 -rn | head -1 | cut -f1)
    if [ -n "${BEST_ID}" ]; then                     # if a best contig was found...
        # extract that one contig as the assembly
        seqkit grep -n -p "${BEST_ID}" "${CONTIGS}" > "${OUTDIR}/${SRR}_consensus.fasta"
        # rename the FASTA header to the sample id
        sed -i "1s/.*/>${SRR}/" "${OUTDIR}/${SRR}_consensus.fasta"
        exit 0                                       # done -- skip the length-based fallback
    fi
fi

# Fallback (no blast available): just take the longest contig.
# sort contigs by length desc, keep the longest
seqkit sort -l -r "${CONTIGS}" 2>/dev/null | seqkit head -n 1 > "${OUTDIR}/${SRR}_consensus.fasta"
# rename the FASTA header to the sample id
sed -i "1s/.*/>${SRR}/" "${OUTDIR}/${SRR}_consensus.fasta"
