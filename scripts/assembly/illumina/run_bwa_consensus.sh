#!/bin/bash
# Reference-mapping assembly: the approach currently used by the production
# pipeline (illumina_u3analysis.sh Step 7) -- bwa mem to HXB2 + bcftools
# consensus. Included here as the baseline to compare SPAdes/SHIVER against.
# Usage: run_bwa_consensus.sh <SRR> <R1.fastq.gz> <R2.fastq.gz> <REF_FASTA> <OUTDIR>
set -uo pipefail
SRR="$1" R1="$2" R2="$3" REF_FASTA="$4" OUTDIR="$5"
mkdir -p "${OUTDIR}"

BAM_OUT="${OUTDIR}/${SRR}.sorted.bam"
VCF_OUT="${OUTDIR}/${SRR}.vcf.gz"
CONSENSUS_OUT="${OUTDIR}/${SRR}_consensus.fasta"

bwa mem -t "${THREADS:-4}" "${REF_FASTA}" "${R1}" "${R2}" 2>"${OUTDIR}/${SRR}_bwa.log" | \
    samtools view -b - | samtools sort -o "${BAM_OUT}" || exit 1
samtools index "${BAM_OUT}" || exit 1

bcftools mpileup -Ou -f "${REF_FASTA}" "${BAM_OUT}" 2>"${OUTDIR}/${SRR}_bcftools.log" | \
    bcftools call -c -Oz -o "${VCF_OUT}" 2>>"${OUTDIR}/${SRR}_bcftools.log" || exit 1
tabix -p vcf "${VCF_OUT}" || exit 1
cat "${REF_FASTA}" | bcftools consensus "${VCF_OUT}" > "${CONSENSUS_OUT}" 2>>"${OUTDIR}/${SRR}_bcftools.log" || exit 1
sed -i "1s/.*/>${SRR}/" "${CONSENSUS_OUT}"
