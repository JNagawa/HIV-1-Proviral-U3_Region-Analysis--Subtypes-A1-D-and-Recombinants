#!/bin/bash
# Reference-mapping assembly: the approach currently used by the production
# pipeline (illumina_u3analysis.sh Step 7) -- bwa mem to HXB2 + bcftools
# consensus. Included here as the baseline to compare SPAdes/SHIVER against.
# Usage: run_bwa_consensus.sh <SRR> <R1.fastq.gz> <R2.fastq.gz> <REF_FASTA> <OUTDIR>
# -u errors on unset vars, pipefail fails a pipe if any stage fails
set -uo pipefail
# positional args: sample id, paired reads, reference, output dir
SRR="$1" R1="$2" R2="$3" REF_FASTA="$4" OUTDIR="$5"
mkdir -p "${OUTDIR}"                                 # ensure the output dir exists

BAM_OUT="${OUTDIR}/${SRR}.sorted.bam"                # sorted alignment output path
VCF_OUT="${OUTDIR}/${SRR}.vcf.gz"                    # called-variants output path
CONSENSUS_OUT="${OUTDIR}/${SRR}_consensus.fasta"     # final consensus FASTA output path

# map reads to HXB2, convert to BAM, sort; bail if any stage fails
bwa mem -t "${THREADS:-4}" "${REF_FASTA}" "${R1}" "${R2}" 2>"${OUTDIR}/${SRR}_bwa.log" | \
    samtools view -b - | samtools sort -o "${BAM_OUT}" || exit 1
samtools index "${BAM_OUT}" || exit 1                # index the sorted BAM (required by mpileup)

# pile up bases and call variants against the reference
bcftools mpileup -Ou -f "${REF_FASTA}" "${BAM_OUT}" 2>"${OUTDIR}/${SRR}_bcftools.log" | \
    bcftools call -c -Oz -o "${VCF_OUT}" 2>>"${OUTDIR}/${SRR}_bcftools.log" || exit 1
# index the VCF so bcftools consensus can read it
tabix -p vcf "${VCF_OUT}" || exit 1
# apply called variants onto the reference to make the consensus
cat "${REF_FASTA}" | bcftools consensus "${VCF_OUT}" > "${CONSENSUS_OUT}" 2>>"${OUTDIR}/${SRR}_bcftools.log" || exit 1
sed -i "1s/.*/>${SRR}/" "${CONSENSUS_OUT}"           # rename the FASTA header to the sample id
