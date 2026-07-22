#!/bin/bash
#SBATCH --job-name=pacbio_download
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# Download the PacBio HIV-SMRTcap subset (the `pacbio` rows of
# scripts/common/subset_samples.tsv) from SRA into data/raw/pacbio/ as
# single-end HiFi FASTQs. Unlike the Illumina/Nanopore arms (already on
# disk), the SMRTcap data is fetched here, so this is its own sbatch job
# rather than folded into download_qc -- these are multi-GB HiFi sets and
# must not run on the login node.
# Usage: sbatch scripts/download_qc/pacbio/download_pacbio.slurm.sh
set -euo pipefail

CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"
if [ -f "$CONDA_SH" ]; then source "$CONDA_SH"; else source "$(conda info --base)/etc/profile.d/conda.sh"; fi
conda activate HIV_U3analysis

cd /etc/ace-data/home/jnagawa/Internship
source scripts/common/lib_compare.sh

RAW_DIR="data/raw/pacbio"
PREFETCH_DIR="data/raw/pacbio/.sra"
mkdir -p "${RAW_DIR}" "${PREFETCH_DIR}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

for SRR in $(subset_accessions pacbio scripts/common/subset_samples.tsv); do
    FINAL="${RAW_DIR}/${SRR}.fastq.gz"
    if [ -s "${FINAL}" ]; then
        echo "=== ${SRR} already present (${FINAL}), skipping ==="
        continue
    fi
    echo "=== prefetch ${SRR} ==="
    # prefetch's default cap is 20G; set 100G explicitly so the subset (and a
    # Revio run, ~18GB, if added later) is never silently skipped.
    prefetch --max-size 100G --output-directory "${PREFETCH_DIR}" "${SRR}"

    echo "=== fasterq-dump ${SRR} (single-end HiFi) ==="
    # HiFi runs are single-end: fasterq-dump yields one <SRR>.fastq. --concatenate-reads
    # keeps any technical/biological split as one stream; then gzip.
    fasterq-dump --threads "${THREADS}" --outdir "${RAW_DIR}" \
        --concatenate-reads "${PREFETCH_DIR}/${SRR}/${SRR}.sra"
    if [ -s "${RAW_DIR}/${SRR}.fastq" ]; then
        gzip -f "${RAW_DIR}/${SRR}.fastq"
    else
        echo "WARNING: fasterq-dump produced no FASTQ for ${SRR}" >&2
    fi
    # prefetched .sra is large and no longer needed once the FASTQ exists
    [ -s "${FINAL}" ] && rm -rf "${PREFETCH_DIR:?}/${SRR}"
done

echo "Done. HiFi FASTQs in ${RAW_DIR}:"
ls -lh "${RAW_DIR}"/*.fastq.gz 2>/dev/null || echo "(none produced)"
