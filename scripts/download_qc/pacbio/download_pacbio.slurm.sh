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
# -e abort on any error, -u error on unset vars, pipefail catch failures anywhere in a pipe
set -euo pipefail

CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"   # usual miniconda hook location
# source it, or fall back to whatever conda is on PATH
if [ -f "$CONDA_SH" ]; then source "$CONDA_SH"; else source "$(conda info --base)/etc/profile.d/conda.sh"; fi
# activate the env with SRA-tools (prefetch/fasterq-dump)
conda activate HIV_U3analysis

# work from the repo root so relative paths below resolve
cd /etc/ace-data/home/jnagawa/Internship
# load shared helpers (needed for subset_accessions)
source scripts/common/lib_compare.sh

RAW_DIR="data/raw/pacbio"                            # where the final HiFi FASTQs land
# scratch dir for the large prefetched .sra files
PREFETCH_DIR="data/raw/pacbio/.sra"
mkdir -p "${RAW_DIR}" "${PREFETCH_DIR}"              # create both dirs if they don't exist
THREADS="${SLURM_CPUS_PER_TASK:-8}"                  # use Slurm's allocated CPUs, else default to 8

# loop over the `pacbio_sra` rows -- the SRA HiFi runs. (`pacbio` is the masked
# local subset, which has no SRA accession to prefetch.)
for SRR in $(subset_accessions pacbio_sra scripts/common/subset_samples.tsv); do
    FINAL="${RAW_DIR}/${SRR}.fastq.gz"               # expected final output for this accession
    if [ -s "${FINAL}" ]; then                       # if it already exists (non-empty)...
        echo "=== ${SRR} already present (${FINAL}), skipping ==="  # ...report it...
        continue                                     # ...and skip the download (idempotent reruns)
    fi
    echo "=== prefetch ${SRR} ==="                   # progress marker
    # prefetch's default cap is 20G; set 100G explicitly so the subset (and a
    # Revio run, ~18GB, if added later) is never silently skipped.
    # download the .sra into the scratch dir with a raised size cap
    prefetch --max-size 100G --output-directory "${PREFETCH_DIR}" "${SRR}"

    echo "=== fasterq-dump ${SRR} (single-end HiFi) ==="  # progress marker
    # HiFi runs are single-end: fasterq-dump yields one <SRR>.fastq. --concatenate-reads
    # keeps any technical/biological split as one stream; then gzip.
    # convert the .sra to a single FASTQ in RAW_DIR
    fasterq-dump --threads "${THREADS}" --outdir "${RAW_DIR}" \
        --concatenate-reads "${PREFETCH_DIR}/${SRR}/${SRR}.sra"
    # if the conversion produced a non-empty FASTQ...
    if [ -s "${RAW_DIR}/${SRR}.fastq" ]; then
        # ...compress it in place (-f overwrites any stale .gz)
        gzip -f "${RAW_DIR}/${SRR}.fastq"
    else
        # otherwise warn (something went wrong)
        echo "WARNING: fasterq-dump produced no FASTQ for ${SRR}" >&2
    fi
    # prefetched .sra is large and no longer needed once the FASTQ exists
    # only delete the .sra once the gzipped FASTQ exists; :? guards against an empty var wiping the
    # wrong path
    [ -s "${FINAL}" ] && rm -rf "${PREFETCH_DIR:?}/${SRR}"
done

echo "Done. HiFi FASTQs in ${RAW_DIR}:"              # final summary header
# list the produced FASTQs, or note if none exist
ls -lh "${RAW_DIR}"/*.fastq.gz 2>/dev/null || echo "(none produced)"
