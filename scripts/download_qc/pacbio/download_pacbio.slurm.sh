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
set -euo pipefail                                    # -e abort on any error, -u error on unset vars, pipefail catch failures anywhere in a pipe

CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"   # usual miniconda hook location
if [ -f "$CONDA_SH" ]; then source "$CONDA_SH"; else source "$(conda info --base)/etc/profile.d/conda.sh"; fi  # source it, or fall back to whatever conda is on PATH
conda activate HIV_U3analysis                        # activate the env with SRA-tools (prefetch/fasterq-dump)

cd /etc/ace-data/home/jnagawa/Internship             # work from the repo root so relative paths below resolve
source scripts/common/lib_compare.sh                 # load shared helpers (needed for subset_accessions)

RAW_DIR="data/raw/pacbio"                            # where the final HiFi FASTQs land
PREFETCH_DIR="data/raw/pacbio/.sra"                  # scratch dir for the large prefetched .sra files
mkdir -p "${RAW_DIR}" "${PREFETCH_DIR}"              # create both dirs if they don't exist
THREADS="${SLURM_CPUS_PER_TASK:-8}"                  # use Slurm's allocated CPUs, else default to 8

for SRR in $(subset_accessions pacbio scripts/common/subset_samples.tsv); do  # loop over just the PacBio accessions in the subset
    FINAL="${RAW_DIR}/${SRR}.fastq.gz"               # expected final output for this accession
    if [ -s "${FINAL}" ]; then                       # if it already exists (non-empty)...
        echo "=== ${SRR} already present (${FINAL}), skipping ==="  # ...report it...
        continue                                     # ...and skip the download (idempotent reruns)
    fi
    echo "=== prefetch ${SRR} ==="                   # progress marker
    # prefetch's default cap is 20G; set 100G explicitly so the subset (and a
    # Revio run, ~18GB, if added later) is never silently skipped.
    prefetch --max-size 100G --output-directory "${PREFETCH_DIR}" "${SRR}"  # download the .sra into the scratch dir with a raised size cap

    echo "=== fasterq-dump ${SRR} (single-end HiFi) ==="  # progress marker
    # HiFi runs are single-end: fasterq-dump yields one <SRR>.fastq. --concatenate-reads
    # keeps any technical/biological split as one stream; then gzip.
    fasterq-dump --threads "${THREADS}" --outdir "${RAW_DIR}" \
        --concatenate-reads "${PREFETCH_DIR}/${SRR}/${SRR}.sra"  # convert the .sra to a single FASTQ in RAW_DIR
    if [ -s "${RAW_DIR}/${SRR}.fastq" ]; then         # if the conversion produced a non-empty FASTQ...
        gzip -f "${RAW_DIR}/${SRR}.fastq"            # ...compress it in place (-f overwrites any stale .gz)
    else
        echo "WARNING: fasterq-dump produced no FASTQ for ${SRR}" >&2  # otherwise warn (something went wrong)
    fi
    # prefetched .sra is large and no longer needed once the FASTQ exists
    [ -s "${FINAL}" ] && rm -rf "${PREFETCH_DIR:?}/${SRR}"  # only delete the .sra once the gzipped FASTQ exists; :? guards against an empty var wiping the wrong path
done

echo "Done. HiFi FASTQs in ${RAW_DIR}:"              # final summary header
ls -lh "${RAW_DIR}"/*.fastq.gz 2>/dev/null || echo "(none produced)"  # list the produced FASTQs, or note if none exist
