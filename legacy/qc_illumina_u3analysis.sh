#!/bin/bash
#SBATCH --job-name=qc_illumina_u3
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --time=12:00:00                 # Estimated processing time (Illumina short reads process faster)
#SBATCH --ntasks=1                      # single multi-threaded task
#SBATCH --cpus-per-task=8               # threads for trimmomatic/fastqc
#SBATCH --mem=16G                       # 16GB is sufficient for short-read QC

# Exit on any command failure within a pipeline
set -o pipefail

##-------DESCRIPTION--------##
## This script downloads and performs quality control on Illumina MiSeq paired-end
## short-read HIV-1 whole genome sequences from BioProject PRJNA207834.
##
## STUDY: "Prevalence and Clinical Impacts of HIV-1 Intersubtype Recombinants
##         in Uganda Revealed by Near-Full-Genome Population and Deep Sequencing"
## ORGANISM: HIV-1 (subtypes A1, D, and A1-D recombinants)
## PLATFORM: Illumina MiSeq, 2x251bp paired-end
## SOURCE: Treatment-naive individuals from rural Mbarara, Uganda
## SAMPLES: 24 metagenomic WGS samples covering near-full-length HIV-1 genome
##          (including 5' LTR region)
##
## QC TOOLS: FastQC (per-base quality), Trimmomatic (adapter/quality trimming),
##           MultiQC (report aggregation)
##
## NOTE: Ensure this script is executed on a system with Anaconda, Miniconda or
## an equivalent Python distribution. Create an environment using the provided
## .yml file to ensure dependencies are installed before running this script.

# Author: Jovita Nagawa

# Activate conda environment
CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"
if [ -f "$CONDA_SH" ]; then
    source "$CONDA_SH"
elif command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
else
    echo "ERROR: Conda not found. Please load conda before running this script." >&2
    exit 1
fi
conda activate HIV_U3analysis

##==========================================================================##
##                     CONFIGURATION & VARIABLES                             ##
##==========================================================================##

BIOPROJECT="PRJNA207834"
THREADS=${SLURM_CPUS_PER_TASK:-8}    # Use SLURM allocation or default to 8

# SRA accessions for BioProject PRJNA207834
# 24 HIV-1 near-full-genome Illumina MiSeq paired-end samples from Uganda
# Subtypes: A1, D, and A1-D intersubtype recombinants (BSRI)
SRR_ACCESSIONS=(
    SRR908430    # AS03-00205
    SRR908431    # AS03-05969
    SRR908432    # AS04-01159
    SRR908433    # AS06-06468
    SRR908434    # AS06-10195
    SRR908435    # AS07-00787
    SRR908436    # AS08-00064
    SRR908437    # AS08-03339
    SRR908438    # AS10-10508
    SRR908439    # AS11-16494
    SRR908440    # AS12-08598
    SRR908441    # AS12-08878
    SRR908442    # MBA1005
    SRR908443    # MBA1089
    SRR908444    # MBA1120
    SRR908445    # MBA1218
    SRR908446    # MBA1256
    SRR908447    # MBA1465
    SRR908448    # MBA1470
    SRR908449    # MBA1478
    SRR908450    # MBA1516
    SRR908451    # MBA1548
    SRR908452    # MBA1549
    SRR908453    # MBA1581
)

# Trimmomatic parameters
TRIM_LEADING=3          # Cut bases from start if below quality
TRIM_TRAILING=3         # Cut bases from end if below quality
TRIM_SLIDINGWINDOW="4:20"  # Sliding window: window_size:quality_threshold
TRIM_MINLEN=50          # Minimum read length after trimming
TRIM_AVGQUAL=20         # Minimum average quality of the read

# Directory structure
BASE_DIR="$(pwd)"
RAW_DIR="${BASE_DIR}/illumina_raw_data"
QC_DIR="${BASE_DIR}/illumina_qc_reports"
FASTQC_PRE_DIR="${QC_DIR}/fastqc_pre"
FASTQC_POST_DIR="${QC_DIR}/fastqc_post"
MULTIQC_DIR="${QC_DIR}/multiqc"
TRIMMED_DIR="${BASE_DIR}/illumina_trimmed_data"
LOG_DIR="${BASE_DIR}/logs"

##==========================================================================##
##                         HELPER FUNCTIONS                                  ##
##==========================================================================##

log_msg() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ILLUMINA] $1"
}

check_exit() {
    if [ $? -ne 0 ]; then
        log_msg "ERROR: $1"
        exit 1
    fi
}

# Locate Trimmomatic adapter file
find_adapter_file() {
    # Search common locations for Trimmomatic adapter files
    local ADAPTER_LOCATIONS=(
        "${CONDA_PREFIX}/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
        "${CONDA_PREFIX}/share/trimmomatic-*/adapters/TruSeq3-PE-2.fa"
        "/usr/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
        "/usr/local/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
    )

    for pattern in "${ADAPTER_LOCATIONS[@]}"; do
        # Use compgen to expand globs safely
        local found
        found=$(compgen -G "${pattern}" 2>/dev/null | head -1)
        if [ -n "${found}" ] && [ -f "${found}" ]; then
            echo "${found}"
            return 0
        fi
    done

    log_msg "WARNING: TruSeq3-PE-2.fa adapter file not found. Skipping adapter trimming."
    echo ""
    return 1
}

##==========================================================================##
##               STEP 1: CREATE DIRECTORY STRUCTURE                          ##
##==========================================================================##

log_msg "========== STEP 1: Setting up directory structure =========="

mkdir -p "${RAW_DIR}" "${FASTQC_PRE_DIR}" "${FASTQC_POST_DIR}" \
         "${MULTIQC_DIR}" "${TRIMMED_DIR}" "${LOG_DIR}"

log_msg "Directory structure created under: ${BASE_DIR}"

##==========================================================================##
##               STEP 2: DOWNLOAD SRA DATA                                   ##
##==========================================================================##

log_msg "========== STEP 2: Downloading SRA data (${#SRR_ACCESSIONS[@]} samples) =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do
    log_msg "--- Processing ${SRR} ---"

    # Skip if paired-end FASTQs already exist
    if [ -f "${RAW_DIR}/${SRR}_1.fastq.gz" ] && [ -f "${RAW_DIR}/${SRR}_2.fastq.gz" ]; then
        log_msg "FASTQs for ${SRR} already exist, skipping download."
        continue
    fi

    # Prefetch SRA file
    log_msg "Prefetching ${SRR}..."
    prefetch "${SRR}" \
        --output-directory "${RAW_DIR}" \
        --max-size 50G \
        --progress \
        2>&1 | tee "${LOG_DIR}/${SRR}_prefetch.log"
    check_exit "prefetch failed for ${SRR}"

    # Validate the downloaded SRA file
    log_msg "Validating ${SRR}..."
    vdb-validate "${RAW_DIR}/${SRR}/${SRR}.sra" 2>&1 | tee "${LOG_DIR}/${SRR}_validate.log"
    if [ $? -ne 0 ]; then
        log_msg "WARNING: Validation failed for ${SRR}, attempting re-download..."
        rm -rf "${RAW_DIR}/${SRR}"
        prefetch "${SRR}" --output-directory "${RAW_DIR}" --max-size 50G --force ALL
        check_exit "Re-download failed for ${SRR}"
    fi

    # Convert SRA to paired-end FASTQ
    log_msg "Converting ${SRR} to paired-end FASTQs..."
    fasterq-dump "${RAW_DIR}/${SRR}/${SRR}.sra" \
        --outdir "${RAW_DIR}" \
        --split-3 \
        --threads "${THREADS}" \
        --progress \
        2>&1 | tee "${LOG_DIR}/${SRR}_fasterq.log"
    check_exit "fasterq-dump failed for ${SRR}"

    # Compress FASTQs to save space
    log_msg "Compressing FASTQs for ${SRR}..."
    gzip -f "${RAW_DIR}/${SRR}_1.fastq" 2>/dev/null
    gzip -f "${RAW_DIR}/${SRR}_2.fastq" 2>/dev/null
    # Also compress unpaired reads if they exist
    gzip -f "${RAW_DIR}/${SRR}.fastq" 2>/dev/null

    # Clean up SRA cache
    rm -rf "${RAW_DIR}/${SRR}"

    log_msg "Completed download for ${SRR}"
done

log_msg "All SRA downloads completed."

##==========================================================================##
##               STEP 3: PRE-TRIMMING QC WITH FASTQC                         ##
##==========================================================================##

log_msg "========== STEP 3: Running pre-trimming FastQC =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do
    R1="${RAW_DIR}/${SRR}_1.fastq.gz"
    R2="${RAW_DIR}/${SRR}_2.fastq.gz"

    if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then
        log_msg "WARNING: Paired FASTQs not found for ${SRR}, skipping pre-QC."
        continue
    fi

    log_msg "Running FastQC (pre-trimming) on ${SRR}..."
    fastqc \
        "${R1}" "${R2}" \
        --outdir "${FASTQC_PRE_DIR}" \
        --threads "${THREADS}" \
        --quiet \
        2>&1 | tee "${LOG_DIR}/${SRR}_fastqc_pre.log"
    check_exit "FastQC (pre-trimming) failed for ${SRR}"

    log_msg "Pre-trimming FastQC completed for ${SRR}"
done

##==========================================================================##
##               STEP 4: ADAPTER & QUALITY TRIMMING WITH TRIMMOMATIC         ##
##==========================================================================##

log_msg "========== STEP 4: Trimming reads with Trimmomatic =========="

# Find adapter file
ADAPTER_FILE=$(find_adapter_file)

for SRR in "${SRR_ACCESSIONS[@]}"; do
    R1="${RAW_DIR}/${SRR}_1.fastq.gz"
    R2="${RAW_DIR}/${SRR}_2.fastq.gz"

    # Output files
    R1_PAIRED="${TRIMMED_DIR}/${SRR}_1_paired.fastq.gz"
    R1_UNPAIRED="${TRIMMED_DIR}/${SRR}_1_unpaired.fastq.gz"
    R2_PAIRED="${TRIMMED_DIR}/${SRR}_2_paired.fastq.gz"
    R2_UNPAIRED="${TRIMMED_DIR}/${SRR}_2_unpaired.fastq.gz"

    if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then
        log_msg "WARNING: Paired FASTQs not found for ${SRR}, skipping trimming."
        continue
    fi

    # Skip if trimmed files already exist
    if [ -f "${R1_PAIRED}" ] && [ -f "${R2_PAIRED}" ]; then
        log_msg "Trimmed FASTQs for ${SRR} already exist, skipping."
        continue
    fi

    log_msg "Running Trimmomatic on ${SRR}..."

    # Build trimmomatic command with or without adapter trimming
    TRIM_STEPS=""
    if [ -n "${ADAPTER_FILE}" ] && [ -f "${ADAPTER_FILE}" ]; then
        TRIM_STEPS="ILLUMINACLIP:${ADAPTER_FILE}:2:30:10:2:True "
    fi
    TRIM_STEPS+="LEADING:${TRIM_LEADING} TRAILING:${TRIM_TRAILING} "
    TRIM_STEPS+="SLIDINGWINDOW:${TRIM_SLIDINGWINDOW} "
    TRIM_STEPS+="AVGQUAL:${TRIM_AVGQUAL} "
    TRIM_STEPS+="MINLEN:${TRIM_MINLEN}"

    trimmomatic PE \
        -threads "${THREADS}" \
        -phred33 \
        -summary "${LOG_DIR}/${SRR}_trimmomatic_summary.txt" \
        "${R1}" "${R2}" \
        "${R1_PAIRED}" "${R1_UNPAIRED}" \
        "${R2_PAIRED}" "${R2_UNPAIRED}" \
        ${TRIM_STEPS} \
        2>&1 | tee "${LOG_DIR}/${SRR}_trimmomatic.log"
    check_exit "Trimmomatic failed for ${SRR}"

    # Report trimming stats
    if [ -f "${LOG_DIR}/${SRR}_trimmomatic_summary.txt" ]; then
        log_msg "Trimmomatic summary for ${SRR}:"
        cat "${LOG_DIR}/${SRR}_trimmomatic_summary.txt" | while read -r line; do
            log_msg "  ${line}"
        done
    fi

    log_msg "Trimming completed for ${SRR}"
done

##==========================================================================##
##               STEP 5: POST-TRIMMING QC WITH FASTQC                        ##
##==========================================================================##

log_msg "========== STEP 5: Running post-trimming FastQC =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do
    R1_PAIRED="${TRIMMED_DIR}/${SRR}_1_paired.fastq.gz"
    R2_PAIRED="${TRIMMED_DIR}/${SRR}_2_paired.fastq.gz"

    if [ ! -f "${R1_PAIRED}" ] || [ ! -f "${R2_PAIRED}" ]; then
        log_msg "WARNING: Trimmed FASTQs not found for ${SRR}, skipping post-QC."
        continue
    fi

    log_msg "Running FastQC (post-trimming) on ${SRR}..."
    fastqc \
        "${R1_PAIRED}" "${R2_PAIRED}" \
        --outdir "${FASTQC_POST_DIR}" \
        --threads "${THREADS}" \
        --quiet \
        2>&1 | tee "${LOG_DIR}/${SRR}_fastqc_post.log"
    check_exit "FastQC (post-trimming) failed for ${SRR}"

    log_msg "Post-trimming FastQC completed for ${SRR}"
done

##==========================================================================##
##               STEP 6: AGGREGATE QC REPORTS WITH MULTIQC                   ##
##==========================================================================##

log_msg "========== STEP 6: Aggregating QC reports with MultiQC =========="

multiqc \
    "${QC_DIR}" "${LOG_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --filename "illumina_qc_report" \
    --title "PRJNA207834 - Illumina QC Summary (HIV-1 Uganda A1/D)" \
    --force \
    2>&1 | tee "${LOG_DIR}/multiqc_illumina.log"
check_exit "MultiQC failed"

log_msg "MultiQC report generated: ${MULTIQC_DIR}/illumina_qc_report.html"

##==========================================================================##
##                          PIPELINE COMPLETE                                ##
##==========================================================================##

log_msg "=========================================="
log_msg "  ILLUMINA RAW DATA QC PIPELINE COMPLETE"
log_msg "=========================================="
log_msg "BioProject:         ${BIOPROJECT}"
log_msg "Study:              HIV-1 intersubtype recombinants in Uganda"
log_msg "Subtypes:           A1, D, and A1-D recombinants"
log_msg "Platform:           Illumina MiSeq, 2x251bp paired-end"
log_msg "Samples processed:  ${#SRR_ACCESSIONS[@]}"
log_msg "Raw data:           ${RAW_DIR}"
log_msg "Trimmed data:       ${TRIMMED_DIR}"
log_msg "QC reports:         ${QC_DIR}"
log_msg "MultiQC summary:    ${MULTIQC_DIR}/illumina_qc_report.html"
log_msg "=========================================="
