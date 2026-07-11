#!/bin/bash
#SBATCH --job-name=qc_u3analysis        # job name
#SBATCH --output=logs/slurm-%j.out      # output file
#SBATCH --error=logs/slurm-%j.err       # error log
#SBATCH --time=24:00:00                 # expected runtime
#SBATCH --ntasks=1                      # single multi-threaded task
#SBATCH --cpus-per-task=8               # number of threads for bwa/samtools/bcftools
#SBATCH --mem=32G                       # adjust based on genome size. (memory per thread x threads) + buffer. -> (1-2GB /t>
                                        # Safe  for viral genomes (1-2GB) X 8 -> 8-16GB + 16GB buffer

# Exit on any command failure within a pipeline
set -o pipefail

##-------DESCRIPTION--------##
#This script is for initial processing and quality control of HIV proviral sequences(longreads) for U3 region analysis.
#It performs all key steps, including, downloading reference genomes and raw sequencing reads, quality checks,
# aligning reads, genome assembly and v calling variants.

## NOTE
# Ensure this script is executed on a system with Anaconda, Miniconda or an equivalent Python distribution.
# Create an environment using the provided .yml file to ensure that all required channels and dependencies are installed /
# before running this script.

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

THREADS=${SLURM_CPUS_PER_TASK:-8}    # Use SLURM allocation or default to 8

# Directory structure
BASE_DIR="$(pwd)"
RAW_DIR="${BASE_DIR}/raw_data"
QC_DIR="${BASE_DIR}/qc_reports"
NANOPLOT_PRE_DIR="${QC_DIR}/nanoplot_pre"
NANOPLOT_POST_DIR="${QC_DIR}/nanoplot_post"
NANOQC_DIR="${QC_DIR}/nanoqc"
NANOSTAT_DIR="${QC_DIR}/nanostat"
MULTIQC_DIR="${QC_DIR}/multiqc"
FILTERED_DIR="${BASE_DIR}/filtered_data"
LOG_DIR="${BASE_DIR}/logs"

# SRA accessions for BioProject PRJNA765218 (NanoHIV - Oxford Nanopore GridION)

BIOPROJECT="PRJNA765218"

# 9 HIV-1 proviral genome samples from Stellenbosch University
SRR_ACCESSIONS=(
    SRR16005710    # 340116_D1P4
    SRR16005711    # 339606_P3G8
    SRR16005712    # 339606_P3G7
    SRR16005713    # 339606_P3D8
    SRR16005714    # 339266_P1C7
    SRR16005715    # 339266_C7P2
    SRR16005716    # 339266_P1C8
    SRR16005717    # 339266_D4P5
    SRR16005718    # 333716_P2D4
)

# Filtering thresholds (NanoFilt)
MIN_QUALITY=10         # Minimum average read quality (Phred score)
MIN_LENGTH=7000        # Minimum read length in bp

##==========================================================================##
##                         HELPER FUNCTIONS                                  ##
##==========================================================================##

log_msg() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

check_exit() {
    if [ $? -ne 0 ]; then
        log_msg "ERROR: $1"
        exit 1
    fi
}

##==========================================================================##
##               STEP 1: CREATE DIRECTORY STRUCTURE                          ##
##==========================================================================##

log_msg "========== STEP 1: Setting up directory structure =========="

mkdir -p "${RAW_DIR}" "${NANOPLOT_PRE_DIR}" "${NANOPLOT_POST_DIR}" \
         "${NANOQC_DIR}" "${NANOSTAT_DIR}" "${MULTIQC_DIR}" \
         "${FILTERED_DIR}" "${LOG_DIR}"

log_msg "Directory structure created under: ${BASE_DIR}"

##==========================================================================##
##               STEP 2: DOWNLOAD SRA DATA                                   ##
##==========================================================================##

log_msg "========== STEP 2: Downloading SRA data (${#SRR_ACCESSIONS[@]} samples) =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do
    log_msg "--- Processing ${SRR} ---"

    # Skip if FASTQ already exists
    if [ -f "${RAW_DIR}/${SRR}.fastq" ] || [ -f "${RAW_DIR}/${SRR}.fastq.gz" ]; then
        log_msg "FASTQ for ${SRR} already exists, skipping download."
        continue
    fi

    # Prefetch SRA file (with retry)
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

    # Convert SRA to FASTQ (single-end for Nanopore)
    log_msg "Converting ${SRR} to FASTQ..."
    fasterq-dump "${RAW_DIR}/${SRR}/${SRR}.sra" \
        --outdir "${RAW_DIR}" \
        --threads "${THREADS}" \
        --progress \
        2>&1 | tee "${LOG_DIR}/${SRR}_fasterq.log"
    check_exit "fasterq-dump failed for ${SRR}"

    # Compress FASTQ to save space
    log_msg "Compressing ${SRR}.fastq..."
    gzip -f "${RAW_DIR}/${SRR}.fastq"

    # Clean up SRA cache to save disk space
    rm -rf "${RAW_DIR}/${SRR}"

    log_msg "Completed download for ${SRR}"
done

log_msg "All SRA downloads completed."

##==========================================================================##
##               STEP 3: PRE-FILTERING QC (Raw Data Assessment)              ##
##==========================================================================##

log_msg "========== STEP 3: Running pre-filtering QC =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do
    FASTQ="${RAW_DIR}/${SRR}.fastq.gz"

    if [ ! -f "${FASTQ}" ]; then
        log_msg "WARNING: ${FASTQ} not found, skipping QC for ${SRR}."
        continue
    fi

    # --- NanoPlot: Comprehensive read length/quality plots ---
    log_msg "Running NanoPlot on ${SRR}..."
    NanoPlot \
        --fastq "${FASTQ}" \
        --outdir "${NANOPLOT_PRE_DIR}/${SRR}" \
        --prefix "${SRR}_pre_" \
        --threads "${THREADS}" \
        --loglength \
        --plots dot \
        --title "${SRR} - Pre-filtering QC" \
        2>&1 | tee "${LOG_DIR}/${SRR}_nanoplot_pre.log"
    check_exit "NanoPlot failed for ${SRR}"

    # --- NanoQC: Per-base quality across read positions ---
    log_msg "Running NanoQC on ${SRR}..."
    nanoQC \
        -o "${NANOQC_DIR}/${SRR}" \
        "${FASTQ}" \
        2>&1 | tee "${LOG_DIR}/${SRR}_nanoqc.log"
    check_exit "NanoQC failed for ${SRR}"

    # --- NanoStat: Quick text summary statistics ---
    log_msg "Running NanoStat on ${SRR}..."
    NanoStat \
        --fastq "${FASTQ}" \
        --outdir "${NANOSTAT_DIR}" \
        --name "${SRR}_pre_stats.txt" \
        --threads "${THREADS}" \
        2>&1 | tee "${LOG_DIR}/${SRR}_nanostat_pre.log"
    check_exit "NanoStat failed for ${SRR}"

    log_msg "Pre-filtering QC completed for ${SRR}"
done

##==========================================================================##
##               STEP 4: READ FILTERING & TRIMMING                           ##
##==========================================================================##

log_msg "========== STEP 4: Filtering and trimming reads =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do
    FASTQ="${RAW_DIR}/${SRR}.fastq.gz"
    TRIMMED="${FILTERED_DIR}/${SRR}_trimmed.fastq.gz"
    FILTERED="${FILTERED_DIR}/${SRR}_filtered.fastq.gz"

    if [ ! -f "${FASTQ}" ]; then
        log_msg "WARNING: ${FASTQ} not found, skipping filtering for ${SRR}."
        continue
    fi

    # Skip if filtered file already exists
    if [ -f "${FILTERED}" ]; then
        log_msg "Filtered FASTQ for ${SRR} already exists, skipping."
        continue
    fi

    # --- Porechop_ABI: Adapter trimming ---
    # Removes adapters and splits chimeric reads
    log_msg "Running Porechop_ABI on ${SRR}..."
    porechop_abi \
        --input "${FASTQ}" \
        --output "${TRIMMED}" \
        --threads "${THREADS}" \
        2>&1 | tee "${LOG_DIR}/${SRR}_porechop.log"

    # If porechop_abi is not available, try porechop
    if [ $? -ne 0 ]; then
        log_msg "Porechop_ABI not available, trying porechop..."
        porechop \
            --input "${FASTQ}" \
            --output "${TRIMMED}" \
            --threads "${THREADS}" \
            2>&1 | tee "${LOG_DIR}/${SRR}_porechop.log"

        # If porechop also fails, use raw file for filtering
        if [ $? -ne 0 ]; then
            log_msg "WARNING: Adapter trimming unavailable. Using raw reads for filtering."
            TRIMMED="${FASTQ}"
        fi
    fi

    # --- NanoFilt: Quality and length filtering ---
    log_msg "Running NanoFilt on ${SRR} (Q>=${MIN_QUALITY}, len>=${MIN_LENGTH})..."
    gunzip -c "${TRIMMED}" | \
        NanoFilt \
            --quality "${MIN_QUALITY}" \
            --length "${MIN_LENGTH}" | \
        gzip > "${FILTERED}"
    check_exit "NanoFilt failed for ${SRR}"

    # Clean up intermediate trimmed file (if different from raw)
    if [ "${TRIMMED}" != "${FASTQ}" ]; then
        rm -f "${TRIMMED}"
    fi

    # Report filtering stats
    RAW_READS=$(zcat "${FASTQ}" | awk 'END{print NR/4}')
    FILT_READS=$(zcat "${FILTERED}" | awk 'END{print NR/4}')
    RETAINED=$(echo "scale=1; ${FILT_READS}*100/${RAW_READS}" | bc)
    log_msg "${SRR}: ${RAW_READS} raw -> ${FILT_READS} filtered (${RETAINED}% retained)"

    log_msg "Filtering completed for ${SRR}"
done

##==========================================================================##
##               STEP 5: POST-FILTERING QC                                   ##
##==========================================================================##

log_msg "========== STEP 5: Running post-filtering QC =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do
    FILTERED="${FILTERED_DIR}/${SRR}_filtered.fastq.gz"

    if [ ! -f "${FILTERED}" ]; then
        log_msg "WARNING: ${FILTERED} not found, skipping post-QC for ${SRR}."
        continue
    fi

    # --- NanoPlot: Post-filtering assessment ---
    log_msg "Running NanoPlot (post-filter) on ${SRR}..."
    NanoPlot \
        --fastq "${FILTERED}" \
        --outdir "${NANOPLOT_POST_DIR}/${SRR}" \
        --prefix "${SRR}_post_" \
        --threads "${THREADS}" \
        --loglength \
        --plots dot \
        --title "${SRR} - Post-filtering QC" \
        2>&1 | tee "${LOG_DIR}/${SRR}_nanoplot_post.log"
    check_exit "NanoPlot (post-filter) failed for ${SRR}"

    # --- NanoStat: Post-filtering summary ---
    log_msg "Running NanoStat (post-filter) on ${SRR}..."
    NanoStat \
        --fastq "${FILTERED}" \
        --outdir "${NANOSTAT_DIR}" \
        --name "${SRR}_post_stats.txt" \
        --threads "${THREADS}" \
        2>&1 | tee "${LOG_DIR}/${SRR}_nanostat_post.log"
    check_exit "NanoStat (post-filter) failed for ${SRR}"

    log_msg "Post-filtering QC completed for ${SRR}"
done

##==========================================================================##
##               STEP 6: AGGREGATE QC REPORTS WITH MULTIQC                   ##
##==========================================================================##

log_msg "========== STEP 6: Aggregating QC reports with MultiQC =========="

multiqc \
    "${QC_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --filename "nanopore_qc_report" \
    --title "PRJNA765218 - Nanopore QC Summary" \
    --force \
    2>&1 | tee "${LOG_DIR}/multiqc.log"
check_exit "MultiQC failed"

log_msg "MultiQC report generated: ${MULTIQC_DIR}/nanopore_qc_report.html"

##==========================================================================##
##                          PIPELINE COMPLETE                                ##
##==========================================================================##

log_msg "=========================================="
log_msg "  RAW DATA QC PIPELINE COMPLETE"
log_msg "=========================================="
log_msg "BioProject:       ${BIOPROJECT}"
log_msg "Samples processed: ${#SRR_ACCESSIONS[@]}"
log_msg "Raw data:          ${RAW_DIR}"
log_msg "Filtered data:     ${FILTERED_DIR}"
log_msg "QC reports:        ${QC_DIR}"
log_msg "MultiQC summary:   ${MULTIQC_DIR}/nanopore_qc_report.html"
log_msg "=========================================="
