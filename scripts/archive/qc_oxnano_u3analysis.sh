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
set -o pipefail                                                    # fail the pipeline if any stage in a pipe fails

##-------DESCRIPTION--------##
#This script is for initial processing and quality control of HIV proviral sequences(longreads) for U3 region analysis.
#It performs all key steps, including, downloading reference genomes and raw sequencing reads, quality checks,
# aligning reads, genome assembly and v calling variants.

## NOTE
# Ensure this script is executed on a system with Anaconda, Miniconda or an equivalent Python distribution.
# Create an environment using the provided .yml file to ensure that all required channels and dependencies are installed /
# before running this script.

# Activate conda environment
CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"                 # usual location of the conda init script
if [ -f "$CONDA_SH" ]; then                                        # if that init script exists...
    source "$CONDA_SH"                                             # ...load it so `conda activate` works
elif command -v conda >/dev/null 2>&1; then                        # otherwise, if conda is already on PATH...
    source "$(conda info --base)/etc/profile.d/conda.sh"           # ...load the init script from conda's base dir
else                                                               # no conda available at all
    echo "ERROR: Conda not found. Please load conda before running this script." >&2  # tell the user on stderr
    exit 1                                                         # bail out since the QC tools live in the env
fi
conda activate HIV_U3analysis                                      # activate the env holding NanoPlot/NanoFilt/porechop/multiqc

##==========================================================================##
##                     CONFIGURATION & VARIABLES                             ##
##==========================================================================##

THREADS=${SLURM_CPUS_PER_TASK:-8}    # Use SLURM allocation or default to 8

# Directory structure
BASE_DIR="$(pwd)"                                                  # project root: assumes the script is launched from it
RAW_DIR="${BASE_DIR}/data/raw/oxnano"                             # where downloaded raw Nanopore FASTQs land
QC_DIR="${BASE_DIR}/results/reports/qc/oxnano"                    # top-level QC output dir
NANOPLOT_PRE_DIR="${QC_DIR}/nanoplot_pre"                        # NanoPlot reports on raw (pre-filter) reads
NANOPLOT_POST_DIR="${QC_DIR}/nanoplot_post"                      # NanoPlot reports on filtered (post-filter) reads
NANOQC_DIR="${QC_DIR}/nanoqc"                                     # nanoQC per-base quality reports
NANOSTAT_DIR="${QC_DIR}/nanostat"                                # NanoStat text summary stats
MULTIQC_DIR="${QC_DIR}/multiqc"                                  # aggregated MultiQC report
FILTERED_DIR="${BASE_DIR}/data/processed/oxnano/filtered"        # filtered/trimmed FASTQ output
LOG_DIR="${BASE_DIR}/logs"                                        # per-step tool logs

# SRA accessions for BioProject PRJNA765218 (NanoHIV - Oxford Nanopore GridION)

BIOPROJECT="PRJNA765218"                                          # NCBI BioProject these Nanopore samples come from

# 9 HIV-1 proviral genome samples from Stellenbosch University
SRR_ACCESSIONS=(                                                  # the 9 run accessions to download and QC
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

log_msg() {                                                       # timestamped progress logger
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"                      # print the message with a date/time stamp
}

check_exit() {                                                    # abort the pipeline if the previous command failed
    if [ $? -ne 0 ]; then                                         # inspect the last command's exit status
        log_msg "ERROR: $1"                                       # report the passed-in error context
        exit 1                                                    # stop the whole pipeline on failure
    fi
}

##==========================================================================##
##               STEP 1: CREATE DIRECTORY STRUCTURE                          ##
##==========================================================================##

log_msg "========== STEP 1: Setting up directory structure =========="  # announce the setup step

mkdir -p "${RAW_DIR}" "${NANOPLOT_PRE_DIR}" "${NANOPLOT_POST_DIR}" \
         "${NANOQC_DIR}" "${NANOSTAT_DIR}" "${MULTIQC_DIR}" \
         "${FILTERED_DIR}" "${LOG_DIR}"                              # create every output dir up front so later steps never fail

log_msg "Directory structure created under: ${BASE_DIR}"          # confirm setup done

##==========================================================================##
##               STEP 2: DOWNLOAD SRA DATA                                   ##
##==========================================================================##

log_msg "========== STEP 2: Downloading SRA data (${#SRR_ACCESSIONS[@]} samples) =========="  # announce download step + count

for SRR in "${SRR_ACCESSIONS[@]}"; do                             # download each accession one by one
    log_msg "--- Processing ${SRR} ---"                           # mark which sample we're on

    # Skip if FASTQ already exists
    if [ -f "${RAW_DIR}/${SRR}.fastq" ] || [ -f "${RAW_DIR}/${SRR}.fastq.gz" ]; then  # already downloaded (either form)?
        log_msg "FASTQ for ${SRR} already exists, skipping download."  # note the skip (makes reruns idempotent)
        continue                                                  # move on to the next accession
    fi

    # Prefetch SRA file (with retry)
    log_msg "Prefetching ${SRR}..."                               # progress marker
    prefetch "${SRR}" \
        --output-directory "${RAW_DIR}" \
        --max-size 50G \
        --progress \
        2>&1 | tee "${LOG_DIR}/${SRR}_prefetch.log"               # download the .sra into RAW_DIR, logging output
    check_exit "prefetch failed for ${SRR}"                       # stop if the download failed

    # Validate the downloaded SRA file
    log_msg "Validating ${SRR}..."                                # progress marker
    vdb-validate "${RAW_DIR}/${SRR}/${SRR}.sra" 2>&1 | tee "${LOG_DIR}/${SRR}_validate.log"  # check the .sra isn't corrupt
    if [ $? -ne 0 ]; then                                         # if validation reported a problem...
        log_msg "WARNING: Validation failed for ${SRR}, attempting re-download..."  # ...warn...
        rm -rf "${RAW_DIR}/${SRR}"                                # ...delete the bad copy...
        prefetch "${SRR}" --output-directory "${RAW_DIR}" --max-size 50G --force ALL  # ...and force a fresh download
        check_exit "Re-download failed for ${SRR}"                # give up if even the retry fails
    fi

    # Convert SRA to FASTQ (single-end for Nanopore)
    log_msg "Converting ${SRR} to FASTQ..."                       # progress marker
    fasterq-dump "${RAW_DIR}/${SRR}/${SRR}.sra" \
        --outdir "${RAW_DIR}" \
        --threads "${THREADS}" \
        --progress \
        2>&1 | tee "${LOG_DIR}/${SRR}_fasterq.log"                # extract a single FASTQ (Nanopore is single-end)
    check_exit "fasterq-dump failed for ${SRR}"                   # stop if conversion failed

    # Compress FASTQ to save space
    log_msg "Compressing ${SRR}.fastq..."                         # progress marker
    gzip -f "${RAW_DIR}/${SRR}.fastq"                             # gzip the FASTQ (downstream tools read .gz)

    # Clean up SRA cache to save disk space
    rm -rf "${RAW_DIR}/${SRR}"                                    # drop the bulky .sra now that the FASTQ exists

    log_msg "Completed download for ${SRR}"                       # per-sample done marker
done

log_msg "All SRA downloads completed."                            # whole download phase finished

##==========================================================================##
##               STEP 3: PRE-FILTERING QC (Raw Data Assessment)              ##
##==========================================================================##

log_msg "========== STEP 3: Running pre-filtering QC =========="  # announce pre-filter QC step

for SRR in "${SRR_ACCESSIONS[@]}"; do                             # QC each sample's raw reads
    FASTQ="${RAW_DIR}/${SRR}.fastq.gz"                            # path to the raw reads

    if [ ! -f "${FASTQ}" ]; then                                  # if the raw FASTQ is missing...
        log_msg "WARNING: ${FASTQ} not found, skipping QC for ${SRR}."  # ...warn...
        continue                                                  # ...and skip this sample
    fi

    # --- NanoPlot: Comprehensive read length/quality plots ---
    log_msg "Running NanoPlot on ${SRR}..."                       # progress marker
    NanoPlot \
        --fastq "${FASTQ}" \
        --outdir "${NANOPLOT_PRE_DIR}/${SRR}" \
        --prefix "${SRR}_pre_" \
        --threads "${THREADS}" \
        --loglength \
        --plots dot \
        --title "${SRR} - Pre-filtering QC" \
        2>&1 | tee "${LOG_DIR}/${SRR}_nanoplot_pre.log"           # length/quality plots for the raw reads
    check_exit "NanoPlot failed for ${SRR}"                       # stop on failure

    # --- NanoQC: Per-base quality across read positions ---
    log_msg "Running NanoQC on ${SRR}..."                         # progress marker
    nanoQC \
        -o "${NANOQC_DIR}/${SRR}" \
        "${FASTQ}" \
        2>&1 | tee "${LOG_DIR}/${SRR}_nanoqc.log"                 # per-base quality across read start/end positions
    check_exit "NanoQC failed for ${SRR}"                         # stop on failure

    # --- NanoStat: Quick text summary statistics ---
    log_msg "Running NanoStat on ${SRR}..."                       # progress marker
    NanoStat \
        --fastq "${FASTQ}" \
        --outdir "${NANOSTAT_DIR}" \
        --name "${SRR}_pre_stats.txt" \
        --threads "${THREADS}" \
        2>&1 | tee "${LOG_DIR}/${SRR}_nanostat_pre.log"           # quick text summary (N reads, N50, mean quality)
    check_exit "NanoStat failed for ${SRR}"                       # stop on failure

    log_msg "Pre-filtering QC completed for ${SRR}"               # per-sample done marker
done

##==========================================================================##
##               STEP 4: READ FILTERING & TRIMMING                           ##
##==========================================================================##

log_msg "========== STEP 4: Filtering and trimming reads =========="  # announce filtering step

for SRR in "${SRR_ACCESSIONS[@]}"; do                             # filter each sample
    FASTQ="${RAW_DIR}/${SRR}.fastq.gz"                            # raw reads input
    TRIMMED="${FILTERED_DIR}/${SRR}_trimmed.fastq.gz"            # intermediate adapter-trimmed reads
    FILTERED="${FILTERED_DIR}/${SRR}_filtered.fastq.gz"         # final quality/length-filtered reads

    if [ ! -f "${FASTQ}" ]; then                                  # if the raw FASTQ is missing...
        log_msg "WARNING: ${FASTQ} not found, skipping filtering for ${SRR}."  # ...warn...
        continue                                                  # ...and skip this sample
    fi

    # Skip if filtered file already exists
    if [ -f "${FILTERED}" ]; then                                 # already filtered on a previous run?
        log_msg "Filtered FASTQ for ${SRR} already exists, skipping."  # note the skip
        continue                                                  # move on
    fi

    # --- Porechop_ABI: Adapter trimming ---
    # Removes adapters and splits chimeric reads
    log_msg "Running Porechop_ABI on ${SRR}..."                   # progress marker
    porechop_abi \
        --input "${FASTQ}" \
        --output "${TRIMMED}" \
        --threads "${THREADS}" \
        2>&1 | tee "${LOG_DIR}/${SRR}_porechop.log"               # trim adapters/split chimeras (ABI infers adapters)

    # If porechop_abi is not available, try porechop
    if [ $? -ne 0 ]; then                                         # porechop_abi failed or isn't installed...
        log_msg "Porechop_ABI not available, trying porechop..."  # ...note the fallback...
        porechop \
            --input "${FASTQ}" \
            --output "${TRIMMED}" \
            --threads "${THREADS}" \
            2>&1 | tee "${LOG_DIR}/${SRR}_porechop.log"           # ...retry with plain porechop

        # If porechop also fails, use raw file for filtering
        if [ $? -ne 0 ]; then                                     # neither trimmer worked...
            log_msg "WARNING: Adapter trimming unavailable. Using raw reads for filtering."  # ...warn...
            TRIMMED="${FASTQ}"                                    # ...feed the raw reads straight into NanoFilt
        fi
    fi

    # --- NanoFilt: Quality and length filtering ---
    log_msg "Running NanoFilt on ${SRR} (Q>=${MIN_QUALITY}, len>=${MIN_LENGTH})..."  # progress marker with thresholds
    gunzip -c "${TRIMMED}" | \
        NanoFilt \
            --quality "${MIN_QUALITY}" \
            --length "${MIN_LENGTH}" | \
        gzip > "${FILTERED}"                                       # decompress, drop low-quality/short reads, recompress result
    check_exit "NanoFilt failed for ${SRR}"                       # stop on failure

    # Clean up intermediate trimmed file (if different from raw)
    if [ "${TRIMMED}" != "${FASTQ}" ]; then                       # only if we actually produced a trimmed file...
        rm -f "${TRIMMED}"                                        # ...delete it (the filtered file is what we keep)
    fi

    # Report filtering stats
    RAW_READS=$(zcat "${FASTQ}" | awk 'END{print NR/4}')          # count raw reads (4 FASTQ lines per read)
    FILT_READS=$(zcat "${FILTERED}" | awk 'END{print NR/4}')      # count reads surviving the filter
    RETAINED=$(echo "scale=1; ${FILT_READS}*100/${RAW_READS}" | bc)  # percent retained, to one decimal
    log_msg "${SRR}: ${RAW_READS} raw -> ${FILT_READS} filtered (${RETAINED}% retained)"  # log the before/after counts

    log_msg "Filtering completed for ${SRR}"                      # per-sample done marker
done

##==========================================================================##
##               STEP 5: POST-FILTERING QC                                   ##
##==========================================================================##

log_msg "========== STEP 5: Running post-filtering QC =========="  # announce post-filter QC step

for SRR in "${SRR_ACCESSIONS[@]}"; do                             # QC each sample's filtered reads
    FILTERED="${FILTERED_DIR}/${SRR}_filtered.fastq.gz"          # filtered reads input

    if [ ! -f "${FILTERED}" ]; then                               # if the filtered FASTQ is missing...
        log_msg "WARNING: ${FILTERED} not found, skipping post-QC for ${SRR}."  # ...warn...
        continue                                                  # ...and skip this sample
    fi

    # --- NanoPlot: Post-filtering assessment ---
    log_msg "Running NanoPlot (post-filter) on ${SRR}..."         # progress marker
    NanoPlot \
        --fastq "${FILTERED}" \
        --outdir "${NANOPLOT_POST_DIR}/${SRR}" \
        --prefix "${SRR}_post_" \
        --threads "${THREADS}" \
        --loglength \
        --plots dot \
        --title "${SRR} - Post-filtering QC" \
        2>&1 | tee "${LOG_DIR}/${SRR}_nanoplot_post.log"          # length/quality plots for the filtered reads (compare vs pre)
    check_exit "NanoPlot (post-filter) failed for ${SRR}"         # stop on failure

    # --- NanoStat: Post-filtering summary ---
    log_msg "Running NanoStat (post-filter) on ${SRR}..."         # progress marker
    NanoStat \
        --fastq "${FILTERED}" \
        --outdir "${NANOSTAT_DIR}" \
        --name "${SRR}_post_stats.txt" \
        --threads "${THREADS}" \
        2>&1 | tee "${LOG_DIR}/${SRR}_nanostat_post.log"          # text summary of the filtered reads
    check_exit "NanoStat (post-filter) failed for ${SRR}"         # stop on failure

    log_msg "Post-filtering QC completed for ${SRR}"              # per-sample done marker
done

##==========================================================================##
##               STEP 6: AGGREGATE QC REPORTS WITH MULTIQC                   ##
##==========================================================================##

log_msg "========== STEP 6: Aggregating QC reports with MultiQC =========="  # announce aggregation step

multiqc \
    "${QC_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --filename "nanopore_qc_report" \
    --title "PRJNA765218 - Nanopore QC Summary" \
    --force \
    2>&1 | tee "${LOG_DIR}/multiqc.log"                           # roll all NanoPlot/NanoStat reports into one HTML
check_exit "MultiQC failed"                                       # stop if aggregation failed

log_msg "MultiQC report generated: ${MULTIQC_DIR}/nanopore_qc_report.html"  # point user at the final report

##==========================================================================##
##                          PIPELINE COMPLETE                                ##
##==========================================================================##

log_msg "=========================================="              # final summary banner
log_msg "  RAW DATA QC PIPELINE COMPLETE"                         # completion headline
log_msg "=========================================="              # banner
log_msg "BioProject:       ${BIOPROJECT}"                         # which BioProject was processed
log_msg "Samples processed: ${#SRR_ACCESSIONS[@]}"                # how many samples went through
log_msg "Raw data:          ${RAW_DIR}"                           # where raw FASTQs live
log_msg "Filtered data:     ${FILTERED_DIR}"                      # where filtered FASTQs live
log_msg "QC reports:        ${QC_DIR}"                            # where QC reports live
log_msg "MultiQC summary:   ${MULTIQC_DIR}/nanopore_qc_report.html"  # the aggregated report path
log_msg "=========================================="              # banner
