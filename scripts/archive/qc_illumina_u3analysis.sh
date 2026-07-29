#!/bin/bash
#SBATCH --job-name=qc_illumina_u3       # SLURM job name shown in the queue
#SBATCH --output=logs/slurm-%j.out      # stdout log path (%j expands to the job ID)
#SBATCH --error=logs/slurm-%j.err       # stderr log path (%j expands to the job ID)
#SBATCH --time=12:00:00                 # Estimated processing time (Illumina short reads process faster)
#SBATCH --ntasks=1                      # single multi-threaded task
#SBATCH --cpus-per-task=8               # threads for trimmomatic/fastqc
#SBATCH --mem=16G                       # 16GB is sufficient for short-read QC

# Exit on any command failure within a pipeline
# fail the pipeline if any stage in a pipe fails
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
# usual location of the conda init script
CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"
if [ -f "$CONDA_SH" ]; then                                        # if that init script exists...
    # ...load it so `conda activate` works
    source "$CONDA_SH"
# otherwise, if conda is already on PATH...
elif command -v conda >/dev/null 2>&1; then
    # ...load the init script from conda's own base dir
    source "$(conda info --base)/etc/profile.d/conda.sh"
else                                                               # no conda available at all
    # tell the user on stderr
    echo "ERROR: Conda not found. Please load conda before running this script." >&2
    # bail out since the tools live in the env
    exit 1
fi
# activate the env holding fastqc/trimmomatic/multiqc
conda activate HIV_U3analysis

##==========================================================================##
##                     CONFIGURATION & VARIABLES                             ##
##==========================================================================##

# NCBI BioProject these Illumina samples come from
BIOPROJECT="PRJNA207834"
THREADS=${SLURM_CPUS_PER_TASK:-8}    # Use SLURM allocation or default to 8

# SRA accessions for BioProject PRJNA207834
# 24 HIV-1 near-full-genome Illumina MiSeq paired-end samples from Uganda
# Subtypes: A1, D, and A1-D intersubtype recombinants (BSRI)
# the 24 run accessions to download and QC
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
# project root: assumes the script is launched from it
BASE_DIR="$(pwd)"
# where downloaded raw FASTQs land
RAW_DIR="${BASE_DIR}/data/raw/illumina"
QC_DIR="${BASE_DIR}/results/reports/qc/illumina"                   # top-level QC output dir
# FastQC reports on raw (pre-trim) reads
FASTQC_PRE_DIR="${QC_DIR}/fastqc_pre"
# FastQC reports on trimmed (post-trim) reads
FASTQC_POST_DIR="${QC_DIR}/fastqc_post"
MULTIQC_DIR="${QC_DIR}/multiqc"                                    # aggregated MultiQC report
TRIMMED_DIR="${BASE_DIR}/data/processed/illumina/trimmed"          # trimmed FASTQ output
LOG_DIR="${BASE_DIR}/logs"                                         # per-step tool logs

##==========================================================================##
##                         HELPER FUNCTIONS                                  ##
##==========================================================================##

# timestamped progress logger prefixed with [ILLUMINA]
log_msg() {
    # print the message with a date/time stamp
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ILLUMINA] $1"
}

# abort the pipeline if the previous command failed
check_exit() {
    # inspect the last command's exit status
    if [ $? -ne 0 ]; then
        # report the passed-in error context
        log_msg "ERROR: $1"
        # stop the whole pipeline on failure
        exit 1
    fi
}

# Locate Trimmomatic adapter file
# search known install paths for the TruSeq3 adapter FASTA
find_adapter_file() {
    # Search common locations for Trimmomatic adapter files
    # candidate paths where the adapter file may live
    local ADAPTER_LOCATIONS=(
        "${CONDA_PREFIX}/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
        "${CONDA_PREFIX}/share/trimmomatic-*/adapters/TruSeq3-PE-2.fa"
        "/usr/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
        "/usr/local/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
    )

    for pattern in "${ADAPTER_LOCATIONS[@]}"; do                   # try each candidate path in turn
        # Use compgen to expand globs safely
        # will hold the first matching real file
        local found
        # glob-expand the pattern, take the first hit
        found=$(compgen -G "${pattern}" 2>/dev/null | head -1)
        if [ -n "${found}" ] && [ -f "${found}" ]; then            # if a real file was found...
            # ...emit its path (function's return value)
            echo "${found}"
            # ...and stop searching, signalling success
            return 0
        fi
    done

    # none found: warn and continue
    log_msg "WARNING: TruSeq3-PE-2.fa adapter file not found. Skipping adapter trimming."
    # emit empty string so callers know there's no adapter
    echo ""
    # signal "not found" to the caller
    return 1
}

##==========================================================================##
##               STEP 1: CREATE DIRECTORY STRUCTURE                          ##
##==========================================================================##

log_msg "========== STEP 1: Setting up directory structure =========="  # announce the setup step

# create every output dir up front so later steps never fail
mkdir -p "${RAW_DIR}" "${FASTQC_PRE_DIR}" "${FASTQC_POST_DIR}" \
         "${MULTIQC_DIR}" "${TRIMMED_DIR}" "${LOG_DIR}"

log_msg "Directory structure created under: ${BASE_DIR}"           # confirm setup done

##==========================================================================##
##               STEP 2: DOWNLOAD SRA DATA                                   ##
##==========================================================================##

# announce download step + count
log_msg "========== STEP 2: Downloading SRA data (${#SRR_ACCESSIONS[@]} samples) =========="

# download each accession one by one
for SRR in "${SRR_ACCESSIONS[@]}"; do
    log_msg "--- Processing ${SRR} ---"                            # mark which sample we're on

    # Skip if paired-end FASTQs already exist
    # both mates already present?
    if [ -f "${RAW_DIR}/${SRR}_1.fastq.gz" ] && [ -f "${RAW_DIR}/${SRR}_2.fastq.gz" ]; then
        # note the skip (makes reruns idempotent)
        log_msg "FASTQs for ${SRR} already exist, skipping download."
        continue                                                   # move on to the next accession
    fi

    # Prefetch SRA file
    log_msg "Prefetching ${SRR}..."                                # progress marker
    # download the .sra into RAW_DIR, logging output
    prefetch "${SRR}" \
        --output-directory "${RAW_DIR}" \
        --max-size 50G \
        --progress \
        2>&1 | tee "${LOG_DIR}/${SRR}_prefetch.log"
    check_exit "prefetch failed for ${SRR}"                        # stop if the download failed

    # Validate the downloaded SRA file
    log_msg "Validating ${SRR}..."                                 # progress marker
    # check the .sra isn't corrupt
    vdb-validate "${RAW_DIR}/${SRR}/${SRR}.sra" 2>&1 | tee "${LOG_DIR}/${SRR}_validate.log"
    # if validation reported a problem...
    if [ $? -ne 0 ]; then
        log_msg "WARNING: Validation failed for ${SRR}, attempting re-download..."  # ...warn...
        rm -rf "${RAW_DIR}/${SRR}"                                 # ...delete the bad copy...
        # ...and force a fresh download
        prefetch "${SRR}" --output-directory "${RAW_DIR}" --max-size 50G --force ALL
        check_exit "Re-download failed for ${SRR}"                 # give up if even the retry fails
    fi

    # Convert SRA to paired-end FASTQ
    log_msg "Converting ${SRR} to paired-end FASTQs..."            # progress marker
    # --split-3 writes _1/_2 (and singletons) FASTQs
    fasterq-dump "${RAW_DIR}/${SRR}/${SRR}.sra" \
        --outdir "${RAW_DIR}" \
        --split-3 \
        --threads "${THREADS}" \
        --progress \
        2>&1 | tee "${LOG_DIR}/${SRR}_fasterq.log"
    check_exit "fasterq-dump failed for ${SRR}"                    # stop if conversion failed

    # Compress FASTQs to save space
    log_msg "Compressing FASTQs for ${SRR}..."                     # progress marker
    # gzip mate 1 (downstream tools read .gz)
    gzip -f "${RAW_DIR}/${SRR}_1.fastq" 2>/dev/null
    gzip -f "${RAW_DIR}/${SRR}_2.fastq" 2>/dev/null                # gzip mate 2
    # Also compress unpaired reads if they exist
    # gzip singleton reads if any were produced
    gzip -f "${RAW_DIR}/${SRR}.fastq" 2>/dev/null

    # Clean up SRA cache
    # drop the bulky .sra now that FASTQs exist
    rm -rf "${RAW_DIR}/${SRR}"

    log_msg "Completed download for ${SRR}"                        # per-sample done marker
done

log_msg "All SRA downloads completed."                             # whole download phase finished

##==========================================================================##
##               STEP 3: PRE-TRIMMING QC WITH FASTQC                         ##
##==========================================================================##

log_msg "========== STEP 3: Running pre-trimming FastQC =========="  # announce pre-trim QC step

for SRR in "${SRR_ACCESSIONS[@]}"; do                              # QC each sample's raw reads
    R1="${RAW_DIR}/${SRR}_1.fastq.gz"                              # path to raw mate 1
    R2="${RAW_DIR}/${SRR}_2.fastq.gz"                              # path to raw mate 2

    if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then                 # if either mate is missing...
        log_msg "WARNING: Paired FASTQs not found for ${SRR}, skipping pre-QC."  # ...warn...
        continue                                                   # ...and skip this sample
    fi

    log_msg "Running FastQC (pre-trimming) on ${SRR}..."           # progress marker
    # per-base quality report on the raw reads
    fastqc \
        "${R1}" "${R2}" \
        --outdir "${FASTQC_PRE_DIR}" \
        --threads "${THREADS}" \
        --quiet \
        2>&1 | tee "${LOG_DIR}/${SRR}_fastqc_pre.log"
    check_exit "FastQC (pre-trimming) failed for ${SRR}"           # stop on failure

    log_msg "Pre-trimming FastQC completed for ${SRR}"             # per-sample done marker
done

##==========================================================================##
##               STEP 4: ADAPTER & QUALITY TRIMMING WITH TRIMMOMATIC         ##
##==========================================================================##

log_msg "========== STEP 4: Trimming reads with Trimmomatic =========="  # announce trimming step

# Find adapter file
# locate the adapter FASTA once, reuse for all samples
ADAPTER_FILE=$(find_adapter_file)

for SRR in "${SRR_ACCESSIONS[@]}"; do                              # trim each sample
    R1="${RAW_DIR}/${SRR}_1.fastq.gz"                              # raw mate 1 input
    R2="${RAW_DIR}/${SRR}_2.fastq.gz"                              # raw mate 2 input

    # Output files
    # mate 1 reads whose partner also survived
    R1_PAIRED="${TRIMMED_DIR}/${SRR}_1_paired.fastq.gz"
    # mate 1 reads whose partner was dropped
    R1_UNPAIRED="${TRIMMED_DIR}/${SRR}_1_unpaired.fastq.gz"
    # mate 2 reads whose partner also survived
    R2_PAIRED="${TRIMMED_DIR}/${SRR}_2_paired.fastq.gz"
    # mate 2 reads whose partner was dropped
    R2_UNPAIRED="${TRIMMED_DIR}/${SRR}_2_unpaired.fastq.gz"

    if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then                 # if either raw mate is missing...
        log_msg "WARNING: Paired FASTQs not found for ${SRR}, skipping trimming."  # ...warn...
        continue                                                   # ...and skip this sample
    fi

    # Skip if trimmed files already exist
    # already trimmed on a previous run?
    if [ -f "${R1_PAIRED}" ] && [ -f "${R2_PAIRED}" ]; then
        log_msg "Trimmed FASTQs for ${SRR} already exist, skipping."  # note the skip
        continue                                                   # move on
    fi

    log_msg "Running Trimmomatic on ${SRR}..."                     # progress marker

    # Build trimmomatic command with or without adapter trimming
    # accumulate the Trimmomatic operation list
    TRIM_STEPS=""
    # only add adapter clipping if we found the FASTA
    if [ -n "${ADAPTER_FILE}" ] && [ -f "${ADAPTER_FILE}" ]; then
        # clip Illumina adapters (seed/palindrome/simple thresholds)
        TRIM_STEPS="ILLUMINACLIP:${ADAPTER_FILE}:2:30:10:2:True "
    fi
    # trim low-quality bases off both read ends
    TRIM_STEPS+="LEADING:${TRIM_LEADING} TRAILING:${TRIM_TRAILING} "
    TRIM_STEPS+="SLIDINGWINDOW:${TRIM_SLIDINGWINDOW} "             # sliding-window quality trimming
    # drop reads below the mean-quality cutoff
    TRIM_STEPS+="AVGQUAL:${TRIM_AVGQUAL} "
    # drop reads shorter than the minimum length
    TRIM_STEPS+="MINLEN:${TRIM_MINLEN}"

    # run paired-end trimming with the assembled step list
    trimmomatic PE \
        -threads "${THREADS}" \
        -phred33 \
        -summary "${LOG_DIR}/${SRR}_trimmomatic_summary.txt" \
        "${R1}" "${R2}" \
        "${R1_PAIRED}" "${R1_UNPAIRED}" \
        "${R2_PAIRED}" "${R2_UNPAIRED}" \
        ${TRIM_STEPS} \
        2>&1 | tee "${LOG_DIR}/${SRR}_trimmomatic.log"
    check_exit "Trimmomatic failed for ${SRR}"                     # stop on failure

    # Report trimming stats
    # if a summary file was written...
    if [ -f "${LOG_DIR}/${SRR}_trimmomatic_summary.txt" ]; then
        log_msg "Trimmomatic summary for ${SRR}:"                  # ...header it in the log...
        # ...read it line by line...
        cat "${LOG_DIR}/${SRR}_trimmomatic_summary.txt" | while read -r line; do
            # ...and echo each stat into the main log
            log_msg "  ${line}"
        done
    fi

    log_msg "Trimming completed for ${SRR}"                        # per-sample done marker
done

##==========================================================================##
##               STEP 5: POST-TRIMMING QC WITH FASTQC                        ##
##==========================================================================##

log_msg "========== STEP 5: Running post-trimming FastQC =========="  # announce post-trim QC step

for SRR in "${SRR_ACCESSIONS[@]}"; do                              # QC each sample's trimmed reads
    R1_PAIRED="${TRIMMED_DIR}/${SRR}_1_paired.fastq.gz"            # trimmed mate 1 (paired) input
    R2_PAIRED="${TRIMMED_DIR}/${SRR}_2_paired.fastq.gz"            # trimmed mate 2 (paired) input

    if [ ! -f "${R1_PAIRED}" ] || [ ! -f "${R2_PAIRED}" ]; then    # if a trimmed mate is missing...
        log_msg "WARNING: Trimmed FASTQs not found for ${SRR}, skipping post-QC."  # ...warn...
        continue                                                   # ...and skip this sample
    fi

    log_msg "Running FastQC (post-trimming) on ${SRR}..."          # progress marker
    # quality report on the trimmed reads (compare vs pre)
    fastqc \
        "${R1_PAIRED}" "${R2_PAIRED}" \
        --outdir "${FASTQC_POST_DIR}" \
        --threads "${THREADS}" \
        --quiet \
        2>&1 | tee "${LOG_DIR}/${SRR}_fastqc_post.log"
    check_exit "FastQC (post-trimming) failed for ${SRR}"          # stop on failure

    log_msg "Post-trimming FastQC completed for ${SRR}"            # per-sample done marker
done

##==========================================================================##
##               STEP 6: AGGREGATE QC REPORTS WITH MULTIQC                   ##
##==========================================================================##

# announce aggregation step
log_msg "========== STEP 6: Aggregating QC reports with MultiQC =========="

# roll every FastQC/Trimmomatic report into one HTML
multiqc \
    "${QC_DIR}" "${LOG_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --filename "illumina_qc_report" \
    --title "PRJNA207834 - Illumina QC Summary (HIV-1 Uganda A1/D)" \
    --force \
    2>&1 | tee "${LOG_DIR}/multiqc_illumina.log"
check_exit "MultiQC failed"                                        # stop if aggregation failed

# point user at the final report
log_msg "MultiQC report generated: ${MULTIQC_DIR}/illumina_qc_report.html"

##==========================================================================##
##                          PIPELINE COMPLETE                                ##
##==========================================================================##

log_msg "=========================================="                         # final summary banner
log_msg "  ILLUMINA RAW DATA QC PIPELINE COMPLETE"                            # completion headline
log_msg "=========================================="                         # banner
# which BioProject was processed
log_msg "BioProject:         ${BIOPROJECT}"
log_msg "Study:              HIV-1 intersubtype recombinants in Uganda"       # study context
log_msg "Subtypes:           A1, D, and A1-D recombinants"                    # subtypes covered
log_msg "Platform:           Illumina MiSeq, 2x251bp paired-end"             # sequencing platform
# how many samples went through
log_msg "Samples processed:  ${#SRR_ACCESSIONS[@]}"
# where raw FASTQs live
log_msg "Raw data:           ${RAW_DIR}"
# where trimmed FASTQs live
log_msg "Trimmed data:       ${TRIMMED_DIR}"
# where QC reports live
log_msg "QC reports:         ${QC_DIR}"
# the aggregated report path
log_msg "MultiQC summary:    ${MULTIQC_DIR}/illumina_qc_report.html"
log_msg "=========================================="                         # banner
