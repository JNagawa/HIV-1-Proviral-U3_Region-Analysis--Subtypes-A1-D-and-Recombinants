#!/bin/bash
#SBATCH --job-name=qc_illumina_u3       # SLURM job name shown in the queue
#SBATCH --output=logs/slurm-%j.out      # stdout log path (%j expands to the job ID)
#SBATCH --error=logs/slurm-%j.err       # stderr log path (%j expands to the job ID)
#SBATCH --time=12:00:00                 # Estimated processing time (Illumina short reads process faster)
#SBATCH --ntasks=1                      # single multi-threaded task
#SBATCH --cpus-per-task=8               # threads for trimmomatic/fastqc
#SBATCH --mem=16G                       # 16GB is sufficient for short-read QC

# Exit on any command failure within a pipeline
set -o pipefail                                                    # fail the pipeline if any stage in a pipe fails

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
CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"                 # usual location of the conda init script
if [ -f "$CONDA_SH" ]; then                                        # if that init script exists...
    source "$CONDA_SH"                                             # ...load it so `conda activate` works
elif command -v conda >/dev/null 2>&1; then                        # otherwise, if conda is already on PATH...
    source "$(conda info --base)/etc/profile.d/conda.sh"           # ...load the init script from conda's own base dir
else                                                               # no conda available at all
    echo "ERROR: Conda not found. Please load conda before running this script." >&2  # tell the user on stderr
    exit 1                                                         # bail out since the tools live in the env
fi
conda activate HIV_U3analysis                                      # activate the env holding fastqc/trimmomatic/multiqc

##==========================================================================##
##                     CONFIGURATION & VARIABLES                             ##
##==========================================================================##

BIOPROJECT="PRJNA207834"                                           # NCBI BioProject these Illumina samples come from
THREADS=${SLURM_CPUS_PER_TASK:-8}    # Use SLURM allocation or default to 8

# SRA accessions for BioProject PRJNA207834
# 24 HIV-1 near-full-genome Illumina MiSeq paired-end samples from Uganda
# Subtypes: A1, D, and A1-D intersubtype recombinants (BSRI)
SRR_ACCESSIONS=(                                                   # the 24 run accessions to download and QC
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
BASE_DIR="$(pwd)"                                                  # project root: assumes the script is launched from it
RAW_DIR="${BASE_DIR}/data/raw/illumina"                            # where downloaded raw FASTQs land
QC_DIR="${BASE_DIR}/results/reports/qc/illumina"                   # top-level QC output dir
FASTQC_PRE_DIR="${QC_DIR}/fastqc_pre"                              # FastQC reports on raw (pre-trim) reads
FASTQC_POST_DIR="${QC_DIR}/fastqc_post"                            # FastQC reports on trimmed (post-trim) reads
MULTIQC_DIR="${QC_DIR}/multiqc"                                    # aggregated MultiQC report
TRIMMED_DIR="${BASE_DIR}/data/processed/illumina/trimmed"          # trimmed FASTQ output
LOG_DIR="${BASE_DIR}/logs"                                         # per-step tool logs

##==========================================================================##
##                         HELPER FUNCTIONS                                  ##
##==========================================================================##

log_msg() {                                                        # timestamped progress logger prefixed with [ILLUMINA]
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ILLUMINA] $1"            # print the message with a date/time stamp
}

check_exit() {                                                     # abort the pipeline if the previous command failed
    if [ $? -ne 0 ]; then                                          # inspect the last command's exit status
        log_msg "ERROR: $1"                                        # report the passed-in error context
        exit 1                                                     # stop the whole pipeline on failure
    fi
}

# Locate Trimmomatic adapter file
find_adapter_file() {                                              # search known install paths for the TruSeq3 adapter FASTA
    # Search common locations for Trimmomatic adapter files
    local ADAPTER_LOCATIONS=(                                      # candidate paths where the adapter file may live
        "${CONDA_PREFIX}/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
        "${CONDA_PREFIX}/share/trimmomatic-*/adapters/TruSeq3-PE-2.fa"
        "/usr/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
        "/usr/local/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
    )

    for pattern in "${ADAPTER_LOCATIONS[@]}"; do                   # try each candidate path in turn
        # Use compgen to expand globs safely
        local found                                                # will hold the first matching real file
        found=$(compgen -G "${pattern}" 2>/dev/null | head -1)     # glob-expand the pattern, take the first hit
        if [ -n "${found}" ] && [ -f "${found}" ]; then            # if a real file was found...
            echo "${found}"                                        # ...emit its path (function's return value)
            return 0                                               # ...and stop searching, signalling success
        fi
    done

    log_msg "WARNING: TruSeq3-PE-2.fa adapter file not found. Skipping adapter trimming."  # none found: warn and continue
    echo ""                                                        # emit empty string so callers know there's no adapter
    return 1                                                       # signal "not found" to the caller
}

##==========================================================================##
##               STEP 1: CREATE DIRECTORY STRUCTURE                          ##
##==========================================================================##

log_msg "========== STEP 1: Setting up directory structure =========="  # announce the setup step

mkdir -p "${RAW_DIR}" "${FASTQC_PRE_DIR}" "${FASTQC_POST_DIR}" \
         "${MULTIQC_DIR}" "${TRIMMED_DIR}" "${LOG_DIR}"             # create every output dir up front so later steps never fail

log_msg "Directory structure created under: ${BASE_DIR}"           # confirm setup done

##==========================================================================##
##               STEP 2: DOWNLOAD SRA DATA                                   ##
##==========================================================================##

log_msg "========== STEP 2: Downloading SRA data (${#SRR_ACCESSIONS[@]} samples) =========="  # announce download step + count

for SRR in "${SRR_ACCESSIONS[@]}"; do                              # download each accession one by one
    log_msg "--- Processing ${SRR} ---"                            # mark which sample we're on

    # Skip if paired-end FASTQs already exist
    if [ -f "${RAW_DIR}/${SRR}_1.fastq.gz" ] && [ -f "${RAW_DIR}/${SRR}_2.fastq.gz" ]; then  # both mates already present?
        log_msg "FASTQs for ${SRR} already exist, skipping download."  # note the skip (makes reruns idempotent)
        continue                                                   # move on to the next accession
    fi

    # Prefetch SRA file
    log_msg "Prefetching ${SRR}..."                                # progress marker
    prefetch "${SRR}" \
        --output-directory "${RAW_DIR}" \
        --max-size 50G \
        --progress \
        2>&1 | tee "${LOG_DIR}/${SRR}_prefetch.log"                # download the .sra into RAW_DIR, logging output
    check_exit "prefetch failed for ${SRR}"                        # stop if the download failed

    # Validate the downloaded SRA file
    log_msg "Validating ${SRR}..."                                 # progress marker
    vdb-validate "${RAW_DIR}/${SRR}/${SRR}.sra" 2>&1 | tee "${LOG_DIR}/${SRR}_validate.log"  # check the .sra isn't corrupt
    if [ $? -ne 0 ]; then                                          # if validation reported a problem...
        log_msg "WARNING: Validation failed for ${SRR}, attempting re-download..."  # ...warn...
        rm -rf "${RAW_DIR}/${SRR}"                                 # ...delete the bad copy...
        prefetch "${SRR}" --output-directory "${RAW_DIR}" --max-size 50G --force ALL  # ...and force a fresh download
        check_exit "Re-download failed for ${SRR}"                 # give up if even the retry fails
    fi

    # Convert SRA to paired-end FASTQ
    log_msg "Converting ${SRR} to paired-end FASTQs..."            # progress marker
    fasterq-dump "${RAW_DIR}/${SRR}/${SRR}.sra" \
        --outdir "${RAW_DIR}" \
        --split-3 \
        --threads "${THREADS}" \
        --progress \
        2>&1 | tee "${LOG_DIR}/${SRR}_fasterq.log"                 # --split-3 writes _1/_2 (and singletons) FASTQs
    check_exit "fasterq-dump failed for ${SRR}"                    # stop if conversion failed

    # Compress FASTQs to save space
    log_msg "Compressing FASTQs for ${SRR}..."                     # progress marker
    gzip -f "${RAW_DIR}/${SRR}_1.fastq" 2>/dev/null                # gzip mate 1 (downstream tools read .gz)
    gzip -f "${RAW_DIR}/${SRR}_2.fastq" 2>/dev/null                # gzip mate 2
    # Also compress unpaired reads if they exist
    gzip -f "${RAW_DIR}/${SRR}.fastq" 2>/dev/null                  # gzip singleton reads if any were produced

    # Clean up SRA cache
    rm -rf "${RAW_DIR}/${SRR}"                                     # drop the bulky .sra now that FASTQs exist

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
    fastqc \
        "${R1}" "${R2}" \
        --outdir "${FASTQC_PRE_DIR}" \
        --threads "${THREADS}" \
        --quiet \
        2>&1 | tee "${LOG_DIR}/${SRR}_fastqc_pre.log"              # per-base quality report on the raw reads
    check_exit "FastQC (pre-trimming) failed for ${SRR}"           # stop on failure

    log_msg "Pre-trimming FastQC completed for ${SRR}"             # per-sample done marker
done

##==========================================================================##
##               STEP 4: ADAPTER & QUALITY TRIMMING WITH TRIMMOMATIC         ##
##==========================================================================##

log_msg "========== STEP 4: Trimming reads with Trimmomatic =========="  # announce trimming step

# Find adapter file
ADAPTER_FILE=$(find_adapter_file)                                  # locate the adapter FASTA once, reuse for all samples

for SRR in "${SRR_ACCESSIONS[@]}"; do                              # trim each sample
    R1="${RAW_DIR}/${SRR}_1.fastq.gz"                              # raw mate 1 input
    R2="${RAW_DIR}/${SRR}_2.fastq.gz"                              # raw mate 2 input

    # Output files
    R1_PAIRED="${TRIMMED_DIR}/${SRR}_1_paired.fastq.gz"            # mate 1 reads whose partner also survived
    R1_UNPAIRED="${TRIMMED_DIR}/${SRR}_1_unpaired.fastq.gz"        # mate 1 reads whose partner was dropped
    R2_PAIRED="${TRIMMED_DIR}/${SRR}_2_paired.fastq.gz"            # mate 2 reads whose partner also survived
    R2_UNPAIRED="${TRIMMED_DIR}/${SRR}_2_unpaired.fastq.gz"        # mate 2 reads whose partner was dropped

    if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then                 # if either raw mate is missing...
        log_msg "WARNING: Paired FASTQs not found for ${SRR}, skipping trimming."  # ...warn...
        continue                                                   # ...and skip this sample
    fi

    # Skip if trimmed files already exist
    if [ -f "${R1_PAIRED}" ] && [ -f "${R2_PAIRED}" ]; then        # already trimmed on a previous run?
        log_msg "Trimmed FASTQs for ${SRR} already exist, skipping."  # note the skip
        continue                                                   # move on
    fi

    log_msg "Running Trimmomatic on ${SRR}..."                     # progress marker

    # Build trimmomatic command with or without adapter trimming
    TRIM_STEPS=""                                                  # accumulate the Trimmomatic operation list
    if [ -n "${ADAPTER_FILE}" ] && [ -f "${ADAPTER_FILE}" ]; then  # only add adapter clipping if we found the FASTA
        TRIM_STEPS="ILLUMINACLIP:${ADAPTER_FILE}:2:30:10:2:True "  # clip Illumina adapters (seed/palindrome/simple thresholds)
    fi
    TRIM_STEPS+="LEADING:${TRIM_LEADING} TRAILING:${TRIM_TRAILING} "  # trim low-quality bases off both read ends
    TRIM_STEPS+="SLIDINGWINDOW:${TRIM_SLIDINGWINDOW} "             # sliding-window quality trimming
    TRIM_STEPS+="AVGQUAL:${TRIM_AVGQUAL} "                         # drop reads below the mean-quality cutoff
    TRIM_STEPS+="MINLEN:${TRIM_MINLEN}"                            # drop reads shorter than the minimum length

    trimmomatic PE \
        -threads "${THREADS}" \
        -phred33 \
        -summary "${LOG_DIR}/${SRR}_trimmomatic_summary.txt" \
        "${R1}" "${R2}" \
        "${R1_PAIRED}" "${R1_UNPAIRED}" \
        "${R2_PAIRED}" "${R2_UNPAIRED}" \
        ${TRIM_STEPS} \
        2>&1 | tee "${LOG_DIR}/${SRR}_trimmomatic.log"             # run paired-end trimming with the assembled step list
    check_exit "Trimmomatic failed for ${SRR}"                     # stop on failure

    # Report trimming stats
    if [ -f "${LOG_DIR}/${SRR}_trimmomatic_summary.txt" ]; then    # if a summary file was written...
        log_msg "Trimmomatic summary for ${SRR}:"                  # ...header it in the log...
        cat "${LOG_DIR}/${SRR}_trimmomatic_summary.txt" | while read -r line; do  # ...read it line by line...
            log_msg "  ${line}"                                    # ...and echo each stat into the main log
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
    fastqc \
        "${R1_PAIRED}" "${R2_PAIRED}" \
        --outdir "${FASTQC_POST_DIR}" \
        --threads "${THREADS}" \
        --quiet \
        2>&1 | tee "${LOG_DIR}/${SRR}_fastqc_post.log"             # quality report on the trimmed reads (compare vs pre)
    check_exit "FastQC (post-trimming) failed for ${SRR}"          # stop on failure

    log_msg "Post-trimming FastQC completed for ${SRR}"            # per-sample done marker
done

##==========================================================================##
##               STEP 6: AGGREGATE QC REPORTS WITH MULTIQC                   ##
##==========================================================================##

log_msg "========== STEP 6: Aggregating QC reports with MultiQC =========="  # announce aggregation step

multiqc \
    "${QC_DIR}" "${LOG_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --filename "illumina_qc_report" \
    --title "PRJNA207834 - Illumina QC Summary (HIV-1 Uganda A1/D)" \
    --force \
    2>&1 | tee "${LOG_DIR}/multiqc_illumina.log"                   # roll every FastQC/Trimmomatic report into one HTML
check_exit "MultiQC failed"                                        # stop if aggregation failed

log_msg "MultiQC report generated: ${MULTIQC_DIR}/illumina_qc_report.html"  # point user at the final report

##==========================================================================##
##                          PIPELINE COMPLETE                                ##
##==========================================================================##

log_msg "=========================================="                         # final summary banner
log_msg "  ILLUMINA RAW DATA QC PIPELINE COMPLETE"                            # completion headline
log_msg "=========================================="                         # banner
log_msg "BioProject:         ${BIOPROJECT}"                                   # which BioProject was processed
log_msg "Study:              HIV-1 intersubtype recombinants in Uganda"       # study context
log_msg "Subtypes:           A1, D, and A1-D recombinants"                    # subtypes covered
log_msg "Platform:           Illumina MiSeq, 2x251bp paired-end"             # sequencing platform
log_msg "Samples processed:  ${#SRR_ACCESSIONS[@]}"                           # how many samples went through
log_msg "Raw data:           ${RAW_DIR}"                                      # where raw FASTQs live
log_msg "Trimmed data:       ${TRIMMED_DIR}"                                  # where trimmed FASTQs live
log_msg "QC reports:         ${QC_DIR}"                                       # where QC reports live
log_msg "MultiQC summary:    ${MULTIQC_DIR}/illumina_qc_report.html"          # the aggregated report path
log_msg "=========================================="                         # banner
