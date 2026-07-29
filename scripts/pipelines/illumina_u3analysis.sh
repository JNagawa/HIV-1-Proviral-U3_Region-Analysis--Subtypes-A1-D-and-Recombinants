#!/bin/bash
#SBATCH --job-name=qc_illumina_u3
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --time=12:00:00                 # Illumina short reads process faster
#SBATCH --ntasks=1                      # single multi-threaded task
#SBATCH --cpus-per-task=8               # threads for trimmomatic/fastqc
#SBATCH --mem=64G                       # sufficient for short-read QC

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

# AUTHOR DETAILS
# No.   NAME                       STUDENT NO.          REG NO.
# 1.    Nagawa Jovita              2500726007           2025/HD07/26007

# Activate conda environment
CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"   # usual miniconda hook location
if [ -f "$CONDA_SH" ]; then                          # if that hook exists...
    source "$CONDA_SH"                               # ...source it directly
elif command -v conda >/dev/null 2>&1; then          # else if conda is already on PATH...
    # ...source the hook from its reported base
    source "$(conda info --base)/etc/profile.d/conda.sh"
else
    # neither worked: fail loudly
    echo "ERROR: Conda not found. Please load conda before running this script." >&2
    exit 1                                           # abort rather than run without the tools
fi
conda activate HIV_U3analysis                        # activate the env holding all pipeline tools

##==========================================================================##
##                     CONFIGURATION & VARIABLES                             ##
##==========================================================================##

BIOPROJECT="PRJNA207834"                             # SRA BioProject this pipeline processes
THREADS=${SLURM_CPUS_PER_TASK:-8}    # Use SLURM allocation or default to 8

# SRA accessions for BioProject PRJNA207834
# 24 HIV-1 near-full-genome Illumina MiSeq paired-end samples from Uganda
# Subtypes: A1, D, and A1-D intersubtype recombinants (BSRI)
SRR_ACCESSIONS=(                                     # the 24 samples processed end-to-end
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
TRIM_SLIDINGWINDOW="4:15"  # Sliding window: window_size:quality_threshold
TRIM_MINLEN=50          # Minimum read length after trimming
TRIM_AVGQUAL=15         # Minimum average quality of the read

# Directory structure
BASE_DIR="$(pwd)"                                    # repo root -- MUST be invoked from here
RAW_DIR="${BASE_DIR}/data/raw/illumina"              # downloaded, untrimmed FASTQs
QC_DIR="${BASE_DIR}/results/reports/qc/illumina"     # all Illumina QC output lives under here
FASTQC_PRE_DIR="${QC_DIR}/fastqc_pre"                # FastQC on raw reads
FASTQC_POST_DIR="${QC_DIR}/fastqc_post"              # FastQC on trimmed reads
MULTIQC_DIR="${QC_DIR}/multiqc"                      # aggregated MultiQC report
TRIMMED_DIR="${BASE_DIR}/data/processed/illumina/trimmed_trimmomatic"  # Trimmomatic's trimmed reads
# fastp --dedup-only pass on trimmed reads
DEDUP_DIR="${BASE_DIR}/data/processed/illumina/dedup_trimmomatic"
KRAKEN2_DIR="${BASE_DIR}/data/processed/illumina/kraken2_trimmomatic"  # Kraken2 host-filtered reads
# size-capped Standard DB (human+bacteria+archaea+viral)
KRAKEN2_DB="${BASE_DIR}/data/reference/kraken2_standard_16gb_db"
LOG_DIR="${BASE_DIR}/logs"                           # per-sample tool logs
# BAMs, VCFs, and per-sample consensus FASTAs
ALIGN_DIR="${BASE_DIR}/data/processed/illumina/alignments"
FILTER_DIR="${BASE_DIR}/data/processed/illumina/filtering"  # biological-filtering outputs (Step 8)
MSA_DIR="${BASE_DIR}/data/processed/illumina/msa"    # multiple-sequence-alignment outputs (Step 9)
SUBTYPE_DIR="${BASE_DIR}/data/processed/illumina/subtyping"  # subtyping outputs (Step 10)
# U3 extraction + motif-mapping outputs (Step 11)
MOTIF_DIR="${BASE_DIR}/data/processed/illumina/motifs"
REF_DIR="${BASE_DIR}/data/reference"                 # reference genome + annotation cache
CHECKPOINT_DIR="${BASE_DIR}/checkpoints"             # resume-state markers

##==========================================================================##
##                         HELPER FUNCTIONS                                  ##
##==========================================================================##

log_msg() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ILLUMINA] $1"  # timestamp + tag + the message argument
}

check_exit() {
    # inspect the exit status of the previous command
    if [ $? -ne 0 ]; then
        log_msg "ERROR: $1"                          # log the caller-supplied failure message
        exit 1                                       # abort the pipeline
    fi
}

# Verify a file exists and is non-empty; logs and returns 1 otherwise.
# Used to gate checkpoint creation instead of touching it unconditionally.
verify_nonempty() {
    # arg 1 = file to check, arg 2 = human label for the message
    local f="$1" label="$2"
    if [ ! -s "${f}" ]; then                         # if the file is missing or zero-byte...
        # ...report which step produced nothing...
        log_msg "ERROR: ${label} produced empty/missing output: ${f}"
        return 1                                     # ...and signal failure to the caller
    fi
    return 0                                          # file is non-empty: success
}

# Cheap integrity check for gzip files (catches truncated/corrupted output,
# e.g. from a crash mid-write).
is_valid_gz() {
    # test gzip integrity silently; exit status is the answer
    gzip -t "$1" >/dev/null 2>&1
}

# True only if the checkpoint file exists AND every output file argument is
# non-empty AND (for *.gz files) passes a gzip integrity check. Replaces bare
# "[ -f "$CHKPT" ]"-style resume checks so a checkpoint next to a truncated,
# corrupted, or zero-byte output is never trusted as "done" again.
is_step_done() {
    local chkpt="$1"                                 # arg 1 = the checkpoint marker file
    # remaining args are the output files to validate
    shift
    [ -f "${chkpt}" ] || return 1                    # no checkpoint means not done
    local f                                          # loop variable over output files
    # every declared output must exist and be intact
    for f in "$@"; do
        # missing/empty output invalidates the checkpoint
        [ -s "${f}" ] || return 1
        case "${f}" in                               # for gzip outputs, also verify integrity
            *.gz) is_valid_gz "${f}" || return 1 ;;  # a truncated .gz counts as not-done
        esac
    done
    return 0                                          # checkpoint present and all outputs valid
}

# Pre-flight dependency check; exits immediately with a clear message instead
# of letting a missing tool fail silently deep inside a pipe (this is what let
# a missing minimap2 mark 9 Nanopore samples "done" with 0-byte output).
require_tool() {
    # abort now if the tool isn't on PATH
    command -v "$1" >/dev/null 2>&1 || { log_msg "ERROR: required tool '$1' not found on PATH. Check conda env HIV_U3analysis."; exit 1; }
}

# Locate Trimmomatic adapter file
find_adapter_file() {
    # Search common locations for Trimmomatic adapter files
    # candidate paths for the TruSeq3-PE-2.fa adapter file
    local ADAPTER_LOCATIONS=(
        "${CONDA_PREFIX}/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
        "${CONDA_PREFIX}/share/trimmomatic-*/adapters/TruSeq3-PE-2.fa"
        "/usr/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
        "/usr/local/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
    )

    for pattern in "${ADAPTER_LOCATIONS[@]}"; do      # try each candidate in order
        # Use compgen to expand globs safely
        local found                                  # holds the first matching path, if any
        # expand the glob, take the first hit
        found=$(compgen -G "${pattern}" 2>/dev/null | head -1)
        if [ -n "${found}" ] && [ -f "${found}" ]; then  # if we got a real, existing file...
            echo "${found}"                          # ...print it (the function's return value)...
            return 0                                 # ...and report success
        fi
    done

    # none found: warn but continue
    log_msg "WARNING: TruSeq3-PE-2.fa adapter file not found. Skipping adapter trimming."
    echo ""                                          # print empty so the caller's $() is empty
    return 1                                          # signal not-found to the caller
}

##==========================================================================##
##               STEP 0: PRE-FLIGHT DEPENDENCY CHECK                         ##
##==========================================================================##

log_msg "========== STEP 0: Checking required tools are on PATH =========="

# every executable this pipeline depends on
for TOOL in prefetch vdb-validate fasterq-dump fastqc trimmomatic multiqc \
            bwa samtools bcftools tabix mafft seqkit; do
    require_tool "${TOOL}"                           # abort immediately if any one is missing
done

log_msg "All required tools found."

##==========================================================================##
##               STEP 1: CREATE DIRECTORY STRUCTURE                          ##
##==========================================================================##

log_msg "========== STEP 1: Setting up directory structure =========="

# create every output/log/checkpoint dir up front
mkdir -p "${RAW_DIR}" "${FASTQC_PRE_DIR}" "${FASTQC_POST_DIR}" \
         "${MULTIQC_DIR}" "${TRIMMED_DIR}" "${LOG_DIR}" \
         "${ALIGN_DIR}" "${FILTER_DIR}" "${MSA_DIR}" \
         "${SUBTYPE_DIR}" "${MOTIF_DIR}" "${REF_DIR}" "${CHECKPOINT_DIR}"

log_msg "Directory structure created under: ${BASE_DIR}"

##==========================================================================##
##               STEP 2: DOWNLOAD SRA DATA                                   ##
##==========================================================================##

log_msg "========== STEP 2: Downloading SRA data (${#SRR_ACCESSIONS[@]} samples) =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do                # download each accession in turn
    log_msg "--- Processing ${SRR} ---"

    # Skip if paired-end FASTQs already exist and are intact (not a truncated
    # partial download from an earlier interrupted/crashed run)
    # both mates present and intact?
    if [ -s "${RAW_DIR}/${SRR}_1.fastq.gz" ] && [ -s "${RAW_DIR}/${SRR}_2.fastq.gz" ] && \
       is_valid_gz "${RAW_DIR}/${SRR}_1.fastq.gz" && is_valid_gz "${RAW_DIR}/${SRR}_2.fastq.gz"; then
        log_msg "FASTQs for ${SRR} already exist, skipping download."  # already done
        continue                                     # skip to the next accession
    fi

    # Prefetch SRA file
    log_msg "Prefetching ${SRR}..."
    # download the .sra, mirroring output to a per-sample log
    prefetch "${SRR}" \
        --output-directory "${RAW_DIR}" \
        --max-size 50G \
        --progress \
        2>&1 | tee "${LOG_DIR}/${SRR}_prefetch.log"
    check_exit "prefetch failed for ${SRR}"          # abort if the download failed

    # Validate the downloaded SRA file
    log_msg "Validating ${SRR}..."
    # integrity-check the downloaded .sra
    vdb-validate "${RAW_DIR}/${SRR}/${SRR}.sra" 2>&1 | tee "${LOG_DIR}/${SRR}_validate.log"
    if [ $? -ne 0 ]; then                            # if validation failed...
        log_msg "WARNING: Validation failed for ${SRR}, attempting re-download..."  # ...warn...
        rm -rf "${RAW_DIR}/${SRR}"                   # ...remove the corrupt download...
        # ...and force a clean re-download
        prefetch "${SRR}" --output-directory "${RAW_DIR}" --max-size 50G --force ALL
        check_exit "Re-download failed for ${SRR}"   # abort if even the retry fails
    fi

    # Convert SRA to paired-end FASTQ
    log_msg "Converting ${SRR} to paired-end FASTQs..."
    # extract paired FASTQs (_1/_2 + orphans) from the .sra
    fasterq-dump "${RAW_DIR}/${SRR}/${SRR}.sra" \
        --outdir "${RAW_DIR}" \
        --split-3 \
        --threads "${THREADS}" \
        --progress \
        2>&1 | tee "${LOG_DIR}/${SRR}_fasterq.log"
    check_exit "fasterq-dump failed for ${SRR}"      # abort if conversion failed

    # Compress FASTQs to save space
    log_msg "Compressing FASTQs for ${SRR}..."
    gzip -f "${RAW_DIR}/${SRR}_1.fastq" 2>/dev/null  # compress the forward mate
    gzip -f "${RAW_DIR}/${SRR}_2.fastq" 2>/dev/null  # compress the reverse mate
    # Also compress unpaired reads if they exist
    # compress orphan reads if present (silently skip if not)
    gzip -f "${RAW_DIR}/${SRR}.fastq" 2>/dev/null

    # Clean up SRA cache
    rm -rf "${RAW_DIR}/${SRR}"                        # delete the .sra cache dir to reclaim space

    log_msg "Completed download for ${SRR}"
done

log_msg "All SRA downloads completed."

##==========================================================================##
##               STEP 3: PRE-TRIMMING QC WITH FASTQC                         ##
##==========================================================================##

log_msg "========== STEP 3: Running pre-trimming FastQC =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do                # QC each sample's raw reads
    R1="${RAW_DIR}/${SRR}_1.fastq.gz"                # forward-mate raw FASTQ
    R2="${RAW_DIR}/${SRR}_2.fastq.gz"                # reverse-mate raw FASTQ

    # skip missing/corrupt input
    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ] || ! is_valid_gz "${R1}" || ! is_valid_gz "${R2}"; then
        log_msg "WARNING: Paired FASTQs not found or corrupted for ${SRR}, skipping pre-QC."
        continue
    fi

    log_msg "Running FastQC (pre-trimming) on ${SRR}..."
    # QC both raw mates, logging to a per-sample file
    fastqc \
        "${R1}" "${R2}" \
        --outdir "${FASTQC_PRE_DIR}" \
        --threads "${THREADS}" \
        --quiet \
        2>&1 | tee "${LOG_DIR}/${SRR}_fastqc_pre.log"
    check_exit "FastQC (pre-trimming) failed for ${SRR}"  # abort on FastQC failure

    log_msg "Pre-trimming FastQC completed for ${SRR}"
done

##==========================================================================##
##               STEP 4: ADAPTER & QUALITY TRIMMING WITH TRIMMOMATIC         ##
##==========================================================================##

log_msg "========== STEP 4: Trimming reads with Trimmomatic =========="

# Find adapter file
# resolve the TruSeq3 adapter path once for all samples
ADAPTER_FILE=$(find_adapter_file)

for SRR in "${SRR_ACCESSIONS[@]}"; do                # trim each sample with Trimmomatic
    R1="${RAW_DIR}/${SRR}_1.fastq.gz"                # forward-mate raw FASTQ
    R2="${RAW_DIR}/${SRR}_2.fastq.gz"                # reverse-mate raw FASTQ

    # Output files
    R1_PAIRED="${TRIMMED_DIR}/${SRR}_1_paired.fastq.gz"    # forward reads whose mate also survived
    R1_UNPAIRED="${TRIMMED_DIR}/${SRR}_1_unpaired.fastq.gz"  # forward reads whose mate was dropped
    R2_PAIRED="${TRIMMED_DIR}/${SRR}_2_paired.fastq.gz"    # reverse reads whose mate also survived
    R2_UNPAIRED="${TRIMMED_DIR}/${SRR}_2_unpaired.fastq.gz"  # reverse reads whose mate was dropped

    # skip missing/corrupt input
    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ] || ! is_valid_gz "${R1}" || ! is_valid_gz "${R2}"; then
        log_msg "WARNING: Paired FASTQs not found or corrupted for ${SRR}, skipping trimming."
        continue
    fi

    # Skip if trimmed files already exist and are intact (a truncated gzip
    # from an earlier crashed Trimmomatic run must NOT be trusted as "done")
    # both paired outputs present and intact?
    if [ -s "${R1_PAIRED}" ] && [ -s "${R2_PAIRED}" ] && \
       is_valid_gz "${R1_PAIRED}" && is_valid_gz "${R2_PAIRED}"; then
        log_msg "Trimmed FASTQs for ${SRR} already exist, skipping."
        continue
    fi

    log_msg "Running Trimmomatic on ${SRR}..."

    # Build trimmomatic command with or without adapter trimming
    # accumulate Trimmomatic's ordered trimming steps here
    TRIM_STEPS=""
    # only add adapter clipping if the adapter file exists
    if [ -n "${ADAPTER_FILE}" ] && [ -f "${ADAPTER_FILE}" ]; then
        # adapter-clip step (recommended TruSeq3 PE params)
        TRIM_STEPS="ILLUMINACLIP:${ADAPTER_FILE}:2:30:10:2:True "
    fi
    # append leading/trailing edge trims
    TRIM_STEPS+="LEADING:${TRIM_LEADING} TRAILING:${TRIM_TRAILING} "
    TRIM_STEPS+="SLIDINGWINDOW:${TRIM_SLIDINGWINDOW} "  # append the 4bp/Q15 sliding-window trim
    # append the whole-read Q15 average-quality gate
    TRIM_STEPS+="AVGQUAL:${TRIM_AVGQUAL} "
    # append the 50bp minimum-length filter (must come last)
    TRIM_STEPS+="MINLEN:${TRIM_MINLEN}"

    # run Trimmomatic PE with the built step list, logging output
    trimmomatic PE \
        -Xmx48g \
        -threads "${THREADS}" \
        -phred33 \
        -summary "${LOG_DIR}/${SRR}_trimmomatic_summary.txt" \
        "${R1}" "${R2}" \
        "${R1_PAIRED}" "${R1_UNPAIRED}" \
        "${R2_PAIRED}" "${R2_UNPAIRED}" \
        ${TRIM_STEPS} \
        2>&1 | tee "${LOG_DIR}/${SRR}_trimmomatic.log"
    check_exit "Trimmomatic failed for ${SRR}"       # abort on Trimmomatic failure

    # Report trimming stats
    # if the summary file was written...
    if [ -f "${LOG_DIR}/${SRR}_trimmomatic_summary.txt" ]; then
        log_msg "Trimmomatic summary for ${SRR}:"
        # ...echo each of its lines into the pipeline log
        cat "${LOG_DIR}/${SRR}_trimmomatic_summary.txt" | while read -r line; do
            log_msg "  ${line}"                      # indent under the summary header
        done
    fi

    log_msg "Trimming completed for ${SRR}"
done

##==========================================================================##
##               STEP 5: POST-TRIMMING QC WITH FASTQC                        ##
##==========================================================================##

log_msg "========== STEP 5: Running post-trimming FastQC =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do                # QC each sample's trimmed reads
    R1_PAIRED="${TRIMMED_DIR}/${SRR}_1_paired.fastq.gz"  # trimmed forward paired reads
    R2_PAIRED="${TRIMMED_DIR}/${SRR}_2_paired.fastq.gz"  # trimmed reverse paired reads

    # skip missing/corrupt input
    if [ ! -s "${R1_PAIRED}" ] || [ ! -s "${R2_PAIRED}" ] || ! is_valid_gz "${R1_PAIRED}" || ! is_valid_gz "${R2_PAIRED}"; then
        log_msg "WARNING: Trimmed FASTQs not found or corrupted for ${SRR}, skipping post-QC."
        continue
    fi

    log_msg "Running FastQC (post-trimming) on ${SRR}..."
    # QC the trimmed reads for a before/after comparison
    fastqc \
        "${R1_PAIRED}" "${R2_PAIRED}" \
        --outdir "${FASTQC_POST_DIR}" \
        --threads "${THREADS}" \
        --quiet \
        2>&1 | tee "${LOG_DIR}/${SRR}_fastqc_post.log"
    check_exit "FastQC (post-trimming) failed for ${SRR}"  # abort on FastQC failure

    log_msg "Post-trimming FastQC completed for ${SRR}"
done

##==========================================================================##
##               STEP 6: AGGREGATE QC REPORTS WITH MULTIQC                   ##
##==========================================================================##

log_msg "========== STEP 6: Aggregating QC reports with MultiQC =========="

# aggregate all QC + Trimmomatic logs into one report; --force overwrites a prior run
multiqc \
    "${QC_DIR}" "${LOG_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --filename "illumina_qc_report" \
    --title "PRJNA207834 - Illumina QC Summary (HIV-1 Uganda A1/D)" \
    --force \
    2>&1 | tee "${LOG_DIR}/multiqc_illumina.log"
check_exit "MultiQC failed"                          # abort if aggregation failed

log_msg "MultiQC report generated: ${MULTIQC_DIR}/illumina_qc_report.html"

##==========================================================================##
##  STEP 6.5: DEDUPLICATION + KRAKEN2 HOST/BACTERIAL CONTAMINATION FILTER    ##
##==========================================================================##
# Per the tools review's recommended chain (fastp -> Kraken2 -> SHIVER):
# Trimmomatic has no native dedup capability, so a separate fastp
# dedup-only pass runs first, then Kraken2 filters host/bacterial reads
# against a size-capped Standard database (bacteria+archaea+viral+human+
# UniVec_Core, ~16GB) -- a purely-viral database can't do this job, since
# it has no human/bacterial genomes to match host contamination against.
# scripts/utils/kraken2_filter_reads.sh keeps unclassified reads and reads
# classified as viral, discarding anything descending from Homo sapiens
# (9606) or Bacteria (2) per the database's own bundled nodes.dmp.
#
# NOTE: existing *_bwa.done checkpoints from a prior run (before this step
# existed) point at Trimmomatic's un-deduplicated, unfiltered output --
# they are NOT automatically invalidated here. Clear the relevant
# checkpoints manually to reprocess already-mapped samples through this
# new filtering step.

log_msg "========== STEP 6.5: Deduplication + Kraken2 contamination filtering =========="

# ensure the dedup and Kraken2 output dirs exist
mkdir -p "${DEDUP_DIR}" "${KRAKEN2_DIR}"

if [ ! -s "${KRAKEN2_DB}/nodes.dmp" ]; then          # if the Kraken2 DB isn't present...
    # ...warn and skip; mapping will use un-filtered reads
    log_msg "WARNING: Kraken2 database not found at ${KRAKEN2_DB} -- skipping dedup+Kraken2 filtering. Reference mapping below will fall back to plain Trimmomatic output."
else
    for SRR in "${SRR_ACCESSIONS[@]}"; do            # dedup + host-filter each sample
        R1_PAIRED="${TRIMMED_DIR}/${SRR}_1_paired.fastq.gz"  # Trimmomatic forward paired reads
        R2_PAIRED="${TRIMMED_DIR}/${SRR}_2_paired.fastq.gz"  # Trimmomatic reverse paired reads
        DEDUP_R1="${DEDUP_DIR}/${SRR}_1.dedup.fastq.gz"  # deduped forward mate
        DEDUP_R2="${DEDUP_DIR}/${SRR}_2.dedup.fastq.gz"  # deduped reverse mate
        KRAKEN_R1="${KRAKEN2_DIR}/${SRR}_1.kraken_filtered.fastq.gz"  # host-removed forward mate
        KRAKEN_R2="${KRAKEN2_DIR}/${SRR}_2.kraken_filtered.fastq.gz"  # host-removed reverse mate

        # skip missing/corrupt trimmed input
        if [ ! -s "${R1_PAIRED}" ] || [ ! -s "${R2_PAIRED}" ] || ! is_valid_gz "${R1_PAIRED}" || ! is_valid_gz "${R2_PAIRED}"; then
            log_msg "WARNING: Trimmomatic output not found or corrupted for ${SRR}, skipping dedup+Kraken2."
            continue
        fi

        # checkpoint + valid outputs means already done
        if is_step_done "${CHECKPOINT_DIR}/${SRR}_kraken2.done" "${KRAKEN_R1}" "${KRAKEN_R2}"; then
            log_msg "Dedup+Kraken2 filtering already completed for ${SRR}, skipping."
            continue
        fi

        # skip dedup if a valid deduped pair already exists
        if [ ! -s "${DEDUP_R1}" ] || [ ! -s "${DEDUP_R2}" ] || ! is_valid_gz "${DEDUP_R1}" || ! is_valid_gz "${DEDUP_R2}"; then
            log_msg "Deduplicating ${SRR}..."
            # run the fastp dedup-only wrapper on Trimmomatic's output
            bash "${BASE_DIR}/scripts/utils/fastp_dedup.sh" \
                "${R1_PAIRED}" "${R2_PAIRED}" "${DEDUP_R1}" "${DEDUP_R2}" \
                "${DEDUP_DIR}/${SRR}_dedup_fastp" \
                2>&1 | tee "${LOG_DIR}/${SRR}_dedup.log"
            check_exit "fastp dedup failed for ${SRR}"  # abort on dedup failure
        fi

        log_msg "Running Kraken2 on ${SRR}..."
        # drop host/bacterial reads from the deduped pair
        bash "${BASE_DIR}/scripts/utils/kraken2_filter_reads.sh" \
            "${DEDUP_R1}" "${DEDUP_R2}" "${KRAKEN2_DIR}" "${SRR}" "${KRAKEN2_DB}" \
            2>&1 | tee "${LOG_DIR}/${SRR}_kraken2.log"
        check_exit "Kraken2 filtering failed for ${SRR}"  # abort on Kraken2 failure
        # sanity-check the filter produced output
        verify_nonempty "${KRAKEN_R1}" "Kraken2-filtered reads for ${SRR}"
        # mark this sample's filtering complete for resumes
        touch "${CHECKPOINT_DIR}/${SRR}_kraken2.done"

        log_msg "Dedup+Kraken2 filtering completed for ${SRR}"
    done
fi

##==========================================================================##
##               STEP 7: REFERENCE MAPPING (BWA & SAMTOOLS)                  ##
##==========================================================================##

log_msg "========== STEP 7: Reference Mapping with BWA-MEM =========="

# HXB2 reference accession (the HIV-1 coordinate standard)
REF_ACC="K03455.1"
REF_FASTA="${REF_DIR}/${REF_ACC}.fasta"              # local path to the reference FASTA

# ensure alignment/reference/checkpoint dirs exist
mkdir -p "${ALIGN_DIR}/bam" "${REF_DIR}" "${CHECKPOINT_DIR}"

# download + index the reference only if it's not already present
if [ ! -s "${REF_FASTA}" ]; then
    log_msg "Downloading HXB2 reference (${REF_ACC})..."
    # fetch the HXB2 FASTA from NCBI
    wget -q -O "${REF_FASTA}" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${REF_ACC}&rettype=fasta&retmode=text"
    verify_nonempty "${REF_FASTA}" "HXB2 reference download"  # ensure the download produced content
    check_exit "HXB2 reference download failed"      # abort if the download failed
    # build the BWA index needed for mapping
    bwa index "${REF_FASTA}" > "${LOG_DIR}/bwa_index.log" 2>&1
    check_exit "bwa index failed on HXB2 reference"  # abort on index failure
    samtools faidx "${REF_FASTA}"                    # build the .fai index needed by bcftools
    check_exit "samtools faidx failed on HXB2 reference"  # abort on faidx failure
fi

MAPPED_COUNT=0                                       # running count of successfully mapped samples
# running count of samples that had input to map
TOTAL_TO_MAP=0

for SRR in "${SRR_ACCESSIONS[@]}"; do                # map each sample to HXB2 and call a consensus
    # Prefer deduplicated, Kraken2-filtered reads (STEP 6.5); fall back to
    # plain Trimmomatic output if filtering was skipped (e.g. Kraken2 DB
    # not present yet).
    # preferred forward input: host-filtered reads
    R1_PAIRED="${KRAKEN2_DIR}/${SRR}_1.kraken_filtered.fastq.gz"
    # preferred reverse input: host-filtered reads
    R2_PAIRED="${KRAKEN2_DIR}/${SRR}_2.kraken_filtered.fastq.gz"
    # if filtering was skipped/absent...
    if [ ! -s "${R1_PAIRED}" ] || [ ! -s "${R2_PAIRED}" ]; then
        # ...fall back to plain trimmed forward reads
        R1_PAIRED="${TRIMMED_DIR}/${SRR}_1_paired.fastq.gz"
        # ...fall back to plain trimmed reverse reads
        R2_PAIRED="${TRIMMED_DIR}/${SRR}_2_paired.fastq.gz"
    fi
    BAM_OUT="${ALIGN_DIR}/bam/${SRR}.sorted.bam"     # sorted alignment output
    VCF_OUT="${ALIGN_DIR}/bam/${SRR}.vcf.gz"         # variant calls used to build the consensus
    CONSENSUS_OUT="${ALIGN_DIR}/${SRR}_consensus.fasta"  # per-sample consensus sequence
    CHKPT="${CHECKPOINT_DIR}/${SRR}_bwa.done"        # resume marker for this sample's mapping

    # no valid input -> nothing to map
    if [ ! -s "${R1_PAIRED}" ] || [ ! -s "${R2_PAIRED}" ] || ! is_valid_gz "${R1_PAIRED}" || ! is_valid_gz "${R2_PAIRED}"; then
        continue
    fi

    # this sample has input, so count it toward the total
    TOTAL_TO_MAP=$((TOTAL_TO_MAP + 1))

    # checkpoint + valid BAM/consensus means already mapped
    if is_step_done "${CHKPT}" "${BAM_OUT}" "${CONSENSUS_OUT}"; then
        log_msg "Mapping already completed for ${SRR}, skipping."
        MAPPED_COUNT=$((MAPPED_COUNT + 1))           # count it as already-done
        continue
    fi

    log_msg "Aligning ${SRR} to HXB2..."
    # align, convert to BAM, and coordinate-sort in one pipe
    bwa mem -t "${THREADS}" "${REF_FASTA}" "${R1_PAIRED}" "${R2_PAIRED}" 2> "${LOG_DIR}/${SRR}_bwa.log" | \
        samtools view -@ "${THREADS}" -b - | \
        samtools sort -@ "${THREADS}" -o "${BAM_OUT}"
    # mapping failed or produced no BAM
    if [ $? -ne 0 ] || ! verify_nonempty "${BAM_OUT}" "alignment for ${SRR}"; then
        log_msg "WARNING: mapping failed for ${SRR}, skipping (not marking done). See ${LOG_DIR}/${SRR}_bwa.log"
        # skip without a checkpoint so it retries next run
        continue
    fi

    # build the BAM index needed by downstream tools
    samtools index "${BAM_OUT}"
    if [ $? -ne 0 ]; then                            # if indexing failed...
        log_msg "WARNING: samtools index failed for ${SRR}, skipping (not marking done)."
        continue                                     # ...skip without marking done
    fi

    # Generate consensus sequence for downstream steps
    # pile up reads and call variants against HXB2
    bcftools mpileup -Ou -f "${REF_FASTA}" "${BAM_OUT}" 2> "${LOG_DIR}/${SRR}_bcftools.log" | \
        bcftools call -c -Oz -o "${VCF_OUT}" 2>> "${LOG_DIR}/${SRR}_bcftools.log"
    # variant calling failed or produced no VCF
    if [ $? -ne 0 ] || ! verify_nonempty "${VCF_OUT}" "VCF for ${SRR}"; then
        log_msg "WARNING: variant calling failed for ${SRR}, skipping (not marking done). See ${LOG_DIR}/${SRR}_bcftools.log"
        continue
    fi

    # index the VCF so bcftools consensus can use it
    tabix -p vcf "${VCF_OUT}" 2>> "${LOG_DIR}/${SRR}_bcftools.log"
    if [ $? -ne 0 ]; then                            # if indexing the VCF failed...
        log_msg "WARNING: tabix indexing failed for ${SRR}, skipping (not marking done)."
        continue
    fi

    # apply the sample's variants onto HXB2 to build its consensus
    cat "${REF_FASTA}" | bcftools consensus "${VCF_OUT}" > "${CONSENSUS_OUT}" 2>> "${LOG_DIR}/${SRR}_bcftools.log"
    # consensus step failed or empty
    if [ $? -ne 0 ] || ! verify_nonempty "${CONSENSUS_OUT}" "consensus for ${SRR}"; then
        log_msg "WARNING: consensus generation failed for ${SRR}, skipping (not marking done). See ${LOG_DIR}/${SRR}_bcftools.log"
        continue
    fi

    # bcftools consensus does not rename the sequence header -- it keeps the
    # reference's own header verbatim, so every sample's consensus would
    # otherwise carry the identical ">${REF_ACC} ..." header. Left unfixed,
    # concatenating all samples for MSA (Step 9) would make every sequence
    # indistinguishable by ID. Rename the header to this sample's accession.
    # rewrite the FASTA header to this sample's accession
    sed -i "1s/.*/>${SRR}/" "${CONSENSUS_OUT}"

    touch "${CHKPT}"                                 # mark mapping+consensus complete for resumes
    MAPPED_COUNT=$((MAPPED_COUNT + 1))               # count this sample as successfully mapped
    log_msg "Mapping and consensus generation completed for ${SRR}"
done

log_msg "Step 7 complete: ${MAPPED_COUNT}/${TOTAL_TO_MAP} samples mapped successfully."

##==========================================================================##
##               STEP 8: BIOLOGICAL FILTERING                                ##
##==========================================================================##

log_msg "========== STEP 8: Biological Filtering (Poplars & HIVSeqinR) =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do                # biological-filter each sample's consensus
    CONSENSUS="${ALIGN_DIR}/${SRR}_consensus.fasta"  # this sample's consensus from Step 7
    CHKPT="${CHECKPOINT_DIR}/${SRR}_filtering.done"  # resume marker for this step

    if [ ! -f "${CONSENSUS}" ]; then continue; fi    # nothing to filter if no consensus exists
    if [ -f "${CHKPT}" ]; then continue; fi          # already filtered -> skip

    log_msg "Running Hypermut 3 (Poplars) & HIVSeqinR on ${SRR}..."
    # Note: These are generalized commands. Adjust paths/params for your cluster environment.
    # poplars hypermut --input "${CONSENSUS}" --output "${FILTER_DIR}/${SRR}_poplars.fasta" || true
    # Rscript /path/to/HIVSeqinR/HIVSeqinR.R --input "${FILTER_DIR}/${SRR}_poplars.fasta" --outdir "${FILTER_DIR}/${SRR}_seqinr" || true

    touch "${CHKPT}"                                 # mark the (currently placeholder) step done
done

##==========================================================================##
##               STEP 9: MULTIPLE SEQUENCE ALIGNMENT                         ##
##==========================================================================##

log_msg "========== STEP 9: Multiple Sequence Alignment (MAFFT) =========="

# all consensuses (plus HXB2) concatenated for alignment
COMBINED_FASTA="${MSA_DIR}/all_filtered_consensus.fasta"
MSA_OUT="${MSA_DIR}/aligned_consensus.fasta"         # the resulting multiple sequence alignment
CHKPT="${CHECKPOINT_DIR}/mafft_alignment.done"       # resume marker for the alignment step

# only align if not already done with valid output
if ! is_step_done "${CHKPT}" "${MSA_OUT}"; then
    log_msg "Combining consensus sequences (with HXB2 as coordinate anchor) and aligning with MAFFT L-INS-i..."
    # HXB2 is prepended so it shares the alignment's column space -- Step 11's
    # U3 extraction anchors on HXB2's own row to convert its real annotated
    # genome coordinates into alignment-column coordinates.
    # concatenate HXB2 + all consensuses (tolerate no matches)
    cat "${REF_FASTA}" "${ALIGN_DIR}"/*_consensus.fasta > "${COMBINED_FASTA}" 2>/dev/null || true
    if [ -s "${COMBINED_FASTA}" ]; then              # only align if there's something to align
        # L-INS-i high-accuracy alignment
        mafft --localpair --maxiterate 1000 --thread "${THREADS}" "${COMBINED_FASTA}" > "${MSA_OUT}" 2> "${LOG_DIR}/mafft.log"
        # only checkpoint if alignment produced output
        if verify_nonempty "${MSA_OUT}" "MAFFT alignment"; then
            touch "${CHKPT}"                         # mark the alignment step complete
        else
            log_msg "WARNING: MAFFT alignment failed, not marking step done. See ${LOG_DIR}/mafft.log"
        fi
    else
        log_msg "No consensus sequences found for MSA."  # nothing was produced upstream
    fi
fi

##==========================================================================##
##               STEP 10: SUBTYPE & RECOMBINATION DETECTION                  ##
##==========================================================================##

log_msg "========== STEP 10: Subtyping (jpHMM, IQ-TREE 2, COMET, REGA) =========="

CHKPT="${CHECKPOINT_DIR}/subtyping.done"             # resume marker for the subtyping step

if [ ! -f "${CHKPT}" ] && [ -f "${MSA_OUT}" ]; then  # only run once, and only if the MSA exists
    log_msg "Running subtyping pipeline..."
    # Generalized commands
    # jpHMM -v HIV -s "${MSA_OUT}" -o "${SUBTYPE_DIR}/jphmm_out" || true
    # iqtree2 -s "${MSA_OUT}" -m MFP -B 1000 -T "${THREADS}" --prefix "${SUBTYPE_DIR}/iqtree" || true
    # comet -i "${MSA_OUT}" -o "${SUBTYPE_DIR}/comet_results.csv" || true
    touch "${CHKPT}"                                 # mark the (currently placeholder) step done
fi

##==========================================================================##
##               STEP 11: U3 EXTRACTION & MOTIF MAPPING                      ##
##==========================================================================##

log_msg "========== STEP 11: U3 Extraction & Motif Mapping =========="

U3_CHKPT="${CHECKPOINT_DIR}/u3_extraction.done"      # resume marker for U3 extraction
# extracted U3 columns still in alignment (gapped) form
U3_GAPPED="${MOTIF_DIR}/U3_aligned.fasta"
U3_EXTRACTED="${MOTIF_DIR}/U3_extracted.fasta"       # final ungapped per-sample U3 sequences
U3_WARNINGS="${MOTIF_DIR}/u3_extraction_warnings.log"  # outlier-check warnings from the extractor

# only extract if not already done with valid output
if ! is_step_done "${U3_CHKPT}" "${U3_EXTRACTED}"; then
    if [ -s "${MSA_OUT}" ]; then                     # can only extract U3 if the MSA exists
        log_msg "Extracting 5' LTR U3 region, anchored to HXB2's real annotated coordinates..."
        # extract U3 by mapping HXB2's annotated coords to alignment columns
        bash "${BASE_DIR}/scripts/utils/extract_u3_by_hxb2_anchor.sh" \
            --alignment "${MSA_OUT}" \
            --hxb2-id "${REF_ACC}" \
            --gb-cache "${REF_DIR}/${REF_ACC}.gb" \
            --out-gapped "${U3_GAPPED}" \
            --out "${U3_EXTRACTED}" \
            --warnings-log "${U3_WARNINGS}" \
            > "${LOG_DIR}/u3_extraction.log" 2>&1
        # only checkpoint if extraction produced output
        if verify_nonempty "${U3_EXTRACTED}" "U3 extraction"; then
            touch "${U3_CHKPT}"                      # mark U3 extraction complete
            if [ -s "${U3_WARNINGS}" ]; then         # surface any outlier warnings to the log
                log_msg "WARNING: some samples flagged by the U3 extraction outlier check, see ${U3_WARNINGS}"
            fi
        else
            log_msg "WARNING: U3 extraction failed, not marking step done. See ${LOG_DIR}/u3_extraction.log"
        fi
    else
        log_msg "No MSA output found, skipping U3 extraction."  # nothing upstream to extract from
    fi
fi

CHKPT="${CHECKPOINT_DIR}/motif_mapping.done"         # resume marker for the motif-mapping step

# only run once, and only if U3 sequences exist
if [ ! -f "${CHKPT}" ] && [ -s "${U3_EXTRACTED}" ]; then
    log_msg "Running motif scanning (FIMO) and G-quadruplex prediction (gquad + pqsfinder) on the extracted U3 regions..."

    # 1. TFBS motif scanning with FIMO against the 6 core JASPAR TFs --
    #    FIMO is the comparison harness's recommended winner (fastest,
    #    cleanest output; see results/motif_mapping/illumina/ease_of_use_notes.md).
    FIMO_OUT="${MOTIF_DIR}/fimo_out"                 # FIMO output dir
    # scan U3 for the 6 core TF binding motifs (p<1e-4)
    fimo --oc "${FIMO_OUT}" --thresh 1e-4 "${REF_DIR}/jaspar/core6_pfms.meme" "${U3_EXTRACTED}" \
        > "${LOG_DIR}/motif_fimo.log" 2>&1
    if [ -s "${FIMO_OUT}/fimo.tsv" ]; then           # FIMO writes hits to fimo.tsv
        log_msg "FIMO motif scan complete: ${FIMO_OUT}/fimo.tsv"
    else
        log_msg "WARNING: FIMO produced no output, see ${LOG_DIR}/motif_fimo.log"
    fi

    # 2. G-quadruplex prediction: gquad (primary, per the tools review) and
    #    pqsfinder (confirmation pass, imperfection-tolerant scoring).
    GQUAD_OUT="${MOTIF_DIR}/gquad_out.gff3"          # gquad G4-prediction output
    bash "${BASE_DIR}/scripts/motif_mapping/illumina/run_gquad.sh" "${U3_EXTRACTED}" "${GQUAD_OUT}" \
        > "${LOG_DIR}/motif_gquad.log" 2>&1          # predict G-quadruplexes in U3 (primary tool)
    if [ -s "${GQUAD_OUT}" ]; then                   # non-empty output means it ran
        log_msg "gquad G-quadruplex prediction complete: ${GQUAD_OUT}"
    else
        log_msg "WARNING: gquad produced no output, see ${LOG_DIR}/motif_gquad.log"
    fi

    PQSFINDER_OUT="${MOTIF_DIR}/pqsfinder_out.gff3"  # pqsfinder G4-prediction output
    # confirmation G4 pass with imperfection-tolerant scoring
    bash "${BASE_DIR}/scripts/motif_mapping/illumina/run_pqsfinder.sh" "${U3_EXTRACTED}" "${PQSFINDER_OUT}" \
        > "${LOG_DIR}/motif_pqsfinder.log" 2>&1
    if [ -s "${PQSFINDER_OUT}" ]; then               # non-empty output means it ran
        log_msg "pqsfinder G-quadruplex prediction complete: ${PQSFINDER_OUT}"
    else
        log_msg "WARNING: pqsfinder produced no output, see ${LOG_DIR}/motif_pqsfinder.log"
    fi

    touch "${CHKPT}"                                 # mark motif mapping complete
fi

##==========================================================================##
##                          PIPELINE COMPLETE                                ##
##==========================================================================##

log_msg "=========================================="
log_msg "  ILLUMINA END-TO-END PIPELINE COMPLETE   "
log_msg "=========================================="
log_msg "BioProject:         ${BIOPROJECT}"
log_msg "Samples processed:  ${#SRR_ACCESSIONS[@]}"
log_msg "Raw data:           ${RAW_DIR}"
log_msg "Alignments:         ${ALIGN_DIR}"
log_msg "MSA:                ${MSA_DIR}"
log_msg "Motifs:             ${MOTIF_DIR}"
log_msg "=========================================="
