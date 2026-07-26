#!/bin/bash
#SBATCH --job-name=qc_u3analysis        # job name
#SBATCH --output=logs/slurm-%j.out      # output file
#SBATCH --error=logs/slurm-%j.err       # error log
#SBATCH --time=24:00:00                 # expected runtime
#SBATCH --ntasks=1                      # single multi-threaded task
#SBATCH --cpus-per-task=8               # number of threads for bwa/samtools/bcftools
#SBATCH --mem=32G                       # adjust based on genome size. (memory per thread x threads) + buffer. -> (1-2GB /thread) /
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

# Author: Jovita Nagawa - MSc Bioinformatics student

# Activate conda environment
CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"   # usual miniconda hook location
if [ -f "$CONDA_SH" ]; then                          # if that hook exists...
    source "$CONDA_SH"                               # ...source it directly
elif command -v conda >/dev/null 2>&1; then          # else if conda is already on PATH...
    source "$(conda info --base)/etc/profile.d/conda.sh"  # ...source the hook from its reported base
else
    echo "ERROR: Conda not found. Please load conda before running this script." >&2  # neither worked: fail loudly
    exit 1                                           # abort rather than run without the tools
fi
conda activate HIV_U3analysis                        # activate the env holding all pipeline tools

##==========================================================================##
##                     CONFIGURATION & VARIABLES                             ##
##==========================================================================##

BIOPROJECT="PRJNA765218"                             # SRA BioProject this pipeline processes
THREADS=${SLURM_CPUS_PER_TASK:-8}    # Use SLURM allocation or default to 8

# SRA accessions for BioProject PRJNA765218 (NanoHIV - Oxford Nanopore GridION)
# 9 HIV-1 proviral genome samples from Stellenbosch University
SRR_ACCESSIONS=(                                     # the 9 Nanopore samples processed end-to-end
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
MIN_QUALITY=7         # Minimum average read quality (Phred score)
MIN_LENGTH=200        # Minimum read length in bp

# Directory structure
BASE_DIR="$(pwd)"                                    # repo root -- MUST be invoked from here
RAW_DIR="${BASE_DIR}/data/raw/oxnano"                # downloaded, unfiltered FASTQs
QC_DIR="${BASE_DIR}/results/reports/qc/oxnano"       # all Nanopore QC output lives under here
NANOPLOT_PRE_DIR="${QC_DIR}/nanoplot_pre"            # NanoPlot on raw reads
NANOPLOT_POST_DIR="${QC_DIR}/nanoplot_post"          # NanoPlot on filtered reads
NANOQC_DIR="${QC_DIR}/nanoqc"                        # nanoQC per-base quality output
NANOSTAT_DIR="${QC_DIR}/nanostat"                    # NanoStat text summaries (pre + post)
MULTIQC_DIR="${QC_DIR}/multiqc"                      # aggregated MultiQC report
FILTERED_DIR="${BASE_DIR}/data/processed/oxnano/filtered"  # adapter-trimmed + quality/length-filtered reads
LOG_DIR="${BASE_DIR}/logs"                           # per-sample tool logs
ALIGN_DIR="${BASE_DIR}/data/processed/oxnano/alignments"  # BAMs, VCFs, and per-sample consensus FASTAs
FILTER_DIR="${BASE_DIR}/data/processed/oxnano/filtering"  # biological-filtering outputs (Step 8)
MSA_DIR="${BASE_DIR}/data/processed/oxnano/msa"      # multiple-sequence-alignment outputs (Step 9)
SUBTYPE_DIR="${BASE_DIR}/data/processed/oxnano/subtyping"  # subtyping outputs (Step 10)
MOTIF_DIR="${BASE_DIR}/data/processed/oxnano/motifs"  # U3 extraction + motif-mapping outputs (Step 11)
REF_DIR="${BASE_DIR}/data/reference"                 # reference genome + annotation cache
CHECKPOINT_DIR="${BASE_DIR}/checkpoints"             # resume-state markers

##==========================================================================##
##                         HELPER FUNCTIONS                                  ##
##==========================================================================##

log_msg() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"         # timestamp + the message argument
}

check_exit() {
    if [ $? -ne 0 ]; then                            # inspect the exit status of the previous command
        log_msg "ERROR: $1"                          # log the caller-supplied failure message
        exit 1                                       # abort the pipeline
    fi
}

# Verify a file exists and is non-empty; logs and returns 1 otherwise.
# Used to gate checkpoint creation instead of touching it unconditionally.
verify_nonempty() {
    local f="$1" label="$2"                          # arg 1 = file to check, arg 2 = human label for the message
    if [ ! -s "${f}" ]; then                         # if the file is missing or zero-byte...
        log_msg "ERROR: ${label} produced empty/missing output: ${f}"  # ...report which step produced nothing...
        return 1                                     # ...and signal failure to the caller
    fi
    return 0                                          # file is non-empty: success
}

# Cheap integrity check for gzip files (catches truncated/corrupted output,
# e.g. from a crash mid-write).
is_valid_gz() {
    gzip -t "$1" >/dev/null 2>&1                      # test gzip integrity silently; exit status is the answer
}

# True only if the checkpoint file exists AND every output file argument is
# non-empty AND (for *.gz files) passes a gzip integrity check. Replaces bare
# "[ -f "$CHKPT" ]"-style resume checks so a checkpoint next to a truncated,
# corrupted, or zero-byte output is never trusted as "done" again.
is_step_done() {
    local chkpt="$1"                                 # arg 1 = the checkpoint marker file
    shift                                            # remaining args are the output files to validate
    [ -f "${chkpt}" ] || return 1                    # no checkpoint means not done
    local f                                          # loop variable over output files
    for f in "$@"; do                                # every declared output must exist and be intact
        [ -s "${f}" ] || return 1                    # missing/empty output invalidates the checkpoint
        case "${f}" in                               # for gzip outputs, also verify integrity
            *.gz) is_valid_gz "${f}" || return 1 ;;  # a truncated .gz counts as not-done
        esac
    done
    return 0                                          # checkpoint present and all outputs valid
}

# Pre-flight dependency check; exits immediately with a clear message instead
# of letting a missing tool fail silently deep inside a pipe (this is exactly
# what let a missing minimap2 mark all 9 samples "done" with 0-byte output).
require_tool() {
    command -v "$1" >/dev/null 2>&1 || { log_msg "ERROR: required tool '$1' not found on PATH. Check conda env HIV_U3analysis."; exit 1; }  # abort now if the tool isn't on PATH
}

##==========================================================================##
##               STEP 0: PRE-FLIGHT DEPENDENCY CHECK                         ##
##==========================================================================##

log_msg "========== STEP 0: Checking required tools are on PATH =========="

for TOOL in prefetch vdb-validate fasterq-dump NanoPlot nanoQC NanoStat \
            NanoFilt multiqc minimap2 samtools bcftools tabix mafft seqkit; do  # every executable this pipeline depends on
    require_tool "${TOOL}"                           # abort immediately if any one is missing
done

# porechop_abi/porechop are checked at call time, not here: Step 4 already
# has a legitimate runtime fallback from porechop_abi to plain porechop, so a
# hard pre-flight requirement here would defeat that fallback.
if ! command -v porechop_abi >/dev/null 2>&1 && ! command -v porechop >/dev/null 2>&1; then  # need at least one adapter trimmer available
    log_msg "ERROR: neither porechop_abi nor porechop found on PATH."
    exit 1                                           # abort if neither variant exists
fi

log_msg "All required tools found."

##==========================================================================##
##               STEP 1: CREATE DIRECTORY STRUCTURE                          ##
##==========================================================================##

log_msg "========== STEP 1: Setting up directory structure =========="

mkdir -p "${RAW_DIR}" "${NANOPLOT_PRE_DIR}" "${NANOPLOT_POST_DIR}" \
         "${NANOQC_DIR}" "${NANOSTAT_DIR}" "${MULTIQC_DIR}" \
         "${FILTERED_DIR}" "${LOG_DIR}" \
         "${ALIGN_DIR}" "${FILTER_DIR}" "${MSA_DIR}" \
         "${SUBTYPE_DIR}" "${MOTIF_DIR}" "${REF_DIR}" "${CHECKPOINT_DIR}"  # create every output/log/checkpoint dir up front

log_msg "Directory structure created under: ${BASE_DIR}"

##==========================================================================##
##               STEP 2: DOWNLOAD SRA DATA                                   ##
##==========================================================================##

log_msg "========== STEP 2: Downloading SRA data (${#SRR_ACCESSIONS[@]} samples) =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do                # download each accession in turn
    log_msg "--- Processing ${SRR} ---"

    # Skip if FASTQ already exists and is intact (not a truncated partial
    # download/compression from an earlier interrupted/crashed run)
    if [ -s "${RAW_DIR}/${SRR}.fastq" ] || \
       { [ -s "${RAW_DIR}/${SRR}.fastq.gz" ] && is_valid_gz "${RAW_DIR}/${SRR}.fastq.gz"; }; then  # plain FASTQ present, or a valid gzipped one?
        log_msg "FASTQ for ${SRR} already exists, skipping download."  # already done
        continue                                     # skip to the next accession
    fi

    # Prefetch SRA file (with retry)
    log_msg "Prefetching ${SRR}..."
    prefetch "${SRR}" \
        --output-directory "${RAW_DIR}" \
        --max-size 50G \
        --progress \
        2>&1 | tee "${LOG_DIR}/${SRR}_prefetch.log"  # download the .sra, mirroring output to a per-sample log
    check_exit "prefetch failed for ${SRR}"          # abort if the download failed

    # Validate the downloaded SRA file
    log_msg "Validating ${SRR}..."
    vdb-validate "${RAW_DIR}/${SRR}/${SRR}.sra" 2>&1 | tee "${LOG_DIR}/${SRR}_validate.log"  # integrity-check the downloaded .sra
    if [ $? -ne 0 ]; then                            # if validation failed...
        log_msg "WARNING: Validation failed for ${SRR}, attempting re-download..."  # ...warn...
        rm -rf "${RAW_DIR}/${SRR}"                   # ...remove the corrupt download...
        prefetch "${SRR}" --output-directory "${RAW_DIR}" --max-size 50G --force ALL  # ...and force a clean re-download
        check_exit "Re-download failed for ${SRR}"   # abort if even the retry fails
    fi

    # Convert SRA to FASTQ (single-end for Nanopore)
    log_msg "Converting ${SRR} to FASTQ..."
    fasterq-dump "${RAW_DIR}/${SRR}/${SRR}.sra" \
        --outdir "${RAW_DIR}" \
        --threads "${THREADS}" \
        --progress \
        2>&1 | tee "${LOG_DIR}/${SRR}_fasterq.log"   # extract a single-end FASTQ from the .sra (Nanopore is single-end)
    check_exit "fasterq-dump failed for ${SRR}"      # abort if conversion failed

    # Compress FASTQ to save space
    log_msg "Compressing ${SRR}.fastq..."
    gzip -f "${RAW_DIR}/${SRR}.fastq"                # compress the FASTQ in place (-f overwrites any stale .gz)

    # Clean up SRA cache to save disk space
    rm -rf "${RAW_DIR}/${SRR}"                        # delete the .sra cache dir to reclaim space

    log_msg "Completed download for ${SRR}"
done

log_msg "All SRA downloads completed."

##==========================================================================##
##               STEP 3: PRE-FILTERING QC (Raw Data Assessment)              ##
##==========================================================================##

log_msg "========== STEP 3: Running pre-filtering QC =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do                # QC each sample's raw reads
    FASTQ="${RAW_DIR}/${SRR}.fastq.gz"               # this sample's raw single-end FASTQ

    if [ ! -s "${FASTQ}" ] || ! is_valid_gz "${FASTQ}"; then  # skip missing/corrupt input
        log_msg "WARNING: ${FASTQ} not found or corrupted, skipping QC for ${SRR}."
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
        2>&1 | tee "${LOG_DIR}/${SRR}_nanoplot_pre.log"  # length/quality plots of the raw reads (baseline)
    check_exit "NanoPlot failed for ${SRR}"          # abort on NanoPlot failure

    # --- NanoQC: Per-base quality across read positions ---
    log_msg "Running NanoQC on ${SRR}..."
    nanoQC \
        -o "${NANOQC_DIR}/${SRR}" \
        "${FASTQ}" \
        2>&1 | tee "${LOG_DIR}/${SRR}_nanoqc.log"    # per-base quality across read positions on raw reads
    check_exit "NanoQC failed for ${SRR}"            # abort on NanoQC failure

    # --- NanoStat: Quick text summary statistics ---
    log_msg "Running NanoStat on ${SRR}..."
    NanoStat \
        --fastq "${FASTQ}" \
        --outdir "${NANOSTAT_DIR}" \
        --name "${SRR}_pre_stats.txt" \
        --threads "${THREADS}" \
        2>&1 | tee "${LOG_DIR}/${SRR}_nanostat_pre.log"  # text summary stats of the raw reads
    check_exit "NanoStat failed for ${SRR}"          # abort on NanoStat failure

    log_msg "Pre-filtering QC completed for ${SRR}"
done

##==========================================================================##
##               STEP 4: READ FILTERING & TRIMMING                           ##
##==========================================================================##

log_msg "========== STEP 4: Filtering and trimming reads =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do                # trim + filter each sample
    FASTQ="${RAW_DIR}/${SRR}.fastq.gz"               # raw input reads
    TRIMMED="${FILTERED_DIR}/${SRR}_trimmed.fastq.gz"  # intermediate adapter-trimmed reads
    FILTERED="${FILTERED_DIR}/${SRR}_filtered.fastq.gz"  # final quality/length-filtered reads

    if [ ! -s "${FASTQ}" ] || ! is_valid_gz "${FASTQ}"; then  # skip missing/corrupt input
        log_msg "WARNING: ${FASTQ} not found or corrupted, skipping filtering for ${SRR}."
        continue
    fi

    # Skip if filtered file already exists and is intact
    if [ -s "${FILTERED}" ] && is_valid_gz "${FILTERED}"; then  # valid filtered output already present?
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
        2>&1 | tee "${LOG_DIR}/${SRR}_porechop.log"  # adapter-trim + de-chimera the raw reads

    # If porechop_abi is not available, try porechop
    if [ $? -ne 0 ]; then                            # if porechop_abi failed/absent...
        log_msg "Porechop_ABI not available, trying porechop..."
        porechop \
            --input "${FASTQ}" \
            --output "${TRIMMED}" \
            --threads "${THREADS}" \
            2>&1 | tee "${LOG_DIR}/${SRR}_porechop.log"  # ...fall back to plain porechop

        # If porechop also fails, use raw file for filtering
        if [ $? -ne 0 ]; then                        # if both trimmers failed...
            log_msg "WARNING: Adapter trimming unavailable. Using raw reads for filtering."
            TRIMMED="${FASTQ}"                       # ...feed the raw reads straight into NanoFilt
        fi
    fi

    # --- NanoFilt: Quality and length filtering ---
    log_msg "Running NanoFilt on ${SRR} (Q>=${MIN_QUALITY}, len>=${MIN_LENGTH})..."
    gunzip -c "${TRIMMED}" | \
        NanoFilt \
            --quality "${MIN_QUALITY}" \
            --length "${MIN_LENGTH}" | \
        gzip > "${FILTERED}"                          # decompress -> quality/length-filter -> recompress (NanoFilt is stdin/stdout)
    check_exit "NanoFilt failed for ${SRR}"          # abort on NanoFilt failure

    # Clean up intermediate trimmed file (if different from raw)
    if [ "${TRIMMED}" != "${FASTQ}" ]; then          # only delete a real trimmed file, never the raw input
        rm -f "${TRIMMED}"                           # drop the intermediate trimmed file to save space
    fi

    # Report filtering stats
    RAW_READS=$(zcat "${FASTQ}" | awk 'END{print NR/4}')   # read count before filtering (lines/4)
    FILT_READS=$(zcat "${FILTERED}" | awk 'END{print NR/4}')  # read count after filtering
    RETAINED=$(echo "scale=1; ${FILT_READS}*100/${RAW_READS}" | bc)  # percent of reads retained
    log_msg "${SRR}: ${RAW_READS} raw -> ${FILT_READS} filtered (${RETAINED}% retained)"  # log the retention summary

    log_msg "Filtering completed for ${SRR}"
done

##==========================================================================##
##               STEP 5: POST-FILTERING QC                                   ##
##==========================================================================##

log_msg "========== STEP 5: Running post-filtering QC =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do                # QC each sample's filtered reads
    FILTERED="${FILTERED_DIR}/${SRR}_filtered.fastq.gz"  # this sample's filtered reads

    if [ ! -s "${FILTERED}" ] || ! is_valid_gz "${FILTERED}"; then  # skip missing/corrupt input
        log_msg "WARNING: ${FILTERED} not found or corrupted, skipping post-QC for ${SRR}."
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
        2>&1 | tee "${LOG_DIR}/${SRR}_nanoplot_post.log"  # length/quality plots of the filtered reads (for before/after)
    check_exit "NanoPlot (post-filter) failed for ${SRR}"  # abort on NanoPlot failure

    # --- NanoStat: Post-filtering summary ---
    log_msg "Running NanoStat (post-filter) on ${SRR}..."
    NanoStat \
        --fastq "${FILTERED}" \
        --outdir "${NANOSTAT_DIR}" \
        --name "${SRR}_post_stats.txt" \
        --threads "${THREADS}" \
        2>&1 | tee "${LOG_DIR}/${SRR}_nanostat_post.log"  # text summary stats of the filtered reads
    check_exit "NanoStat (post-filter) failed for ${SRR}"  # abort on NanoStat failure

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
    2>&1 | tee "${LOG_DIR}/multiqc.log"              # aggregate all per-sample QC into one report; --force overwrites a prior run
check_exit "MultiQC failed"                          # abort if aggregation failed

log_msg "MultiQC report generated: ${MULTIQC_DIR}/nanopore_qc_report.html"

##==========================================================================##
##               STEP 7: REFERENCE MAPPING (MINIMAP2 & SAMTOOLS)             ##
##==========================================================================##

log_msg "========== STEP 7: Reference Mapping with Minimap2 =========="

REF_ACC="K03455.1"                                   # HXB2 reference accession (the HIV-1 coordinate standard)
REF_FASTA="${REF_DIR}/${REF_ACC}.fasta"              # local path to the reference FASTA

mkdir -p "${ALIGN_DIR}/bam" "${REF_DIR}" "${CHECKPOINT_DIR}"  # ensure alignment/reference/checkpoint dirs exist

if [ ! -s "${REF_FASTA}" ]; then                     # download + index the reference only if it's not already present
    log_msg "Downloading HXB2 reference (${REF_ACC})..."
    wget -q -O "${REF_FASTA}" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${REF_ACC}&rettype=fasta&retmode=text"  # fetch the HXB2 FASTA from NCBI
    verify_nonempty "${REF_FASTA}" "HXB2 reference download"  # ensure the download produced content
    check_exit "HXB2 reference download failed"      # abort if the download failed
    samtools faidx "${REF_FASTA}"                    # build the .fai index needed by bcftools (minimap2 indexes on the fly)
    check_exit "samtools faidx failed on HXB2 reference"  # abort on faidx failure
fi

MAPPED_COUNT=0                                       # running count of successfully mapped samples
TOTAL_TO_MAP=0                                       # running count of samples that had input to map

for SRR in "${SRR_ACCESSIONS[@]}"; do                # map each sample to HXB2 and call a consensus
    FILTERED="${FILTERED_DIR}/${SRR}_filtered.fastq.gz"  # filtered reads to align
    BAM_OUT="${ALIGN_DIR}/bam/${SRR}.sorted.bam"     # sorted alignment output
    VCF_OUT="${ALIGN_DIR}/bam/${SRR}.vcf.gz"         # variant calls used to build the consensus
    CONSENSUS_OUT="${ALIGN_DIR}/${SRR}_consensus.fasta"  # per-sample consensus sequence
    CHKPT="${CHECKPOINT_DIR}/${SRR}_minimap2.done"   # resume marker for this sample's mapping

    if [ ! -s "${FILTERED}" ] || ! is_valid_gz "${FILTERED}"; then  # no valid input -> nothing to map
        continue
    fi

    TOTAL_TO_MAP=$((TOTAL_TO_MAP + 1))               # this sample has input, so count it toward the total

    if is_step_done "${CHKPT}" "${BAM_OUT}" "${CONSENSUS_OUT}"; then  # checkpoint + valid BAM/consensus means already mapped
        log_msg "Mapping already completed for ${SRR}, skipping."
        MAPPED_COUNT=$((MAPPED_COUNT + 1))           # count it as already-done
        continue
    fi

    log_msg "Aligning ${SRR} to HXB2..."
    minimap2 -ax map-ont -t "${THREADS}" "${REF_FASTA}" "${FILTERED}" 2> "${LOG_DIR}/${SRR}_minimap2.log" | \
        samtools view -@ "${THREADS}" -b - | \
        samtools sort -@ "${THREADS}" -o "${BAM_OUT}"  # long-read (map-ont) align, to BAM, coordinate-sort in one pipe
    if [ $? -ne 0 ] || ! verify_nonempty "${BAM_OUT}" "alignment for ${SRR}"; then  # mapping failed or produced no BAM
        log_msg "WARNING: mapping failed for ${SRR}, skipping (not marking done). See ${LOG_DIR}/${SRR}_minimap2.log"
        continue                                     # skip without a checkpoint so it retries next run
    fi

    samtools index "${BAM_OUT}"                      # build the BAM index needed by downstream tools
    if [ $? -ne 0 ]; then                            # if indexing failed...
        log_msg "WARNING: samtools index failed for ${SRR}, skipping (not marking done)."
        continue                                     # ...skip without marking done
    fi

    # Generate consensus sequence for downstream steps
    bcftools mpileup -Ou -f "${REF_FASTA}" "${BAM_OUT}" 2> "${LOG_DIR}/${SRR}_bcftools.log" | \
        bcftools call -c -Oz -o "${VCF_OUT}" 2>> "${LOG_DIR}/${SRR}_bcftools.log"  # pile up reads and call variants against HXB2
    if [ $? -ne 0 ] || ! verify_nonempty "${VCF_OUT}" "VCF for ${SRR}"; then  # variant calling failed or produced no VCF
        log_msg "WARNING: variant calling failed for ${SRR}, skipping (not marking done). See ${LOG_DIR}/${SRR}_bcftools.log"
        continue
    fi

    tabix -p vcf "${VCF_OUT}" 2>> "${LOG_DIR}/${SRR}_bcftools.log"  # index the VCF so bcftools consensus can use it
    if [ $? -ne 0 ]; then                            # if indexing the VCF failed...
        log_msg "WARNING: tabix indexing failed for ${SRR}, skipping (not marking done)."
        continue
    fi

    cat "${REF_FASTA}" | bcftools consensus "${VCF_OUT}" > "${CONSENSUS_OUT}" 2>> "${LOG_DIR}/${SRR}_bcftools.log"  # apply the sample's variants onto HXB2 to build its consensus
    if [ $? -ne 0 ] || ! verify_nonempty "${CONSENSUS_OUT}" "consensus for ${SRR}"; then  # consensus step failed or empty
        log_msg "WARNING: consensus generation failed for ${SRR}, skipping (not marking done). See ${LOG_DIR}/${SRR}_bcftools.log"
        continue
    fi

    # bcftools consensus does not rename the sequence header -- it keeps the
    # reference's own header verbatim, so every sample's consensus would
    # otherwise carry the identical ">${REF_ACC} ..." header. Left unfixed,
    # concatenating all samples for MSA (Step 9) would make every sequence
    # indistinguishable by ID. Rename the header to this sample's accession.
    sed -i "1s/.*/>${SRR}/" "${CONSENSUS_OUT}"       # rewrite the FASTA header to this sample's accession

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
    # poplars hypermut --input "${CONSENSUS}" --output "${FILTER_DIR}/${SRR}_poplars.fasta" || true
    # Rscript /path/to/HIVSeqinR/HIVSeqinR.R --input "${FILTER_DIR}/${SRR}_poplars.fasta" --outdir "${FILTER_DIR}/${SRR}_seqinr" || true

    touch "${CHKPT}"                                 # mark the (currently placeholder) step done
done

##==========================================================================##
##               STEP 9: MULTIPLE SEQUENCE ALIGNMENT                         ##
##==========================================================================##

log_msg "========== STEP 9: Multiple Sequence Alignment (MAFFT) =========="

COMBINED_FASTA="${MSA_DIR}/all_filtered_consensus.fasta"  # all consensuses (plus HXB2) concatenated for alignment
MSA_OUT="${MSA_DIR}/aligned_consensus.fasta"         # the resulting multiple sequence alignment
CHKPT="${CHECKPOINT_DIR}/mafft_alignment.done"       # resume marker for the alignment step

if ! is_step_done "${CHKPT}" "${MSA_OUT}"; then      # only align if not already done with valid output
    log_msg "Combining consensus sequences (with HXB2 as coordinate anchor) and aligning with MAFFT L-INS-i..."
    # HXB2 is prepended so it shares the alignment's column space -- Step 11's
    # U3 extraction anchors on HXB2's own row to convert its real annotated
    # genome coordinates into alignment-column coordinates.
    cat "${REF_FASTA}" "${ALIGN_DIR}"/*_consensus.fasta > "${COMBINED_FASTA}" 2>/dev/null || true  # concatenate HXB2 + all consensuses (tolerate no matches)
    if [ -s "${COMBINED_FASTA}" ]; then              # only align if there's something to align
        mafft --localpair --maxiterate 1000 --thread "${THREADS}" "${COMBINED_FASTA}" > "${MSA_OUT}" 2> "${LOG_DIR}/mafft.log"  # L-INS-i high-accuracy alignment
        if verify_nonempty "${MSA_OUT}" "MAFFT alignment"; then  # only checkpoint if alignment produced output
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
U3_GAPPED="${MOTIF_DIR}/U3_aligned.fasta"            # extracted U3 columns still in alignment (gapped) form
U3_EXTRACTED="${MOTIF_DIR}/U3_extracted.fasta"       # final ungapped per-sample U3 sequences
U3_WARNINGS="${MOTIF_DIR}/u3_extraction_warnings.log"  # outlier-check warnings from the extractor

if ! is_step_done "${U3_CHKPT}" "${U3_EXTRACTED}"; then  # only extract if not already done with valid output
    if [ -s "${MSA_OUT}" ]; then                     # can only extract U3 if the MSA exists
        log_msg "Extracting 5' LTR U3 region, anchored to HXB2's real annotated coordinates..."
        bash "${BASE_DIR}/scripts/utils/extract_u3_by_hxb2_anchor.sh" \
            --alignment "${MSA_OUT}" \
            --hxb2-id "${REF_ACC}" \
            --gb-cache "${REF_DIR}/${REF_ACC}.gb" \
            --out-gapped "${U3_GAPPED}" \
            --out "${U3_EXTRACTED}" \
            --warnings-log "${U3_WARNINGS}" \
            > "${LOG_DIR}/u3_extraction.log" 2>&1    # extract U3 by mapping HXB2's annotated coords to alignment columns
        if verify_nonempty "${U3_EXTRACTED}" "U3 extraction"; then  # only checkpoint if extraction produced output
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

if [ ! -f "${CHKPT}" ] && [ -s "${U3_EXTRACTED}" ]; then  # only run once, and only if U3 sequences exist
    log_msg "Motif scanning and G-quadruplex prediction are pending tool_comparison results (not yet wired in)."
    # 1. Motif scanning with FIMO (requires JASPAR database)
    # fimo --oc "${MOTIF_DIR}/fimo_out" JASPAR2024_CORE_vertebrates_non-redundant_pfms.meme "${U3_EXTRACTED}" || true

    # 2. G-quadruplex prediction (pqsfinder/gquad in R)
    # Rscript scripts/predict_g4.R "${U3_EXTRACTED}" "${MOTIF_DIR}" || true
fi

##==========================================================================##
##                          PIPELINE COMPLETE                                ##
##==========================================================================##

log_msg "=========================================="
log_msg "  NANOPORE END-TO-END PIPELINE COMPLETE   "
log_msg "=========================================="
log_msg "BioProject:        ${BIOPROJECT}"
log_msg "Samples processed: ${#SRR_ACCESSIONS[@]}"
log_msg "Raw data:          ${RAW_DIR}"
log_msg "Alignments:        ${ALIGN_DIR}"
log_msg "MSA:               ${MSA_DIR}"
log_msg "Motifs:            ${MOTIF_DIR}"
log_msg "=========================================="
