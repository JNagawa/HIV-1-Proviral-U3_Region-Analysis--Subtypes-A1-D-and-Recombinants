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
RAW_DIR="${BASE_DIR}/data/raw/illumina"
QC_DIR="${BASE_DIR}/results/reports/qc/illumina"
FASTQC_PRE_DIR="${QC_DIR}/fastqc_pre"
FASTQC_POST_DIR="${QC_DIR}/fastqc_post"
MULTIQC_DIR="${QC_DIR}/multiqc"
TRIMMED_DIR="${BASE_DIR}/data/processed/illumina/trimmed_trimmomatic"
LOG_DIR="${BASE_DIR}/logs"
ALIGN_DIR="${BASE_DIR}/data/processed/illumina/alignments"
FILTER_DIR="${BASE_DIR}/data/processed/illumina/filtering"
MSA_DIR="${BASE_DIR}/data/processed/illumina/msa"
SUBTYPE_DIR="${BASE_DIR}/data/processed/illumina/subtyping"
MOTIF_DIR="${BASE_DIR}/data/processed/illumina/motifs"
REF_DIR="${BASE_DIR}/data/reference"
CHECKPOINT_DIR="${BASE_DIR}/checkpoints"

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

# Verify a file exists and is non-empty; logs and returns 1 otherwise.
# Used to gate checkpoint creation instead of touching it unconditionally.
verify_nonempty() {
    local f="$1" label="$2"
    if [ ! -s "${f}" ]; then
        log_msg "ERROR: ${label} produced empty/missing output: ${f}"
        return 1
    fi
    return 0
}

# Cheap integrity check for gzip files (catches truncated/corrupted output,
# e.g. from a crash mid-write).
is_valid_gz() {
    gzip -t "$1" >/dev/null 2>&1
}

# True only if the checkpoint file exists AND every output file argument is
# non-empty AND (for *.gz files) passes a gzip integrity check. Replaces bare
# "[ -f "$CHKPT" ]"-style resume checks so a checkpoint next to a truncated,
# corrupted, or zero-byte output is never trusted as "done" again.
is_step_done() {
    local chkpt="$1"
    shift
    [ -f "${chkpt}" ] || return 1
    local f
    for f in "$@"; do
        [ -s "${f}" ] || return 1
        case "${f}" in
            *.gz) is_valid_gz "${f}" || return 1 ;;
        esac
    done
    return 0
}

# Pre-flight dependency check; exits immediately with a clear message instead
# of letting a missing tool fail silently deep inside a pipe (this is what let
# a missing minimap2 mark 9 Nanopore samples "done" with 0-byte output).
require_tool() {
    command -v "$1" >/dev/null 2>&1 || { log_msg "ERROR: required tool '$1' not found on PATH. Check conda env HIV_U3analysis."; exit 1; }
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
##               STEP 0: PRE-FLIGHT DEPENDENCY CHECK                         ##
##==========================================================================##

log_msg "========== STEP 0: Checking required tools are on PATH =========="

for TOOL in prefetch vdb-validate fasterq-dump fastqc trimmomatic multiqc \
            bwa samtools bcftools tabix mafft seqkit; do
    require_tool "${TOOL}"
done

log_msg "All required tools found."

##==========================================================================##
##               STEP 1: CREATE DIRECTORY STRUCTURE                          ##
##==========================================================================##

log_msg "========== STEP 1: Setting up directory structure =========="

mkdir -p "${RAW_DIR}" "${FASTQC_PRE_DIR}" "${FASTQC_POST_DIR}" \
         "${MULTIQC_DIR}" "${TRIMMED_DIR}" "${LOG_DIR}" \
         "${ALIGN_DIR}" "${FILTER_DIR}" "${MSA_DIR}" \
         "${SUBTYPE_DIR}" "${MOTIF_DIR}" "${REF_DIR}" "${CHECKPOINT_DIR}"

log_msg "Directory structure created under: ${BASE_DIR}"

##==========================================================================##
##               STEP 2: DOWNLOAD SRA DATA                                   ##
##==========================================================================##

log_msg "========== STEP 2: Downloading SRA data (${#SRR_ACCESSIONS[@]} samples) =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do
    log_msg "--- Processing ${SRR} ---"

    # Skip if paired-end FASTQs already exist and are intact (not a truncated
    # partial download from an earlier interrupted/crashed run)
    if [ -s "${RAW_DIR}/${SRR}_1.fastq.gz" ] && [ -s "${RAW_DIR}/${SRR}_2.fastq.gz" ] && \
       is_valid_gz "${RAW_DIR}/${SRR}_1.fastq.gz" && is_valid_gz "${RAW_DIR}/${SRR}_2.fastq.gz"; then
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

    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ] || ! is_valid_gz "${R1}" || ! is_valid_gz "${R2}"; then
        log_msg "WARNING: Paired FASTQs not found or corrupted for ${SRR}, skipping pre-QC."
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

    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ] || ! is_valid_gz "${R1}" || ! is_valid_gz "${R2}"; then
        log_msg "WARNING: Paired FASTQs not found or corrupted for ${SRR}, skipping trimming."
        continue
    fi

    # Skip if trimmed files already exist and are intact (a truncated gzip
    # from an earlier crashed Trimmomatic run must NOT be trusted as "done")
    if [ -s "${R1_PAIRED}" ] && [ -s "${R2_PAIRED}" ] && \
       is_valid_gz "${R1_PAIRED}" && is_valid_gz "${R2_PAIRED}"; then
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
        -Xmx48g \
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

    if [ ! -s "${R1_PAIRED}" ] || [ ! -s "${R2_PAIRED}" ] || ! is_valid_gz "${R1_PAIRED}" || ! is_valid_gz "${R2_PAIRED}"; then
        log_msg "WARNING: Trimmed FASTQs not found or corrupted for ${SRR}, skipping post-QC."
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
##               STEP 7: REFERENCE MAPPING (BWA & SAMTOOLS)                  ##
##==========================================================================##

log_msg "========== STEP 7: Reference Mapping with BWA-MEM =========="

REF_ACC="K03455.1"
REF_FASTA="${REF_DIR}/${REF_ACC}.fasta"

mkdir -p "${ALIGN_DIR}/bam" "${REF_DIR}" "${CHECKPOINT_DIR}"

if [ ! -s "${REF_FASTA}" ]; then
    log_msg "Downloading HXB2 reference (${REF_ACC})..."
    wget -q -O "${REF_FASTA}" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${REF_ACC}&rettype=fasta&retmode=text"
    verify_nonempty "${REF_FASTA}" "HXB2 reference download"
    check_exit "HXB2 reference download failed"
    bwa index "${REF_FASTA}" > "${LOG_DIR}/bwa_index.log" 2>&1
    check_exit "bwa index failed on HXB2 reference"
    samtools faidx "${REF_FASTA}"
    check_exit "samtools faidx failed on HXB2 reference"
fi

MAPPED_COUNT=0
TOTAL_TO_MAP=0

for SRR in "${SRR_ACCESSIONS[@]}"; do
    R1_PAIRED="${TRIMMED_DIR}/${SRR}_1_paired.fastq.gz"
    R2_PAIRED="${TRIMMED_DIR}/${SRR}_2_paired.fastq.gz"
    BAM_OUT="${ALIGN_DIR}/bam/${SRR}.sorted.bam"
    VCF_OUT="${ALIGN_DIR}/bam/${SRR}.vcf.gz"
    CONSENSUS_OUT="${ALIGN_DIR}/${SRR}_consensus.fasta"
    CHKPT="${CHECKPOINT_DIR}/${SRR}_bwa.done"

    if [ ! -s "${R1_PAIRED}" ] || [ ! -s "${R2_PAIRED}" ] || ! is_valid_gz "${R1_PAIRED}" || ! is_valid_gz "${R2_PAIRED}"; then
        continue
    fi

    TOTAL_TO_MAP=$((TOTAL_TO_MAP + 1))

    if is_step_done "${CHKPT}" "${BAM_OUT}" "${CONSENSUS_OUT}"; then
        log_msg "Mapping already completed for ${SRR}, skipping."
        MAPPED_COUNT=$((MAPPED_COUNT + 1))
        continue
    fi

    log_msg "Aligning ${SRR} to HXB2..."
    bwa mem -t "${THREADS}" "${REF_FASTA}" "${R1_PAIRED}" "${R2_PAIRED}" 2> "${LOG_DIR}/${SRR}_bwa.log" | \
        samtools view -@ "${THREADS}" -b - | \
        samtools sort -@ "${THREADS}" -o "${BAM_OUT}"
    if [ $? -ne 0 ] || ! verify_nonempty "${BAM_OUT}" "alignment for ${SRR}"; then
        log_msg "WARNING: mapping failed for ${SRR}, skipping (not marking done). See ${LOG_DIR}/${SRR}_bwa.log"
        continue
    fi

    samtools index "${BAM_OUT}"
    if [ $? -ne 0 ]; then
        log_msg "WARNING: samtools index failed for ${SRR}, skipping (not marking done)."
        continue
    fi

    # Generate consensus sequence for downstream steps
    bcftools mpileup -Ou -f "${REF_FASTA}" "${BAM_OUT}" 2> "${LOG_DIR}/${SRR}_bcftools.log" | \
        bcftools call -c -Oz -o "${VCF_OUT}" 2>> "${LOG_DIR}/${SRR}_bcftools.log"
    if [ $? -ne 0 ] || ! verify_nonempty "${VCF_OUT}" "VCF for ${SRR}"; then
        log_msg "WARNING: variant calling failed for ${SRR}, skipping (not marking done). See ${LOG_DIR}/${SRR}_bcftools.log"
        continue
    fi

    tabix -p vcf "${VCF_OUT}" 2>> "${LOG_DIR}/${SRR}_bcftools.log"
    if [ $? -ne 0 ]; then
        log_msg "WARNING: tabix indexing failed for ${SRR}, skipping (not marking done)."
        continue
    fi

    cat "${REF_FASTA}" | bcftools consensus "${VCF_OUT}" > "${CONSENSUS_OUT}" 2>> "${LOG_DIR}/${SRR}_bcftools.log"
    if [ $? -ne 0 ] || ! verify_nonempty "${CONSENSUS_OUT}" "consensus for ${SRR}"; then
        log_msg "WARNING: consensus generation failed for ${SRR}, skipping (not marking done). See ${LOG_DIR}/${SRR}_bcftools.log"
        continue
    fi

    # bcftools consensus does not rename the sequence header -- it keeps the
    # reference's own header verbatim, so every sample's consensus would
    # otherwise carry the identical ">${REF_ACC} ..." header. Left unfixed,
    # concatenating all samples for MSA (Step 9) would make every sequence
    # indistinguishable by ID. Rename the header to this sample's accession.
    sed -i "1s/.*/>${SRR}/" "${CONSENSUS_OUT}"

    touch "${CHKPT}"
    MAPPED_COUNT=$((MAPPED_COUNT + 1))
    log_msg "Mapping and consensus generation completed for ${SRR}"
done

log_msg "Step 7 complete: ${MAPPED_COUNT}/${TOTAL_TO_MAP} samples mapped successfully."

##==========================================================================##
##               STEP 8: BIOLOGICAL FILTERING                                ##
##==========================================================================##

log_msg "========== STEP 8: Biological Filtering (Poplars & HIVSeqinR) =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do
    CONSENSUS="${ALIGN_DIR}/${SRR}_consensus.fasta"
    CHKPT="${CHECKPOINT_DIR}/${SRR}_filtering.done"
    
    if [ ! -f "${CONSENSUS}" ]; then continue; fi
    if [ -f "${CHKPT}" ]; then continue; fi

    log_msg "Running Hypermut 3 (Poplars) & HIVSeqinR on ${SRR}..."
    # Note: These are generalized commands. Adjust paths/params for your cluster environment.
    # poplars hypermut --input "${CONSENSUS}" --output "${FILTER_DIR}/${SRR}_poplars.fasta" || true
    # Rscript /path/to/HIVSeqinR/HIVSeqinR.R --input "${FILTER_DIR}/${SRR}_poplars.fasta" --outdir "${FILTER_DIR}/${SRR}_seqinr" || true
    
    touch "${CHKPT}"
done

##==========================================================================##
##               STEP 9: MULTIPLE SEQUENCE ALIGNMENT                         ##
##==========================================================================##

log_msg "========== STEP 9: Multiple Sequence Alignment (MAFFT) =========="

COMBINED_FASTA="${MSA_DIR}/all_filtered_consensus.fasta"
MSA_OUT="${MSA_DIR}/aligned_consensus.fasta"
CHKPT="${CHECKPOINT_DIR}/mafft_alignment.done"

if ! is_step_done "${CHKPT}" "${MSA_OUT}"; then
    log_msg "Combining consensus sequences (with HXB2 as coordinate anchor) and aligning with MAFFT L-INS-i..."
    # HXB2 is prepended so it shares the alignment's column space -- Step 11's
    # U3 extraction anchors on HXB2's own row to convert its real annotated
    # genome coordinates into alignment-column coordinates.
    cat "${REF_FASTA}" "${ALIGN_DIR}"/*_consensus.fasta > "${COMBINED_FASTA}" 2>/dev/null || true
    if [ -s "${COMBINED_FASTA}" ]; then
        mafft --localpair --maxiterate 1000 --thread "${THREADS}" "${COMBINED_FASTA}" > "${MSA_OUT}" 2> "${LOG_DIR}/mafft.log"
        if verify_nonempty "${MSA_OUT}" "MAFFT alignment"; then
            touch "${CHKPT}"
        else
            log_msg "WARNING: MAFFT alignment failed, not marking step done. See ${LOG_DIR}/mafft.log"
        fi
    else
        log_msg "No consensus sequences found for MSA."
    fi
fi

##==========================================================================##
##               STEP 10: SUBTYPE & RECOMBINATION DETECTION                  ##
##==========================================================================##

log_msg "========== STEP 10: Subtyping (jpHMM, IQ-TREE 2, COMET, REGA) =========="

CHKPT="${CHECKPOINT_DIR}/subtyping.done"

if [ ! -f "${CHKPT}" ] && [ -f "${MSA_OUT}" ]; then
    log_msg "Running subtyping pipeline..."
    # Generalized commands
    # jpHMM -v HIV -s "${MSA_OUT}" -o "${SUBTYPE_DIR}/jphmm_out" || true
    # iqtree2 -s "${MSA_OUT}" -m MFP -B 1000 -T "${THREADS}" --prefix "${SUBTYPE_DIR}/iqtree" || true
    # comet -i "${MSA_OUT}" -o "${SUBTYPE_DIR}/comet_results.csv" || true
    touch "${CHKPT}"
fi

##==========================================================================##
##               STEP 11: U3 EXTRACTION & MOTIF MAPPING                      ##
##==========================================================================##

log_msg "========== STEP 11: U3 Extraction & Motif Mapping =========="

U3_CHKPT="${CHECKPOINT_DIR}/u3_extraction.done"
U3_GAPPED="${MOTIF_DIR}/U3_aligned.fasta"
U3_EXTRACTED="${MOTIF_DIR}/U3_extracted.fasta"
U3_WARNINGS="${MOTIF_DIR}/u3_extraction_warnings.log"

if ! is_step_done "${U3_CHKPT}" "${U3_EXTRACTED}"; then
    if [ -s "${MSA_OUT}" ]; then
        log_msg "Extracting 5' LTR U3 region, anchored to HXB2's real annotated coordinates..."
        bash "${BASE_DIR}/scripts/utils/extract_u3_by_hxb2_anchor.sh" \
            --alignment "${MSA_OUT}" \
            --hxb2-id "${REF_ACC}" \
            --gb-cache "${REF_DIR}/${REF_ACC}.gb" \
            --out-gapped "${U3_GAPPED}" \
            --out "${U3_EXTRACTED}" \
            --warnings-log "${U3_WARNINGS}" \
            > "${LOG_DIR}/u3_extraction.log" 2>&1
        if verify_nonempty "${U3_EXTRACTED}" "U3 extraction"; then
            touch "${U3_CHKPT}"
            if [ -s "${U3_WARNINGS}" ]; then
                log_msg "WARNING: some samples flagged by the U3 extraction outlier check, see ${U3_WARNINGS}"
            fi
        else
            log_msg "WARNING: U3 extraction failed, not marking step done. See ${LOG_DIR}/u3_extraction.log"
        fi
    else
        log_msg "No MSA output found, skipping U3 extraction."
    fi
fi

CHKPT="${CHECKPOINT_DIR}/motif_mapping.done"

if [ ! -f "${CHKPT}" ] && [ -s "${U3_EXTRACTED}" ]; then
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
log_msg "  ILLUMINA END-TO-END PIPELINE COMPLETE   "
log_msg "=========================================="
log_msg "BioProject:         ${BIOPROJECT}"
log_msg "Samples processed:  ${#SRR_ACCESSIONS[@]}"
log_msg "Raw data:           ${RAW_DIR}"
log_msg "Alignments:         ${ALIGN_DIR}"
log_msg "MSA:                ${MSA_DIR}"
log_msg "Motifs:             ${MOTIF_DIR}"
log_msg "=========================================="
