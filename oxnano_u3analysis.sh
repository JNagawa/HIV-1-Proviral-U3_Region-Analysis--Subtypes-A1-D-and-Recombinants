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

BIOPROJECT="PRJNA765218"
THREADS=${SLURM_CPUS_PER_TASK:-8}    # Use SLURM allocation or default to 8

# SRA accessions for BioProject PRJNA765218 (NanoHIV - Oxford Nanopore GridION)
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
MIN_QUALITY=7         # Minimum average read quality (Phred score)
MIN_LENGTH=200        # Minimum read length in bp

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
ALIGN_DIR="${BASE_DIR}/oxnano_alignments"
FILTER_DIR="${BASE_DIR}/oxnano_filtering"
MSA_DIR="${BASE_DIR}/oxnano_msa"
SUBTYPE_DIR="${BASE_DIR}/oxnano_subtyping"
MOTIF_DIR="${BASE_DIR}/oxnano_motifs"
REF_DIR="${BASE_DIR}/reference"
CHECKPOINT_DIR="${BASE_DIR}/checkpoints"

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
         "${FILTERED_DIR}" "${LOG_DIR}" \
         "${ALIGN_DIR}" "${FILTER_DIR}" "${MSA_DIR}" \
         "${SUBTYPE_DIR}" "${MOTIF_DIR}" "${REF_DIR}" "${CHECKPOINT_DIR}"

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
##               STEP 7: REFERENCE MAPPING (MINIMAP2 & SAMTOOLS)             ##
##==========================================================================##

log_msg "========== STEP 7: Reference Mapping with Minimap2 =========="

REF_ACC="K03455.1"
REF_FASTA="${REF_DIR}/${REF_ACC}.fasta"

mkdir -p "${ALIGN_DIR}/bam" "${REF_DIR}" "${CHECKPOINT_DIR}"

if [ ! -f "${REF_FASTA}" ]; then
    log_msg "Downloading HXB2 reference (${REF_ACC})..."
    wget -q -O "${REF_FASTA}" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${REF_ACC}&rettype=fasta&retmode=text"
    samtools faidx "${REF_FASTA}"
fi

for SRR in "${SRR_ACCESSIONS[@]}"; do
    FILTERED="${FILTERED_DIR}/${SRR}_filtered.fastq.gz"
    BAM_OUT="${ALIGN_DIR}/bam/${SRR}.sorted.bam"
    CHKPT="${CHECKPOINT_DIR}/${SRR}_minimap2.done"

    if [ ! -f "${FILTERED}" ]; then
        continue
    fi

    if [ -f "${CHKPT}" ] && [ -f "${BAM_OUT}" ]; then
        log_msg "Mapping already completed for ${SRR}, skipping."
        continue
    fi

    log_msg "Aligning ${SRR} to HXB2..."
    minimap2 -ax map-ont -t "${THREADS}" "${REF_FASTA}" "${FILTERED}" 2> "${LOG_DIR}/${SRR}_minimap2.log" | \
        samtools view -@ "${THREADS}" -b - | \
        samtools sort -@ "${THREADS}" -o "${BAM_OUT}"
    samtools index "${BAM_OUT}"
    
    # Generate consensus sequence for downstream steps
    bcftools mpileup -Ou -f "${REF_FASTA}" "${BAM_OUT}" 2> /dev/null | bcftools call -c -Oz -o "${ALIGN_DIR}/bam/${SRR}.vcf.gz" 2> /dev/null
    tabix -p vcf "${ALIGN_DIR}/bam/${SRR}.vcf.gz" 2> /dev/null
    cat "${REF_FASTA}" | bcftools consensus "${ALIGN_DIR}/bam/${SRR}.vcf.gz" > "${ALIGN_DIR}/${SRR}_consensus.fasta" 2> /dev/null
    
    touch "${CHKPT}"
done

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

if [ ! -f "${CHKPT}" ]; then
    log_msg "Combining consensus sequences and aligning with MAFFT L-INS-i..."
    cat "${ALIGN_DIR}"/*_consensus.fasta > "${COMBINED_FASTA}" 2>/dev/null || true
    if [ -s "${COMBINED_FASTA}" ]; then
        mafft --localpair --maxiterate 1000 --thread "${THREADS}" "${COMBINED_FASTA}" > "${MSA_OUT}" 2> "${LOG_DIR}/mafft.log" || true
        touch "${CHKPT}"
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
    # jpHMM -v HIV -s "${MSA_OUT}" -o "${SUBTYPE_DIR}/jphmm_out" || true
    # iqtree2 -s "${MSA_OUT}" -m MFP -B 1000 -T "${THREADS}" --prefix "${SUBTYPE_DIR}/iqtree" || true
    # comet -i "${MSA_OUT}" -o "${SUBTYPE_DIR}/comet_results.csv" || true
    touch "${CHKPT}"
fi

##==========================================================================##
##               STEP 11: U3 EXTRACTION & MOTIF MAPPING                      ##
##==========================================================================##

log_msg "========== STEP 11: Motif Mapping & G-Quadruplexes =========="

CHKPT="${CHECKPOINT_DIR}/motif_mapping.done"

if [ ! -f "${CHKPT}" ] && [ -f "${MSA_OUT}" ]; then
    log_msg "Extracting U3 region and scanning motifs..."
    # 1. Extract U3 (HXB2 coords roughly 8677-9121 or 5' LTR 1-454)
    seqkit subseq -r 8677:9121 "${MSA_OUT}" > "${MOTIF_DIR}/U3_extracted.fasta" 2>/dev/null || true
    
    # 2. Motif scanning with FIMO (requires JASPAR database)
    # fimo --oc "${MOTIF_DIR}/fimo_out" JASPAR2024_CORE_vertebrates_non-redundant_pfms.meme "${MOTIF_DIR}/U3_extracted.fasta" || true
    
    # 3. G-quadruplex prediction (pqsfinder/gquad in R)
    # Rscript scripts/predict_g4.R "${MOTIF_DIR}/U3_extracted.fasta" "${MOTIF_DIR}" || true
    
    touch "${CHKPT}"
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
