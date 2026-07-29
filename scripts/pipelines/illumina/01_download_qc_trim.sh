#!/bin/bash
#SBATCH --job-name=illumina_step1_download_qc_trim
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --time=14:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# Exit the script if any command in a pipe (e.g. "cmd1 | cmd2") fails, not
# just the last one -- otherwise a failed prefetch/fastqc/etc. piped into
# `tee` would be masked by tee's own (near-always successful) exit code.
set -o pipefail

##-------DESCRIPTION--------##
## Step 1 of the (in-progress) per-step Illumina pipeline breakdown -- see
## scripts/pipelines/illumina_u3analysis.sh for the still-current monolithic
## end-to-end script; this replaces that script's steps 1-6 (dir setup,
## download, pre-trim FastQC, trimming, post-trim FastQC, MultiQC) and adds
## fastp as a second, real trimming track alongside Trimmomatic (previously
## fastp was only used for a side-by-side comparison in
## scripts/download_qc_illumina/, never as production output).
##
## STUDY: BioProject PRJNA207834 (HIV-1 Uganda A1/D, Illumina MiSeq 2x251bp PE)
##
## Produces, per sample:
##   - raw paired FASTQs                          (data/raw/illumina/)
##   - pre-trim FastQC                             (results/reports/qc/illumina/fastqc_pre/)
##   - Trimmomatic-trimmed reads + summary         (data/processed/illumina/trimmed_trimmomatic/)
##   - fastp-trimmed reads + JSON/HTML report      (data/processed/illumina/trimmed_fastp/)
##   - post-trim FastQC on the Trimmomatic output  (results/reports/qc/illumina/fastqc_post/)
## Then two MultiQC reports:
##   - fastqc_trimmomatic_report.html  (fastqc_pre + fastqc_post + Trimmomatic logs)
##   - fastp_report.html               (fastp's own JSON reports)
##
## QUALITY-THRESHOLD RATIONALE (see also the CONFIGURATION section below):
## Phred quality score Q is defined as Q = -10*log10(P_error) (Ewing & Green
## 1998) -- Q15 corresponds to a 1-in-32 (~97%) per-base call accuracy, Q20
## to 1-in-100 (99%). Q15 is used throughout below (sliding-window trim AND
## whole-read average), matching Trimmomatic's own textbook
## SLIDINGWINDOW:4:15 example.
##
## NOTE: Ensure this script is executed on a system with Anaconda, Miniconda or
## an equivalent Python distribution. Create an environment using the provided
## .yml file to ensure dependencies are installed before running this script.

# AUTHOR DETAILS
# No.   NAME                       STUDENT NO.          REG NO.
# 1.    Nagawa Jovita              2500726007           2025/HD07/26007

# ---------------------------------------------------------------------------
# Activate the conda environment that has all the tools this script needs
# (prefetch/fasterq-dump, fastqc, trimmomatic, fastp, multiqc). Tries the
# usual miniconda3 install location first, falls back to whatever `conda`
# is already on PATH, and refuses to continue silently if neither works.
# ---------------------------------------------------------------------------
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

BIOPROJECT="PRJNA207834"                             # SRA BioProject this step processes
# Use however many CPUs Slurm actually allocated this job; if run outside
# Slurm (no SLURM_CPUS_PER_TASK set), fall back to 8.
# thread count for every multi-threaded tool below
THREADS=${SLURM_CPUS_PER_TASK:-8}

# SRA accessions for BioProject PRJNA207834
# 24 HIV-1 near-full-genome Illumina MiSeq paired-end samples from Uganda
# Subtypes: A1, D, and A1-D intersubtype recombinants (BSRI)
SRR_ACCESSIONS=(                                     # the 24 samples this step downloads and trims
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

# --- Trimming thresholds shared by BOTH trimmers (see rationale above) ---
# LEADING/TRAILING:3 -- Trimmomatic's own manual example value: cut a base
# from the very start/end of a read the moment its quality drops below Q3.
# Q3 is barely above "no confidence at all" -- this is deliberately lenient,
# just removing the occasional genuinely-unusable edge base; the real
# quality bar is enforced by SLIDINGWINDOW and AVGQUAL below, not this.
# trim leading bases below Q3 (lenient edge cleanup)
TRIM_LEADING=3
# trim trailing bases below Q3 (lenient edge cleanup)
TRIM_TRAILING=3
# SLIDINGWINDOW:4:15 -- slide a 4-base window along the read; once the
# window's mean quality falls below Q15 (~97% accuracy), trim from there to
# the read's end. 4bp is Trimmomatic's own recommended window size; Q15 is
# Trimmomatic's own manual example threshold.
# 4bp window, trim once its mean quality drops below Q15
TRIM_SLIDINGWINDOW="4:15"
# MINLEN:50 -- discard a read/pair if trimming leaves it under 50bp. Reads
# are 251bp to start, and the reference (HXB2, K03455.1) is only ~9.7kb, so
# very short leftover fragments map ambiguously (multi-mapping) with
# BWA-MEM; 50bp keeps enough sequence for confident, unique placement while
# still retaining reads that only lost their tail to quality trimming.
# drop reads shorter than 50bp after trimming (avoids multi-mapping)
TRIM_MINLEN=50
# AVGQUAL:15 -- independent of the sliding window, also require the read's
# OVERALL average quality to be >=Q15 after trimming. A read can pass a
# 4bp sliding window check while still being poor quality on average; this
# is the final whole-read quality gate.
# whole-read mean-quality floor (Q15), applied to both trimmers
TRIM_AVGQUAL=15

# Directory structure (matches scripts/pipelines/illumina_u3analysis.sh for
# everything except the trimmed-reads directories, so this step's output is
# otherwise a drop-in replacement for that script's steps 1-6)
# repo root -- MUST be invoked from here (see README)
BASE_DIR="$(pwd)"
RAW_DIR="${BASE_DIR}/data/raw/illumina"                            # downloaded, untrimmed FASTQs
# all Illumina QC output lives under here
QC_DIR="${BASE_DIR}/results/reports/qc/illumina"
FASTQC_PRE_DIR="${QC_DIR}/fastqc_pre"                              # FastQC on raw reads
# FastQC on Trimmomatic's trimmed reads
FASTQC_POST_DIR="${QC_DIR}/fastqc_post"
# the two aggregated MultiQC reports
MULTIQC_DIR="${QC_DIR}/multiqc"
# Trimmomatic's trimmed reads
TRIMMED_TRIMMOMATIC_DIR="${BASE_DIR}/data/processed/illumina/trimmed_trimmomatic"
# fastp's trimmed reads
TRIMMED_FASTP_DIR="${BASE_DIR}/data/processed/illumina/trimmed_fastp"
# fastp --dedup-only pass on Trimmomatic's output
DEDUP_TRIMMOMATIC_DIR="${BASE_DIR}/data/processed/illumina/dedup_trimmomatic"
# Kraken2-filtered, Trimmomatic track
KRAKEN2_TRIMMOMATIC_DIR="${BASE_DIR}/data/processed/illumina/kraken2_trimmomatic"
# Kraken2-filtered, fastp track (production input to assembly)
KRAKEN2_FASTP_DIR="${BASE_DIR}/data/processed/illumina/kraken2_fastp"
# size-capped Standard DB (human+bacteria+archaea+viral)
KRAKEN2_DB="${BASE_DIR}/data/reference/kraken2_standard_16gb_db"
# per-sample tool logs (shared across all steps)
LOG_DIR="${BASE_DIR}/logs"
# resume-state markers (shared across all steps)
CHECKPOINT_DIR="${BASE_DIR}/checkpoints"

##==========================================================================##
##                         HELPER FUNCTIONS                                  ##
##==========================================================================##

# Timestamped log line, tagged so it's easy to grep this step's output out
# of a combined Slurm log.
log_msg() {
    # timestamp + step tag + the message argument
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ILLUMINA/STEP1] $1"
}

# Call immediately after a command: if that command's exit code ($?) was
# non-zero, log the given message and abort the whole script.
check_exit() {
    # inspect the exit status of the previous command
    if [ $? -ne 0 ]; then
        log_msg "ERROR: $1"                          # log the caller-supplied failure message
        # abort the pipeline (this step's failures are fatal)
        exit 1
    fi
}

# gzip -t verifies the archive's integrity without decompressing it to disk
# -- used everywhere below to distinguish a genuinely-complete output file
# from a truncated one left behind by a crashed/killed earlier run.
is_valid_gz() {
    # test gzip integrity silently; exit status is the answer
    gzip -t "$1" >/dev/null 2>&1
}

# Fail fast with a clear message if a required tool is missing, rather than
# discovering it 20 minutes into a run when a pipe silently produces
# nothing.
require_tool() {
    # abort now if the tool isn't on PATH
    command -v "$1" >/dev/null 2>&1 || { log_msg "ERROR: required tool '$1' not found on PATH. Check conda env HIV_U3analysis."; exit 1; }
}

# Trimmomatic needs an explicit path to its TruSeq adapter FASTA; conda
# installs it in different places depending on version, so probe the
# common locations rather than hardcoding one.
find_adapter_file() {
    # candidate paths for the TruSeq3-PE-2.fa adapter file
    local ADAPTER_LOCATIONS=(
        "${CONDA_PREFIX}/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
        "${CONDA_PREFIX}/share/trimmomatic-*/adapters/TruSeq3-PE-2.fa"
        "/usr/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
        "/usr/local/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
    )
    for pattern in "${ADAPTER_LOCATIONS[@]}"; do      # try each candidate in order
        local found                                  # holds the first matching path, if any
        # compgen -G expands the glob safely (empty result instead of an error
        # if nothing matches), so this loop never crashes on a missing path.
        # expand the glob, take the first hit
        found=$(compgen -G "${pattern}" 2>/dev/null | head -1)
        if [ -n "${found}" ] && [ -f "${found}" ]; then  # if we got a real, existing file...
            echo "${found}"                          # ...print it (the function's return value)...
            return 0                                 # ...and report success
        fi
    done
    # none found: warn but continue
    log_msg "WARNING: TruSeq3-PE-2.fa adapter file not found. Trimmomatic will skip adapter trimming."
    echo ""                                          # print empty so the caller's $() is empty
    return 1                                          # signal not-found to the caller
}

##==========================================================================##
##               STEP 0: PRE-FLIGHT DEPENDENCY CHECK                         ##
##==========================================================================##
# Check every tool this script calls is actually reachable BEFORE doing any
# work, so a missing tool fails immediately and clearly instead of partway
# through a 24-sample loop.

log_msg "========== Checking required tools are on PATH =========="

# every executable this step depends on
for TOOL in prefetch vdb-validate fasterq-dump fastqc trimmomatic fastp multiqc; do
    require_tool "${TOOL}"                           # abort immediately if any one is missing
done

log_msg "All required tools found."

##==========================================================================##
##               1a: DIRECTORY STRUCTURE                                     ##
##==========================================================================##
# Create every output directory this step writes to. mkdir -p is a no-op on
# directories that already exist, so this is always safe to re-run.

# create all output/log/checkpoint dirs in one call
mkdir -p "${RAW_DIR}" "${FASTQC_PRE_DIR}" "${FASTQC_POST_DIR}" "${MULTIQC_DIR}" \
         "${TRIMMED_TRIMMOMATIC_DIR}" "${TRIMMED_FASTP_DIR}" \
         "${DEDUP_TRIMMOMATIC_DIR}" "${KRAKEN2_TRIMMOMATIC_DIR}" "${KRAKEN2_FASTP_DIR}" \
         "${LOG_DIR}" "${CHECKPOINT_DIR}"

log_msg "Directory structure ready under: ${BASE_DIR}"

##==========================================================================##
##               1b: DOWNLOAD SRA DATA                                       ##
##==========================================================================##
# For each accession: skip if already downloaded intact, else fetch the .sra
# file, validate it, convert to paired FASTQ, compress, and clean up the
# large intermediate .sra cache.

log_msg "========== Downloading SRA data (${#SRR_ACCESSIONS[@]} samples) =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do                # process each accession in turn
    log_msg "--- Processing ${SRR} ---"

    # Resume support: if both paired FASTQs already exist AND pass a gzip
    # integrity check, this sample is done -- don't re-download it.
    # both mates present and intact?
    if [ -s "${RAW_DIR}/${SRR}_1.fastq.gz" ] && [ -s "${RAW_DIR}/${SRR}_2.fastq.gz" ] && \
       is_valid_gz "${RAW_DIR}/${SRR}_1.fastq.gz" && is_valid_gz "${RAW_DIR}/${SRR}_2.fastq.gz"; then
        log_msg "FASTQs for ${SRR} already exist, skipping download."  # already done
        continue                                     # skip to the next accession
    fi

    # Step 1: pull the .sra file itself from NCBI's SRA archive.
    log_msg "Prefetching ${SRR}..."
    # download the .sra, mirroring output to a per-sample log
    prefetch "${SRR}" \
        --output-directory "${RAW_DIR}" \
        --max-size 50G \
        --progress \
        2>&1 | tee "${LOG_DIR}/${SRR}_prefetch.log"
    check_exit "prefetch failed for ${SRR}"          # abort if the download failed

    # Step 2: verify the .sra download isn't corrupted/truncated before
    # spending time converting it; if it fails, force a clean re-download
    # once rather than converting a broken file.
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

    # Step 3: convert the validated .sra into paired-end FASTQ files.
    # --split-3 separates mates into _1/_2 files plus a third file for any
    # unpaired/orphan reads (rather than interleaving or dropping them).
    log_msg "Converting ${SRR} to paired-end FASTQs..."
    # extract paired FASTQs (_1/_2 + orphans) from the .sra
    fasterq-dump "${RAW_DIR}/${SRR}/${SRR}.sra" \
        --outdir "${RAW_DIR}" \
        --split-3 \
        --threads "${THREADS}" \
        --progress \
        2>&1 | tee "${LOG_DIR}/${SRR}_fasterq.log"
    check_exit "fasterq-dump failed for ${SRR}"      # abort if conversion failed

    # Step 4: fasterq-dump writes plain (uncompressed) FASTQ -- compress to
    # save disk space before moving on. -f overwrites without prompting;
    # 2>/dev/null on the orphan-reads file since it may legitimately not
    # exist for every sample.
    log_msg "Compressing FASTQs for ${SRR}..."
    gzip -f "${RAW_DIR}/${SRR}_1.fastq" 2>/dev/null  # compress the forward mate
    gzip -f "${RAW_DIR}/${SRR}_2.fastq" 2>/dev/null  # compress the reverse mate
    # compress orphan reads if present (silently skip if not)
    gzip -f "${RAW_DIR}/${SRR}.fastq" 2>/dev/null

    # Step 5: the raw .sra cache is no longer needed once FASTQs exist --
    # it's large and would otherwise sit around consuming disk space.
    rm -rf "${RAW_DIR}/${SRR}"                        # delete the .sra cache dir to reclaim space

    log_msg "Completed download for ${SRR}"
done

log_msg "All SRA downloads completed."

##==========================================================================##
##               1c: PRE-TRIMMING QC WITH FASTQC                             ##
##==========================================================================##
# Run FastQC on the untouched raw reads, so the pre/post-trim MultiQC report
# can show exactly what trimming changed (adapter content, per-base quality,
# etc.) for each sample.

log_msg "========== Running pre-trimming FastQC =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do                # QC each sample's raw reads
    R1="${RAW_DIR}/${SRR}_1.fastq.gz"                # forward-mate raw FASTQ
    R2="${RAW_DIR}/${SRR}_2.fastq.gz"                # reverse-mate raw FASTQ

    # Can't QC a sample whose download didn't succeed/is corrupted.
    # skip missing/corrupt input
    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ] || ! is_valid_gz "${R1}" || ! is_valid_gz "${R2}"; then
        log_msg "WARNING: Paired FASTQs not found or corrupted for ${SRR}, skipping pre-QC."
        continue
    fi

    # Resume support: FastQC's own HTML report existing means this sample
    # is already done.
    # both reports already present?
    if [ -s "${FASTQC_PRE_DIR}/${SRR}_1_fastqc.html" ] && [ -s "${FASTQC_PRE_DIR}/${SRR}_2_fastqc.html" ]; then
        log_msg "Pre-trim FastQC for ${SRR} already exists, skipping."
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
##               1d: TRIMMING -- TRIMMOMATIC                                 ##
##==========================================================================##
# Adapter, leading/trailing, sliding-window, and average-quality trimming
# with Trimmomatic -- see the CONFIGURATION section above for why each
# threshold was chosen.

log_msg "========== Trimming reads with Trimmomatic =========="

# Locate the adapter FASTA once, outside the per-sample loop -- it's the
# same file for every sample.
# resolve the TruSeq3 adapter path once for all samples
ADAPTER_FILE=$(find_adapter_file)

for SRR in "${SRR_ACCESSIONS[@]}"; do                # trim each sample with Trimmomatic
    R1="${RAW_DIR}/${SRR}_1.fastq.gz"                # forward-mate raw FASTQ
    R2="${RAW_DIR}/${SRR}_2.fastq.gz"                # reverse-mate raw FASTQ

    # Trimmomatic PE produces 4 output files: paired + unpaired reads for
    # each of the two mates (a read can lose its partner during trimming if
    # the partner is discarded but this one survives).
    # forward reads whose mate also survived
    R1_PAIRED="${TRIMMED_TRIMMOMATIC_DIR}/${SRR}_1_paired.fastq.gz"
    # forward reads whose mate was dropped
    R1_UNPAIRED="${TRIMMED_TRIMMOMATIC_DIR}/${SRR}_1_unpaired.fastq.gz"
    # reverse reads whose mate also survived
    R2_PAIRED="${TRIMMED_TRIMMOMATIC_DIR}/${SRR}_2_paired.fastq.gz"
    # reverse reads whose mate was dropped
    R2_UNPAIRED="${TRIMMED_TRIMMOMATIC_DIR}/${SRR}_2_unpaired.fastq.gz"

    # skip missing/corrupt input
    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ] || ! is_valid_gz "${R1}" || ! is_valid_gz "${R2}"; then
        log_msg "WARNING: Paired FASTQs not found or corrupted for ${SRR}, skipping Trimmomatic."
        continue
    fi

    # Resume support: valid paired output already present means this
    # sample's Trimmomatic run already succeeded.
    # both paired outputs present and intact?
    if [ -s "${R1_PAIRED}" ] && [ -s "${R2_PAIRED}" ] && \
       is_valid_gz "${R1_PAIRED}" && is_valid_gz "${R2_PAIRED}"; then
        log_msg "Trimmomatic output for ${SRR} already exists, skipping."
        continue
    fi

    log_msg "Running Trimmomatic on ${SRR}..."

    # Build the trimming step list: adapter clipping only if we actually
    # found an adapter file, then the four quality/length filters from the
    # CONFIGURATION section, in Trimmomatic's required order.
    # accumulate Trimmomatic's ordered trimming steps here
    TRIM_STEPS=""
    # only add adapter clipping if the adapter file exists
    if [ -n "${ADAPTER_FILE}" ] && [ -f "${ADAPTER_FILE}" ]; then
        # 2:30:10:2:True -- Trimmomatic's own recommended ILLUMINACLIP
        # parameters for paired-end TruSeq3 data (seed mismatches:
        # palindrome/simple clip score thresholds: min adapter length:
        # keep both reads of a palindrome match).
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

    # -Xmx48g: the bioconda wrapper's default 1GB Java heap crashes on
    # datasets this size; -phred33 matches modern Illumina quality encoding;
    # -summary writes a per-sample trimming-stats file MultiQC can parse.
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

    log_msg "Trimmomatic completed for ${SRR}"
done

##==========================================================================##
##               1e: TRIMMING -- FASTP                                       ##
##==========================================================================##
# Independent second trimming track on the SAME raw reads, tuned to enforce
# an equivalent quality bar to Trimmomatic above (not fastp's own more
# lenient defaults), so the two tools' outputs are a fair comparison rather
# than one being trimmed harder than the other by accident.

log_msg "========== Trimming reads with fastp =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do                # trim each sample independently with fastp
    R1="${RAW_DIR}/${SRR}_1.fastq.gz"                # forward-mate raw FASTQ
    R2="${RAW_DIR}/${SRR}_2.fastq.gz"                # reverse-mate raw FASTQ

    R1_OUT="${TRIMMED_FASTP_DIR}/${SRR}_1.trimmed.fastq.gz"  # fastp's trimmed forward mate
    R2_OUT="${TRIMMED_FASTP_DIR}/${SRR}_2.trimmed.fastq.gz"  # fastp's trimmed reverse mate
    # machine-readable report -- what the fastp MultiQC report is built from
    JSON_OUT="${TRIMMED_FASTP_DIR}/${SRR}_fastp.json"
    HTML_OUT="${TRIMMED_FASTP_DIR}/${SRR}_fastp.html"   # human-readable per-sample report

    # skip missing/corrupt input
    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ] || ! is_valid_gz "${R1}" || ! is_valid_gz "${R2}"; then
        log_msg "WARNING: Paired FASTQs not found or corrupted for ${SRR}, skipping fastp."
        continue
    fi

    # Resume support: same pattern as the Trimmomatic loop above.
    # both outputs present and intact?
    if [ -s "${R1_OUT}" ] && [ -s "${R2_OUT}" ] && is_valid_gz "${R1_OUT}" && is_valid_gz "${R2_OUT}"; then
        log_msg "fastp output for ${SRR} already exists, skipping."
        continue
    fi

    log_msg "Running fastp on ${SRR}..."
    # trim + dedup with the shared Q15/len50 bar; emit JSON/HTML reports
    fastp \
        -i "${R1}" -I "${R2}" \
        -o "${R1_OUT}" -O "${R2_OUT}" \
        --cut_right --cut_right_window_size 4 --cut_right_mean_quality "${TRIM_AVGQUAL}" \
        --average_qual "${TRIM_AVGQUAL}" \
        --length_required "${TRIM_MINLEN}" \
        --dedup \
        --json "${JSON_OUT}" --html "${HTML_OUT}" \
        --thread "${THREADS}" \
        2>&1 | tee "${LOG_DIR}/${SRR}_fastp.log"
    check_exit "fastp failed for ${SRR}"             # abort on fastp failure
    # Notes on the flags above:
    #   --cut_right ... : fastp's equivalent of Trimmomatic's SLIDINGWINDOW
    #     -- slides a 4bp window from 5' to 3' and trims from the first
    #     window whose mean quality drops below Q15 onward. Off by default
    #     in fastp, so must be explicitly enabled for parity.
    #   --average_qual  : same whole-read Q15 average-quality gate as
    #     Trimmomatic's AVGQUAL.
    #   --length_required : same MINLEN:50 floor as Trimmomatic.
    #   --dedup : PCR/optical duplicate removal, off by default in fastp.
    #     Trimmomatic has no native dedup capability, so its own track gets
    #     a separate fastp dedup-only pass later in this script instead.
    #   Adapter trimming is fastp's default behavior (auto-detected from
    #   read overlap for paired-end data) -- no flag needed to enable it,
    #   unlike Trimmomatic which needs an explicit ILLUMINACLIP + adapter
    #   FASTA.

    log_msg "fastp completed for ${SRR}"
done

##==========================================================================##
##               1f: POST-TRIMMING QC WITH FASTQC (on Trimmomatic output)    ##
##==========================================================================##
# Re-run FastQC on the trimmed reads so pre-vs-post comparison in MultiQC
# shows the effect of trimming. Scoped to Trimmomatic's output specifically
# (fastp's own JSON/HTML report already serves this purpose for its track,
# aggregated separately below).

log_msg "========== Running post-trimming FastQC =========="

for SRR in "${SRR_ACCESSIONS[@]}"; do                # QC each sample's Trimmomatic output
    R1_PAIRED="${TRIMMED_TRIMMOMATIC_DIR}/${SRR}_1_paired.fastq.gz"  # trimmed forward paired reads
    R2_PAIRED="${TRIMMED_TRIMMOMATIC_DIR}/${SRR}_2_paired.fastq.gz"  # trimmed reverse paired reads

    # skip missing/corrupt input
    if [ ! -s "${R1_PAIRED}" ] || [ ! -s "${R2_PAIRED}" ] || ! is_valid_gz "${R1_PAIRED}" || ! is_valid_gz "${R2_PAIRED}"; then
        log_msg "WARNING: Trimmomatic output not found or corrupted for ${SRR}, skipping post-QC."
        continue
    fi

    # both post-trim reports already present?
    if [ -s "${FASTQC_POST_DIR}/${SRR}_1_paired_fastqc.html" ] && [ -s "${FASTQC_POST_DIR}/${SRR}_2_paired_fastqc.html" ]; then
        log_msg "Post-trim FastQC for ${SRR} already exists, skipping."
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
##               1g: MULTIQC -- TWO REPORTS                                  ##
##==========================================================================##
# Two separate reports rather than one combined report: FastQC + Trimmomatic
# summaries logically belong together (one trimming track's full QC story),
# while fastp's report is a self-contained, independent second opinion on
# the same raw data -- keeping them apart avoids one giant report where it's
# unclear which trimmer a given sample's numbers came from.

log_msg "========== Aggregating FastQC + Trimmomatic QC with MultiQC =========="

# aggregate pre/post FastQC + Trimmomatic logs into one report
multiqc \
    "${FASTQC_PRE_DIR}" "${FASTQC_POST_DIR}" "${LOG_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --filename "fastqc_trimmomatic_report" \
    --title "PRJNA207834 - Illumina FastQC + Trimmomatic Summary (HIV-1 Uganda A1/D)" \
    --force \
    2>&1 | tee "${LOG_DIR}/multiqc_fastqc_trimmomatic.log"
check_exit "MultiQC (FastQC + Trimmomatic) failed"   # abort if aggregation failed

log_msg "MultiQC report generated: ${MULTIQC_DIR}/fastqc_trimmomatic_report.html"

log_msg "========== Aggregating fastp reports with MultiQC =========="

# aggregate fastp's own JSON reports into a separate report
multiqc \
    "${TRIMMED_FASTP_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --filename "fastp_report" \
    --title "PRJNA207834 - Illumina fastp Summary (HIV-1 Uganda A1/D)" \
    --force \
    2>&1 | tee "${LOG_DIR}/multiqc_fastp.log"
check_exit "MultiQC (fastp) failed"                  # abort if aggregation failed

log_msg "MultiQC report generated: ${MULTIQC_DIR}/fastp_report.html"

##==========================================================================##
##      1h: DEDUPLICATION + KRAKEN2 HOST/BACTERIAL CONTAMINATION FILTER      ##
##==========================================================================##
# Per the tools review's recommended chain (fastp -> Kraken2 -> SHIVER):
# fastp's own --dedup above already deduplicates its track; Trimmomatic has
# no native dedup capability, so its track gets a separate fastp
# dedup-only pass here first. Both tracks then go through Kraken2 against a
# size-capped Standard database (bacteria+archaea+viral+human+UniVec_Core,
# ~16GB) -- a purely-viral database can't do this, since it has no
# human/bacterial genomes to match host contamination against, the
# dominant contamination source in HIV proviral sequencing.
#
# scripts/utils/kraken2_filter_reads.sh keeps unclassified reads and reads
# classified as viral, discarding anything whose assigned taxID descends
# from Homo sapiens (9606) or Bacteria (2), resolved via the database's own
# bundled nodes.dmp.

log_msg "========== STEP 1h: Deduplication + Kraken2 contamination filtering =========="

if [ ! -s "${KRAKEN2_DB}/nodes.dmp" ]; then          # if the Kraken2 DB isn't present...
    # ...warn and skip this whole block
    log_msg "WARNING: Kraken2 database not found at ${KRAKEN2_DB} -- skipping dedup+Kraken2 filtering. Download it first (see data/reference/kraken2_standard_16gb_db/SOURCE.md)."
else
    for SRR in "${SRR_ACCESSIONS[@]}"; do            # process both trimmer tracks per sample
        # --- Trimmomatic track: dedup first (fastp --dedup-only pass), then Kraken2 ---
        # Trimmomatic forward paired reads
        TRIMMOMATIC_R1="${TRIMMED_TRIMMOMATIC_DIR}/${SRR}_1_paired.fastq.gz"
        # Trimmomatic reverse paired reads
        TRIMMOMATIC_R2="${TRIMMED_TRIMMOMATIC_DIR}/${SRR}_2_paired.fastq.gz"
        # deduped forward mate for this track
        DEDUP_R1="${DEDUP_TRIMMOMATIC_DIR}/${SRR}_1.dedup.fastq.gz"
        # deduped reverse mate for this track
        DEDUP_R2="${DEDUP_TRIMMOMATIC_DIR}/${SRR}_2.dedup.fastq.gz"

        # only proceed with intact Trimmomatic output
        if [ -s "${TRIMMOMATIC_R1}" ] && [ -s "${TRIMMOMATIC_R2}" ] && is_valid_gz "${TRIMMOMATIC_R1}" && is_valid_gz "${TRIMMOMATIC_R2}"; then
            # skip dedup if a valid deduped pair already exists
            if [ ! -s "${DEDUP_R1}" ] || [ ! -s "${DEDUP_R2}" ] || ! is_valid_gz "${DEDUP_R1}" || ! is_valid_gz "${DEDUP_R2}"; then
                log_msg "Deduplicating Trimmomatic output for ${SRR}..."
                # run the fastp dedup-only wrapper on Trimmomatic's output
                bash "${BASE_DIR}/scripts/utils/fastp_dedup.sh" \
                    "${TRIMMOMATIC_R1}" "${TRIMMOMATIC_R2}" "${DEDUP_R1}" "${DEDUP_R2}" \
                    "${DEDUP_TRIMMOMATIC_DIR}/${SRR}_dedup_fastp" \
                    2>&1 | tee "${LOG_DIR}/${SRR}_dedup_trimmomatic.log"
                # abort on dedup failure
                check_exit "fastp dedup (Trimmomatic track) failed for ${SRR}"
            fi

            # only run Kraken2 if the deduped pair exists
            if [ -s "${DEDUP_R1}" ] && [ -s "${DEDUP_R2}" ]; then
                # host-removed forward mate
                KRAKEN_R1_OUT="${KRAKEN2_TRIMMOMATIC_DIR}/${SRR}_1.kraken_filtered.fastq.gz"
                # host-removed reverse mate
                KRAKEN_R2_OUT="${KRAKEN2_TRIMMOMATIC_DIR}/${SRR}_2.kraken_filtered.fastq.gz"
                # skip if already filtered
                if [ ! -s "${KRAKEN_R1_OUT}" ] || [ ! -s "${KRAKEN_R2_OUT}" ]; then
                    log_msg "Running Kraken2 (Trimmomatic track) on ${SRR}..."
                    # drop host/bacterial reads from the deduped Trimmomatic pair
                    bash "${BASE_DIR}/scripts/utils/kraken2_filter_reads.sh" \
                        "${DEDUP_R1}" "${DEDUP_R2}" "${KRAKEN2_TRIMMOMATIC_DIR}" "${SRR}" "${KRAKEN2_DB}" \
                        2>&1 | tee "${LOG_DIR}/${SRR}_kraken2_trimmomatic.log"
                    # abort on Kraken2 failure
                    check_exit "Kraken2 (Trimmomatic track) failed for ${SRR}"
                fi
            fi
        else
            # nothing to filter for this track
            log_msg "WARNING: Trimmomatic output not found for ${SRR}, skipping its dedup+Kraken2 pass."
        fi

        # --- fastp track: already deduplicated above, straight to Kraken2 ---
        # fastp trimmed+deduped forward mate
        FASTP_R1="${TRIMMED_FASTP_DIR}/${SRR}_1.trimmed.fastq.gz"
        # fastp trimmed+deduped reverse mate
        FASTP_R2="${TRIMMED_FASTP_DIR}/${SRR}_2.trimmed.fastq.gz"
        # only proceed with intact fastp output
        if [ -s "${FASTP_R1}" ] && [ -s "${FASTP_R2}" ] && is_valid_gz "${FASTP_R1}" && is_valid_gz "${FASTP_R2}"; then
            # host-removed forward mate (fastp track)
            KRAKEN_R1_OUT="${KRAKEN2_FASTP_DIR}/${SRR}_1.kraken_filtered.fastq.gz"
            # host-removed reverse mate (fastp track)
            KRAKEN_R2_OUT="${KRAKEN2_FASTP_DIR}/${SRR}_2.kraken_filtered.fastq.gz"
            # skip if already filtered
            if [ ! -s "${KRAKEN_R1_OUT}" ] || [ ! -s "${KRAKEN_R2_OUT}" ]; then
                log_msg "Running Kraken2 (fastp track) on ${SRR}..."
                # drop host/bacterial reads from the fastp pair (production input)
                bash "${BASE_DIR}/scripts/utils/kraken2_filter_reads.sh" \
                    "${FASTP_R1}" "${FASTP_R2}" "${KRAKEN2_FASTP_DIR}" "${SRR}" "${KRAKEN2_DB}" \
                    2>&1 | tee "${LOG_DIR}/${SRR}_kraken2_fastp.log"
                check_exit "Kraken2 (fastp track) failed for ${SRR}"  # abort on Kraken2 failure
            fi
        else
            # nothing to filter for this track
            log_msg "WARNING: fastp output not found for ${SRR}, skipping its Kraken2 pass."
        fi

        log_msg "Deduplication + Kraken2 filtering completed for ${SRR}"
    done
fi

##==========================================================================##
##                          STEP 1 COMPLETE                                  ##
##==========================================================================##

log_msg "=========================================="
log_msg "  STEP 1 COMPLETE (download, QC, trim)     "
log_msg "=========================================="
log_msg "BioProject:                 ${BIOPROJECT}"
log_msg "Samples processed:         ${#SRR_ACCESSIONS[@]}"
log_msg "Raw data:                  ${RAW_DIR}"
log_msg "Trimmomatic output:        ${TRIMMED_TRIMMOMATIC_DIR}"
log_msg "fastp output:              ${TRIMMED_FASTP_DIR}"
log_msg "Deduplicated (Trimmomatic): ${DEDUP_TRIMMOMATIC_DIR}"
log_msg "Kraken2-filtered (Trimmomatic): ${KRAKEN2_TRIMMOMATIC_DIR}"
log_msg "Kraken2-filtered (fastp, recommended production input): ${KRAKEN2_FASTP_DIR}"
log_msg "FastQC+Trimmomatic MultiQC: ${MULTIQC_DIR}/fastqc_trimmomatic_report.html"
log_msg "fastp MultiQC:              ${MULTIQC_DIR}/fastp_report.html"
log_msg "=========================================="
