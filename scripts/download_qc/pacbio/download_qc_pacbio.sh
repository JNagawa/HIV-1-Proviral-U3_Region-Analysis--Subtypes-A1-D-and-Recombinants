#!/bin/bash
#SBATCH --job-name=pb_download_qc
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# 64G because Kraken2 loads the whole 16GB standard DB into RAM and is run
# three times per sample (once after each filter); 24h because this is five
# raw HiFi runs (~12GB gzipped) through NanoPlot + 3 filters + 3 Kraken2
# passes + 2 dedup tools, and Kraken2 on 10-15kb reads is the long pole --
# merely downloading this same data took 4h40m (job 111580).
# This script is for the qc analysis of single_end SMRTcap HiFi/CCS fastq sequences
# It involves four key steps and two or more tools are used for each step
# for comparison to guide choice of which tool is best for this kind of data.
#   1. NanoPlot QC (long-read analogue of FastQC; review's PacBio QC pick)
#   2. length/quality FILTER, comparing NanoFilt vs fastp vs chopper
#   3. Kraken2 host/bacterial removal (single-end) after each filter
#   4. DEDUP, comparing fastp --dedup vs seqkit rmdup
# For step 3, all the three filters are held to the same HiFi bar (LEN_MIN/Q_MIN) so the
# comparison reflects the tools, not mismatched settings
# Usage: sbatch scripts/download_qc/pacbio/download_qc_pacbio.sh
#   Submit from the repo root -- the --output/--error paths above are relative to
#   the submitting directory, so logs/ must exist where you run sbatch.
#   Still works unchanged as `./download_qc_pacbio.sh` for a quick direct run, and
#   via `sbatch scripts/utils/run_comparison_step.slurm.sh <this script>`; the
#   #SBATCH lines are ordinary comments in both of those cases.
# -u errors on unset vars, pipefail fails a pipe if any stage fails; no -e so one tool failing
# doesn't kill the whole comparison
set -uo pipefail

# Activate the tool env ourselves so this is submittable on its own, not only
# through run_comparison_step.slurm.sh (which does the same thing). Skipped when
# the env is already active, so a direct run in an activated shell is untouched.
if [ "${CONDA_DEFAULT_ENV:-}" != "HIV_U3analysis" ]; then
    CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"  # usual miniconda hook location
    # source it, or fall back to whatever conda is on PATH
    if [ -f "$CONDA_SH" ]; then source "$CONDA_SH"; else source "$(conda info --base)/etc/profile.d/conda.sh"; fi
    conda activate HIV_U3analysis                    # env holding NanoPlot/fastp/chopper/Kraken2/seqkit
fi

# absolute path of this script's own dir, paths work no matter where it's launched from
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
# repo top-level (three dirs up), the base for every other path below
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
# load shared helpers: measure_and_run, parse_time_metrics, append_summary_row, subset_accessions
source "${REPO_ROOT}/scripts/common/lib_compare.sh"

RAW_DIR="${REPO_ROOT}/data/raw/pacbio"                            # input: downloaded HiFi FASTQs
# output: all QC results + summary go here
RESULTS_DIR="${REPO_ROOT}/results/download_qc/pacbio"
# the single TSV every tool appends a timing/validity row to
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
# Kraken2 database used for host/bacterial read removal
KRAKEN2_DB="${REPO_ROOT}/data/reference/kraken2_standard_16gb_db"
# create the results dir (and parents) if it doesn't exist
mkdir -p "${RESULTS_DIR}"
# seed the manual notes file from the template only on the first run
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] \
    || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" \
          "${RESULTS_DIR}/ease_of_use_notes.md"

# HiFi reads are long (~10-15kb) and already high-accuracy (Q20+ by
# definition of "HiFi"), so filtering is light: keep reads >=1000bp (enough
# to carry a chunk of the ~9kb provirus) and mean quality >=Q20. All three
# filters use exactly these two thresholds.
LEN_MIN=1000                                         # minimum read length (bp) kept by every filter
# minimum mean read quality (Phred) kept by every filter
Q_MIN=20
# thread count; honour an externally-set THREADS, else Slurm's allocation, else 4
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
# export so child scripts/tools (e.g. the Kraken2 wrapper) inherit it
export THREADS

reads_in_fastq() {  # count reads in a (gzipped) fastq: lines/4
    local f="$1"                                     # arg 1 = path to the fastq(.gz) to count
    # if the file is missing/empty, report 0 and bail (avoids a zcat error)
    [ -s "${f}" ] || { echo 0; return; }
    # 4 lines per fastq record, so total lines / 4 = read count
    echo $(( $(zcat "${f}" 2>/dev/null | wc -l) / 4 ))
}

# Loop over the `pacbio_sra` rows, NOT `pacbio`: the four `pacbio` samples are the
# host-N-masked reads that enter the harness at proviral extraction already QC'd, so
# they have no raw FASTQ here. The QC-tool comparison needs the raw SRA HiFi runs.
for SRR in $(subset_accessions pacbio_sra "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do
    READS="${RAW_DIR}/${SRR}.fastq.gz"               # expected raw HiFi FASTQ for this accession
    if [ ! -s "${READS}" ]; then                     # if that FASTQ is missing or empty...
        # ...warn on stderr...
        echo "WARNING: HiFi FASTQ missing for ${SRR} (run download_pacbio.slurm.sh first), skipping." >&2
        # ...and skip to the next accession instead of crashing
        continue
    fi

    # --- NanoPlot pre-filter QC: one baseline per sample. ---
    # per-sample output dir for the pre-filter QC report
    NP_PRE="${RESULTS_DIR}/nanoplot_pre_out/${SRR}"
    # skip if the report already exists (makes reruns idempotent)
    if [ ! -s "${NP_PRE}/NanoPlot-report.html" ]; then
        echo "=== NanoPlot (pre-filter) on ${SRR} ==="  # progress marker in the log
        mkdir -p "${NP_PRE}"                         # ensure the per-sample output dir exists
        # QC the raw reads; --tsv_stats gives machine-readable stats; redirect stdout+stderr to the
        # log
        NanoPlot --fastq "${READS}" --outdir "${NP_PRE}" --prefix "${SRR}_pre_" \
            --threads "${THREADS}" --tsv_stats > "${RESULTS_DIR}/nanoplot_pre_${SRR}.log" 2>&1
    fi

    # run each filter tool on the same input for a head-to-head comparison
    for TOOL in nanofilt fastp chopper; do
        # chopper is optional...
        if [ "${TOOL}" = "chopper" ] && ! command -v chopper >/dev/null 2>&1; then
            # ...note its absence...
            echo "NOTE: chopper not installed, skipping (see HIV_U3analysis_env.yml)." >&2
            continue                                 # ...and skip it if the binary isn't on PATH
        fi
        OUTDIR="${RESULTS_DIR}/${TOOL}_out"          # per-tool output dir (e.g. fastp_out/)
        mkdir -p "${OUTDIR}"                         # create it if needed
        # this tool's filtered-reads output for this sample
        FILT="${OUTDIR}/${SRR}.filtered.fastq.gz"
        # file where measure_and_run records wallclock/RSS
        TIMELOG="${RESULTS_DIR}/${TOOL}_${SRR}.time"
        LOG="${RESULTS_DIR}/${TOOL}_${SRR}.log"      # captured stdout+stderr of the tool

        echo "=== ${TOOL} filter on ${SRR} ==="      # progress marker in the log
        # dispatch to the right command per tool (same thresholds, tool-specific syntax)
        case "${TOOL}" in
            nanofilt)
                # NanoFilt reads stdin, so decompress in, filter, gzip out
                measure_and_run "${TIMELOG}" -- \
                    bash -c 'gunzip -c "$1" | NanoFilt -q "$2" -l "$3" | gzip > "$4"' _ \
                    "${READS}" "${Q_MIN}" "${LEN_MIN}" "${FILT}" > "${LOG}" 2>&1 ;;
            fastp)
                # fastp works on the .gz directly; adapter trimming off (HiFi is adapter-clean);
                # emit JSON/HTML reports
                measure_and_run "${TIMELOG}" -- \
                    fastp -i "${READS}" -o "${FILT}" \
                        --disable_adapter_trimming \
                        --average_qual "${Q_MIN}" --length_required "${LEN_MIN}" \
                        --json "${OUTDIR}/${SRR}_fastp.json" --html "${OUTDIR}/${SRR}_fastp.html" \
                        --thread "${THREADS}" > "${LOG}" 2>&1 ;;
            chopper)
                # chopper is also stdin/stdout, so same decompress|filter|gzip pattern as NanoFilt
                measure_and_run "${TIMELOG}" -- \
                    bash -c 'gunzip -c "$1" | chopper -q "$2" -l "$3" --threads "$4" | gzip > "$5"' _ \
                    "${READS}" "${Q_MIN}" "${LEN_MIN}" "${THREADS}" "${FILT}" > "${LOG}" 2>&1 ;;
        esac
        # capture the tool's exit status before $? is overwritten
        EXIT_CODE=$?
        # set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file
        parse_time_metrics "${TIMELOG}"

        N_FILT=$(reads_in_fastq "${FILT}")           # how many reads survived this filter
        VALID=0; METRIC="n/a"                        # assume invalid until proven otherwise
        # valid = non-empty, not-corrupt gzip, and at least one read
        if [ -s "${FILT}" ] && gzip -t "${FILT}" 2>/dev/null && [ "${N_FILT}" -gt 0 ]; then
            # mark valid and record the human-readable key metric
            VALID=1; METRIC="${N_FILT} reads passed (>=${LEN_MIN}bp, >=Q${Q_MIN})"
        fi
        # write this filter's row to summary.tsv
        append_summary_row "download_qc_pacbio" "${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"

        # --- Kraken2 host/bacterial removal on this filter's output (single-end). ---
        if [ ! -d "${KRAKEN2_DB}" ]; then            # if the Kraken2 DB isn't present...
            # ...note it and skip host removal
            echo "NOTE: Kraken2 DB ${KRAKEN2_DB} absent, skipping host removal for ${TOOL}/${SRR}." >&2
        # only run Kraken2 when the filter actually produced usable reads
        elif [ "${VALID}" = "1" ]; then
            KOUT="${RESULTS_DIR}/kraken2_${TOOL}_out"          # per-filter Kraken2 output dir
            KTIME="${RESULTS_DIR}/kraken2_${TOOL}_${SRR}.time" # timing file for this Kraken2 run
            KLOG="${RESULTS_DIR}/kraken2_${TOOL}_${SRR}.log"   # log for this Kraken2 run
            echo "=== Kraken2 host removal after ${TOOL} on ${SRR} ==="  # progress marker
            # single-end Kraken2 wrapper: classify reads and drop host/bacterial ones
            measure_and_run "${KTIME}" -- \
                "${REPO_ROOT}/scripts/utils/kraken2_filter_reads_se.sh" \
                "${FILT}" "${KOUT}" "${SRR}" "${KRAKEN2_DB}" > "${KLOG}" 2>&1
            KEXIT=$?                                  # capture Kraken2's exit status
            # refresh WALLCLOCK_SEC / PEAK_RSS_MB from the Kraken2 timing file
            parse_time_metrics "${KTIME}"
            # the host-removed reads the wrapper produces
            KFILT="${KOUT}/${SRR}.kraken_filtered.fastq.gz"
            N_KRAK=$(reads_in_fastq "${KFILT}")      # reads remaining after host removal
            KVALID=0; KMETRIC="n/a"                  # default to invalid until checked
            # valid if the output exists and has reads
            [ -s "${KFILT}" ] && [ "${N_KRAK}" -gt 0 ] && { KVALID=1; KMETRIC="${N_KRAK} reads retained after host removal"; }
            # record the Kraken2 row
            append_summary_row "download_qc_pacbio" "kraken2_after_${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${KEXIT}" "${KVALID}" "${KMETRIC}"
        fi
    done

    # --- Dedup comparison: fastp --dedup vs seqkit rmdup, on the fastp-filtered
    # + Kraken2-cleaned reads (the canonical clean input for assembly). ---
    # preferred dedup input: fastp-filtered + host-removed reads
    CLEAN="${RESULTS_DIR}/kraken2_fastp_out/${SRR}.kraken_filtered.fastq.gz"
    # fall back to just fastp-filtered reads if Kraken2 was skipped
    [ -s "${CLEAN}" ] || CLEAN="${RESULTS_DIR}/fastp_out/${SRR}.filtered.fastq.gz"
    if [ -s "${CLEAN}" ]; then                       # only dedup if we actually have a clean input
        # read count before dedup, for the "N duplicates removed" metric
        N_BEFORE=$(reads_in_fastq "${CLEAN}")
        # compare the two dedup tools on identical input
        for DTOOL in fastp seqkit; do
            # per-tool dedup output dir (created if needed)
            DOUT="${RESULTS_DIR}/dedup_${DTOOL}_out"; mkdir -p "${DOUT}"
            DFILE="${DOUT}/${SRR}.dedup.fastq.gz"    # deduplicated reads output
            DTIME="${RESULTS_DIR}/dedup_${DTOOL}_${SRR}.time"  # timing file for this dedup run
            DLOG="${RESULTS_DIR}/dedup_${DTOOL}_${SRR}.log"    # log for this dedup run
            echo "=== dedup (${DTOOL}) on ${SRR} ==="  # progress marker
            case "${DTOOL}" in                        # tool-specific dedup command
                fastp)
                    # fastp in dedup-only mode: all other filtering disabled so it only removes
                    # duplicates
                    measure_and_run "${DTIME}" -- \
                        fastp -i "${CLEAN}" -o "${DFILE}" --dedup \
                            --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering \
                            --json "${DOUT}/${SRR}_dedup.json" --html "${DOUT}/${SRR}_dedup.html" \
                            --thread "${THREADS}" > "${DLOG}" 2>&1 ;;
                seqkit)
                    # seqkit rmdup -s dedups by sequence content
                    measure_and_run "${DTIME}" -- \
                        bash -c 'seqkit rmdup -s "$1" -o "$2"' _ "${CLEAN}" "${DFILE}" > "${DLOG}" 2>&1 ;;
            esac
            DEXIT=$?                                  # capture the dedup tool's exit status
            parse_time_metrics "${DTIME}"            # refresh timing/RSS from this dedup run
            N_AFTER=$(reads_in_fastq "${DFILE}")     # read count after dedup
            DVALID=0; DMETRIC="n/a"                  # default to invalid until checked
            # valid if output has reads; report kept vs removed
            [ -s "${DFILE}" ] && [ "${N_AFTER}" -gt 0 ] && { DVALID=1; DMETRIC="${N_AFTER}/${N_BEFORE} reads kept ($((N_BEFORE-N_AFTER)) duplicates removed)"; }
            # record the dedup row
            append_summary_row "download_qc_pacbio" "dedup_${DTOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${DEXIT}" "${DVALID}" "${DMETRIC}"
        done
    fi
done

# final confirmation pointing the user at the results table
echo "Done. See ${SUMMARY_TSV}"
