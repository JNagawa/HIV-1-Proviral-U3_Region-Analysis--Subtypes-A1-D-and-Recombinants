#!/bin/bash
# This script is for the qc analysis of single_end SMRTcap HiFi/CCS fastq sequences
# It involves four key steps and two or more tools are used for each step
# for comparison to guide choice of which tool is best for this kind of data.
#   1. NanoPlot QC (long-read analogue of FastQC; review's PacBio QC pick)
#   2. length/quality FILTER, comparing NanoFilt vs fastp vs chopper
#   3. Kraken2 host/bacterial removal (single-end) after each filter
#   4. DEDUP, comparing fastp --dedup vs seqkit rmdup
# For step 3, all the three filters are held to the same HiFi bar (LEN_MIN/Q_MIN) so the
# comparison reflects the tools, not mismatched settings
# Usage: ./download_qc_pacbio.sh   (run via sbatch scripts/utils/run_comparison_step.slurm.sh)
set -uo pipefail                                     # -u errors on unset vars, pipefail fails a pipe if any stage fails; no -e so one tool failing doesn't kill the whole comparison

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"           # absolute path of this script's own dir, so paths work no matter where it's launched from
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"     # repo top-level (three dirs up), the base for every other path below
source "${REPO_ROOT}/scripts/common/lib_compare.sh"  # load shared helpers: measure_and_run, parse_time_metrics, append_summary_row, subset_accessions

RAW_DIR="${REPO_ROOT}/data/raw/pacbio"                            # input: downloaded HiFi FASTQs live here
RESULTS_DIR="${REPO_ROOT}/results/download_qc/pacbio"             # output: all QC results + summary go here
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"                          # the single TSV every tool appends a timing/validity row to
KRAKEN2_DB="${REPO_ROOT}/data/reference/kraken2_standard_16gb_db" # Kraken2 database used for host/bacterial read removal
mkdir -p "${RESULTS_DIR}"                                                                        # create the results dir (and parents) if it doesn't exist
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"  # seed the manual notes file from the template only on the first run

# HiFi reads are long (~10-15kb) and already high-accuracy (Q20+ by
# definition of "HiFi"), so filtering is light: keep reads >=1000bp (enough
# to carry a chunk of the ~9kb provirus) and mean quality >=Q20. All three
# filters use exactly these two thresholds.
LEN_MIN=1000                                         # minimum read length (bp) kept by every filter
Q_MIN=20                                             # minimum mean read quality (Phred) kept by every filter
THREADS="${THREADS:-4}"                              # thread count; honour an externally-set THREADS, otherwise default to 4
export THREADS                                       # export so child scripts/tools (e.g. the Kraken2 wrapper) inherit it

reads_in_fastq() {  # count reads in a (gzipped) fastq: lines/4
    local f="$1"                                     # arg 1 = path to the fastq(.gz) to count
    [ -s "${f}" ] || { echo 0; return; }             # if the file is missing/empty, report 0 and bail (avoids a zcat error)
    echo $(( $(zcat "${f}" 2>/dev/null | wc -l) / 4 ))  # 4 lines per fastq record, so total lines / 4 = read count
}

for SRR in $(subset_accessions pacbio "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do  # loop over just the PacBio accessions chosen for this comparison
    READS="${RAW_DIR}/${SRR}.fastq.gz"               # expected raw HiFi FASTQ for this accession
    if [ ! -s "${READS}" ]; then                     # if that FASTQ is missing or empty...
        echo "WARNING: HiFi FASTQ missing for ${SRR} (run download_pacbio.slurm.sh first), skipping." >&2  # ...warn on stderr...
        continue                                     # ...and skip to the next accession instead of crashing
    fi

    # --- NanoPlot pre-filter QC: one baseline per sample. ---
    NP_PRE="${RESULTS_DIR}/nanoplot_pre_out/${SRR}"  # per-sample output dir for the pre-filter QC report
    if [ ! -s "${NP_PRE}/NanoPlot-report.html" ]; then  # skip if the report already exists (makes reruns idempotent)
        echo "=== NanoPlot (pre-filter) on ${SRR} ==="  # progress marker in the log
        mkdir -p "${NP_PRE}"                         # ensure the per-sample output dir exists
        NanoPlot --fastq "${READS}" --outdir "${NP_PRE}" --prefix "${SRR}_pre_" \
            --threads "${THREADS}" --tsv_stats > "${RESULTS_DIR}/nanoplot_pre_${SRR}.log" 2>&1  # QC the raw reads; --tsv_stats gives machine-readable stats; redirect stdout+stderr to the log
    fi

    for TOOL in nanofilt fastp chopper; do           # run each filter tool on the same input for a head-to-head comparison
        if [ "${TOOL}" = "chopper" ] && ! command -v chopper >/dev/null 2>&1; then  # chopper is optional...
            echo "NOTE: chopper not installed, skipping (see HIV_U3analysis_env.yml)." >&2  # ...note its absence...
            continue                                 # ...and skip it if the binary isn't on PATH
        fi
        OUTDIR="${RESULTS_DIR}/${TOOL}_out"          # per-tool output dir (e.g. fastp_out/)
        mkdir -p "${OUTDIR}"                         # create it if needed
        FILT="${OUTDIR}/${SRR}.filtered.fastq.gz"    # this tool's filtered-reads output for this sample
        TIMELOG="${RESULTS_DIR}/${TOOL}_${SRR}.time" # file where measure_and_run records wallclock/RSS
        LOG="${RESULTS_DIR}/${TOOL}_${SRR}.log"      # captured stdout+stderr of the tool

        echo "=== ${TOOL} filter on ${SRR} ==="      # progress marker in the log
        case "${TOOL}" in                            # dispatch to the right command per tool (same thresholds, tool-specific syntax)
            nanofilt)
                measure_and_run "${TIMELOG}" -- \
                    bash -c 'gunzip -c "$1" | NanoFilt -q "$2" -l "$3" | gzip > "$4"' _ \
                    "${READS}" "${Q_MIN}" "${LEN_MIN}" "${FILT}" > "${LOG}" 2>&1 ;;  # NanoFilt reads stdin, so decompress in, filter, gzip out
            fastp)
                measure_and_run "${TIMELOG}" -- \
                    fastp -i "${READS}" -o "${FILT}" \
                        --disable_adapter_trimming \
                        --average_qual "${Q_MIN}" --length_required "${LEN_MIN}" \
                        --json "${OUTDIR}/${SRR}_fastp.json" --html "${OUTDIR}/${SRR}_fastp.html" \
                        --thread "${THREADS}" > "${LOG}" 2>&1 ;;  # fastp works on the .gz directly; adapter trimming off (HiFi is adapter-clean); emit JSON/HTML reports
            chopper)
                measure_and_run "${TIMELOG}" -- \
                    bash -c 'gunzip -c "$1" | chopper -q "$2" -l "$3" --threads "$4" | gzip > "$5"' _ \
                    "${READS}" "${Q_MIN}" "${LEN_MIN}" "${THREADS}" "${FILT}" > "${LOG}" 2>&1 ;;  # chopper is also stdin/stdout, so same decompress|filter|gzip pattern as NanoFilt
        esac
        EXIT_CODE=$?                                 # capture the tool's exit status before $? is overwritten
        parse_time_metrics "${TIMELOG}"              # set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file

        N_FILT=$(reads_in_fastq "${FILT}")           # how many reads survived this filter
        VALID=0; METRIC="n/a"                        # assume invalid until proven otherwise
        if [ -s "${FILT}" ] && gzip -t "${FILT}" 2>/dev/null && [ "${N_FILT}" -gt 0 ]; then  # valid = non-empty, not-corrupt gzip, and at least one read
            VALID=1; METRIC="${N_FILT} reads passed (>=${LEN_MIN}bp, >=Q${Q_MIN})"  # mark valid and record the human-readable key metric
        fi
        append_summary_row "download_qc_pacbio" "${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"  # write this filter's row to summary.tsv

        # --- Kraken2 host/bacterial removal on this filter's output (single-end). ---
        if [ ! -d "${KRAKEN2_DB}" ]; then            # if the Kraken2 DB isn't present...
            echo "NOTE: Kraken2 DB ${KRAKEN2_DB} absent, skipping host removal for ${TOOL}/${SRR}." >&2  # ...note it and skip host removal
        elif [ "${VALID}" = "1" ]; then              # only run Kraken2 when the filter actually produced usable reads
            KOUT="${RESULTS_DIR}/kraken2_${TOOL}_out"          # per-filter Kraken2 output dir
            KTIME="${RESULTS_DIR}/kraken2_${TOOL}_${SRR}.time" # timing file for this Kraken2 run
            KLOG="${RESULTS_DIR}/kraken2_${TOOL}_${SRR}.log"   # log for this Kraken2 run
            echo "=== Kraken2 host removal after ${TOOL} on ${SRR} ==="  # progress marker
            measure_and_run "${KTIME}" -- \
                "${REPO_ROOT}/scripts/utils/kraken2_filter_reads_se.sh" \
                "${FILT}" "${KOUT}" "${SRR}" "${KRAKEN2_DB}" > "${KLOG}" 2>&1  # single-end Kraken2 wrapper: classify reads and drop host/bacterial ones
            KEXIT=$?                                  # capture Kraken2's exit status
            parse_time_metrics "${KTIME}"            # refresh WALLCLOCK_SEC / PEAK_RSS_MB from the Kraken2 timing file
            KFILT="${KOUT}/${SRR}.kraken_filtered.fastq.gz"  # the host-removed reads the wrapper produces
            N_KRAK=$(reads_in_fastq "${KFILT}")      # reads remaining after host removal
            KVALID=0; KMETRIC="n/a"                  # default to invalid until checked
            [ -s "${KFILT}" ] && [ "${N_KRAK}" -gt 0 ] && { KVALID=1; KMETRIC="${N_KRAK} reads retained after host removal"; }  # valid if the output exists and has reads
            append_summary_row "download_qc_pacbio" "kraken2_after_${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${KEXIT}" "${KVALID}" "${KMETRIC}"  # record the Kraken2 row
        fi
    done

    # --- Dedup comparison: fastp --dedup vs seqkit rmdup, on the fastp-filtered
    # + Kraken2-cleaned reads (the canonical clean input for assembly). ---
    CLEAN="${RESULTS_DIR}/kraken2_fastp_out/${SRR}.kraken_filtered.fastq.gz"  # preferred dedup input: fastp-filtered + host-removed reads
    [ -s "${CLEAN}" ] || CLEAN="${RESULTS_DIR}/fastp_out/${SRR}.filtered.fastq.gz"  # fall back to just fastp-filtered reads if Kraken2 was skipped
    if [ -s "${CLEAN}" ]; then                       # only dedup if we actually have a clean input
        N_BEFORE=$(reads_in_fastq "${CLEAN}")        # read count before dedup, for the "N duplicates removed" metric
        for DTOOL in fastp seqkit; do                # compare the two dedup tools on identical input
            DOUT="${RESULTS_DIR}/dedup_${DTOOL}_out"; mkdir -p "${DOUT}"  # per-tool dedup output dir (created if needed)
            DFILE="${DOUT}/${SRR}.dedup.fastq.gz"    # deduplicated reads output
            DTIME="${RESULTS_DIR}/dedup_${DTOOL}_${SRR}.time"  # timing file for this dedup run
            DLOG="${RESULTS_DIR}/dedup_${DTOOL}_${SRR}.log"    # log for this dedup run
            echo "=== dedup (${DTOOL}) on ${SRR} ==="  # progress marker
            case "${DTOOL}" in                        # tool-specific dedup command
                fastp)
                    measure_and_run "${DTIME}" -- \
                        fastp -i "${CLEAN}" -o "${DFILE}" --dedup \
                            --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering \
                            --json "${DOUT}/${SRR}_dedup.json" --html "${DOUT}/${SRR}_dedup.html" \
                            --thread "${THREADS}" > "${DLOG}" 2>&1 ;;  # fastp in dedup-only mode: all other filtering disabled so it only removes duplicates
                seqkit)
                    measure_and_run "${DTIME}" -- \
                        bash -c 'seqkit rmdup -s "$1" -o "$2"' _ "${CLEAN}" "${DFILE}" > "${DLOG}" 2>&1 ;;  # seqkit rmdup -s dedups by sequence content
            esac
            DEXIT=$?                                  # capture the dedup tool's exit status
            parse_time_metrics "${DTIME}"            # refresh timing/RSS from this dedup run
            N_AFTER=$(reads_in_fastq "${DFILE}")     # read count after dedup
            DVALID=0; DMETRIC="n/a"                  # default to invalid until checked
            [ -s "${DFILE}" ] && [ "${N_AFTER}" -gt 0 ] && { DVALID=1; DMETRIC="${N_AFTER}/${N_BEFORE} reads kept ($((N_BEFORE-N_AFTER)) duplicates removed)"; }  # valid if output has reads; report kept vs removed
            append_summary_row "download_qc_pacbio" "dedup_${DTOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${DEXIT}" "${DVALID}" "${DMETRIC}"  # record the dedup row
        done
    fi
done

echo "Done. See ${SUMMARY_TSV}"                      # final confirmation pointing the user at the results table
