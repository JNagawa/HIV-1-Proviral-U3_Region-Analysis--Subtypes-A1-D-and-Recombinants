#!/bin/bash
# Download-QC step of the PacBio (HIV-SMRTcap) tool-comparison harness. The
# SMRTcap data is already HiFi/CCS (the SMRT Link `ccs` step is done before
# SRA deposition), so this starts from single-end HiFi FASTQs and does:
#   1. NanoPlot QC (long-read analogue of FastQC; review's PacBio QC pick)
#   2. length/quality FILTER, comparing NanoFilt vs fastp vs chopper
#   3. Kraken2 host/bacterial removal (single-end) after each filter
#   4. DEDUP, comparing fastp --dedup vs seqkit rmdup
# All three filters are held to the same HiFi bar (LEN_MIN/Q_MIN) so the
# comparison reflects the tools, not mismatched settings -- same principle as
# the Illumina download_qc harness.
# Usage: ./download_qc_pacbio.sh   (run via sbatch scripts/utils/run_comparison_step.slurm.sh)
set -uo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
source "${REPO_ROOT}/scripts/common/lib_compare.sh"

RAW_DIR="${REPO_ROOT}/data/raw/pacbio"
RESULTS_DIR="${REPO_ROOT}/results/download_qc/pacbio"
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
KRAKEN2_DB="${REPO_ROOT}/data/reference/kraken2_standard_16gb_db"
mkdir -p "${RESULTS_DIR}"
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

# HiFi reads are long (~10-15kb) and already high-accuracy (Q20+ by
# definition of "HiFi"), so filtering is light: keep reads >=1000bp (enough
# to carry a chunk of the ~9kb provirus) and mean quality >=Q20. All three
# filters use exactly these two thresholds.
LEN_MIN=1000
Q_MIN=20
THREADS="${THREADS:-4}"
export THREADS

reads_in_fastq() {  # count reads in a (gzipped) fastq: lines/4
    local f="$1"
    [ -s "${f}" ] || { echo 0; return; }
    echo $(( $(zcat "${f}" 2>/dev/null | wc -l) / 4 ))
}

for SRR in $(subset_accessions pacbio "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do
    READS="${RAW_DIR}/${SRR}.fastq.gz"
    if [ ! -s "${READS}" ]; then
        echo "WARNING: HiFi FASTQ missing for ${SRR} (run download_pacbio.slurm.sh first), skipping." >&2
        continue
    fi

    # --- NanoPlot pre-filter QC: one baseline per sample. ---
    NP_PRE="${RESULTS_DIR}/nanoplot_pre_out/${SRR}"
    if [ ! -s "${NP_PRE}/NanoPlot-report.html" ]; then
        echo "=== NanoPlot (pre-filter) on ${SRR} ==="
        mkdir -p "${NP_PRE}"
        NanoPlot --fastq "${READS}" --outdir "${NP_PRE}" --prefix "${SRR}_pre_" \
            --threads "${THREADS}" --tsv_stats > "${RESULTS_DIR}/nanoplot_pre_${SRR}.log" 2>&1
    fi

    for TOOL in nanofilt fastp chopper; do
        if [ "${TOOL}" = "chopper" ] && ! command -v chopper >/dev/null 2>&1; then
            echo "NOTE: chopper not installed, skipping (see HIV_U3analysis_env.yml)." >&2
            continue
        fi
        OUTDIR="${RESULTS_DIR}/${TOOL}_out"
        mkdir -p "${OUTDIR}"
        FILT="${OUTDIR}/${SRR}.filtered.fastq.gz"
        TIMELOG="${RESULTS_DIR}/${TOOL}_${SRR}.time"
        LOG="${RESULTS_DIR}/${TOOL}_${SRR}.log"

        echo "=== ${TOOL} filter on ${SRR} ==="
        case "${TOOL}" in
            nanofilt)
                measure_and_run "${TIMELOG}" -- \
                    bash -c 'gunzip -c "$1" | NanoFilt -q "$2" -l "$3" | gzip > "$4"' _ \
                    "${READS}" "${Q_MIN}" "${LEN_MIN}" "${FILT}" > "${LOG}" 2>&1 ;;
            fastp)
                measure_and_run "${TIMELOG}" -- \
                    fastp -i "${READS}" -o "${FILT}" \
                        --disable_adapter_trimming \
                        --average_qual "${Q_MIN}" --length_required "${LEN_MIN}" \
                        --json "${OUTDIR}/${SRR}_fastp.json" --html "${OUTDIR}/${SRR}_fastp.html" \
                        --thread "${THREADS}" > "${LOG}" 2>&1 ;;
            chopper)
                measure_and_run "${TIMELOG}" -- \
                    bash -c 'gunzip -c "$1" | chopper -q "$2" -l "$3" --threads "$4" | gzip > "$5"' _ \
                    "${READS}" "${Q_MIN}" "${LEN_MIN}" "${THREADS}" "${FILT}" > "${LOG}" 2>&1 ;;
        esac
        EXIT_CODE=$?
        parse_time_metrics "${TIMELOG}"

        N_FILT=$(reads_in_fastq "${FILT}")
        VALID=0; METRIC="n/a"
        if [ -s "${FILT}" ] && gzip -t "${FILT}" 2>/dev/null && [ "${N_FILT}" -gt 0 ]; then
            VALID=1; METRIC="${N_FILT} reads passed (>=${LEN_MIN}bp, >=Q${Q_MIN})"
        fi
        append_summary_row "download_qc_pacbio" "${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"

        # --- Kraken2 host/bacterial removal on this filter's output (single-end). ---
        if [ ! -d "${KRAKEN2_DB}" ]; then
            echo "NOTE: Kraken2 DB ${KRAKEN2_DB} absent, skipping host removal for ${TOOL}/${SRR}." >&2
        elif [ "${VALID}" = "1" ]; then
            KOUT="${RESULTS_DIR}/kraken2_${TOOL}_out"
            KTIME="${RESULTS_DIR}/kraken2_${TOOL}_${SRR}.time"
            KLOG="${RESULTS_DIR}/kraken2_${TOOL}_${SRR}.log"
            echo "=== Kraken2 host removal after ${TOOL} on ${SRR} ==="
            measure_and_run "${KTIME}" -- \
                "${REPO_ROOT}/scripts/utils/kraken2_filter_reads_se.sh" \
                "${FILT}" "${KOUT}" "${SRR}" "${KRAKEN2_DB}" > "${KLOG}" 2>&1
            KEXIT=$?
            parse_time_metrics "${KTIME}"
            KFILT="${KOUT}/${SRR}.kraken_filtered.fastq.gz"
            N_KRAK=$(reads_in_fastq "${KFILT}")
            KVALID=0; KMETRIC="n/a"
            [ -s "${KFILT}" ] && [ "${N_KRAK}" -gt 0 ] && { KVALID=1; KMETRIC="${N_KRAK} reads retained after host removal"; }
            append_summary_row "download_qc_pacbio" "kraken2_after_${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${KEXIT}" "${KVALID}" "${KMETRIC}"
        fi
    done

    # --- Dedup comparison: fastp --dedup vs seqkit rmdup, on the fastp-filtered
    # + Kraken2-cleaned reads (the canonical clean input for assembly). ---
    CLEAN="${RESULTS_DIR}/kraken2_fastp_out/${SRR}.kraken_filtered.fastq.gz"
    [ -s "${CLEAN}" ] || CLEAN="${RESULTS_DIR}/fastp_out/${SRR}.filtered.fastq.gz"
    if [ -s "${CLEAN}" ]; then
        N_BEFORE=$(reads_in_fastq "${CLEAN}")
        for DTOOL in fastp seqkit; do
            DOUT="${RESULTS_DIR}/dedup_${DTOOL}_out"; mkdir -p "${DOUT}"
            DFILE="${DOUT}/${SRR}.dedup.fastq.gz"
            DTIME="${RESULTS_DIR}/dedup_${DTOOL}_${SRR}.time"
            DLOG="${RESULTS_DIR}/dedup_${DTOOL}_${SRR}.log"
            echo "=== dedup (${DTOOL}) on ${SRR} ==="
            case "${DTOOL}" in
                fastp)
                    measure_and_run "${DTIME}" -- \
                        fastp -i "${CLEAN}" -o "${DFILE}" --dedup \
                            --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering \
                            --json "${DOUT}/${SRR}_dedup.json" --html "${DOUT}/${SRR}_dedup.html" \
                            --thread "${THREADS}" > "${DLOG}" 2>&1 ;;
                seqkit)
                    measure_and_run "${DTIME}" -- \
                        bash -c 'seqkit rmdup -s "$1" -o "$2"' _ "${CLEAN}" "${DFILE}" > "${DLOG}" 2>&1 ;;
            esac
            DEXIT=$?
            parse_time_metrics "${DTIME}"
            N_AFTER=$(reads_in_fastq "${DFILE}")
            DVALID=0; DMETRIC="n/a"
            [ -s "${DFILE}" ] && [ "${N_AFTER}" -gt 0 ] && { DVALID=1; DMETRIC="${N_AFTER}/${N_BEFORE} reads kept ($((N_BEFORE-N_AFTER)) duplicates removed)"; }
            append_summary_row "download_qc_pacbio" "dedup_${DTOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${DEXIT}" "${DVALID}" "${DMETRIC}"
        done
    fi
done

echo "Done. See ${SUMMARY_TSV}"
