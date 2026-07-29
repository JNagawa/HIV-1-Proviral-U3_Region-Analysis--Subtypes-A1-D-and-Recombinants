#!/bin/bash
# Download-QC step of the tool-comparison harness, Nanopore side. Unlike
# download_qc/illumina (fastp vs Trimmomatic), the tools-review document
# offers no competing alternative to Porechop_ABI (adapter trimming) or
# NanoFilt (quality/length filtering) for Nanopore data -- so this is a
# validity check of the single recommended chain, not a "tool A vs tool B"
# comparison. Runs NanoPlot/NanoQC/NanoStat (pre) -> Porechop_ABI -> NanoFilt
# -> NanoPlot/NanoStat (post) -> MultiQC on the Nanopore sample subset.
# Usage: ./download_qc_oxfordnano.sh
# -u errors on unset vars, pipefail fails a pipe if any stage fails; no -e so one tool failing
# doesn't kill the whole run
set -uo pipefail

# absolute path of this script's own dir, so paths work regardless of launch location
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
# repo top-level (three dirs up), the base for every other path
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
# load shared helpers: measure_and_run, parse_time_metrics, append_summary_row, subset_accessions
source "${REPO_ROOT}/scripts/common/lib_compare.sh"

RAW_DIR="${REPO_ROOT}/data/raw/oxnano"               # input: downloaded Nanopore FASTQs live here
# output: all QC results + summary go here
RESULTS_DIR="${REPO_ROOT}/results/download_qc/oxfordnano"
# the single TSV every run appends a timing/validity row to
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
# create the results dir (and parents) if it doesn't exist
mkdir -p "${RESULTS_DIR}"
# seed the manual notes file from the template only on the first run
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

# Same thresholds as the production pipeline (scripts/pipelines/oxnano_u3analysis.sh) --
# MIN_QUALITY=7 is deliberately more permissive than Illumina's Q20 bar
# (see download_qc/illumina): Nanopore's per-base error rate is inherently
# higher, so a Q20 floor here would discard nearly all reads. MIN_LENGTH=200
# keeps short/degraded fragments out while still retaining partial reads
# (a stricter near-full-length threshold is a legacy/archived option --
# see scripts/archive/README.md's "Open discrepancy" note).
# NanoFilt mean-quality floor (Q7, permissive for noisier Nanopore reads)
MIN_QUALITY=7
# NanoFilt minimum read length (bp), drops short/degraded fragments
MIN_LENGTH=200
# thread count; honour an externally-set THREADS, else default to 4
THREADS="${THREADS:-4}"
export THREADS                                       # export so child scripts/tools inherit it

# loop over just the Nanopore accessions chosen for this comparison
for SRR in $(subset_accessions nanopore "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do
    # expected raw single-end FASTQ for this accession
    FASTQ="${RAW_DIR}/${SRR}.fastq.gz"
    if [ ! -s "${FASTQ}" ]; then                     # if it's missing or empty...
        echo "WARNING: raw FASTQ missing for ${SRR}, skipping." >&2  # ...warn on stderr...
        continue                                     # ...and skip to the next accession
    fi

    OUTDIR="${RESULTS_DIR}/porechop_nanofilt_out"    # output dir for the trimmed+filtered reads
    # timing file for the porechop+nanofilt run
    TIMELOG="${RESULTS_DIR}/porechop_nanofilt_${SRR}.time"
    LOG="${RESULTS_DIR}/porechop_nanofilt_${SRR}.log"  # captured stdout+stderr of the run
    mkdir -p "${OUTDIR}"                             # ensure the output dir exists
    # intermediate adapter-trimmed reads (removed later)
    TRIMMED="${OUTDIR}/${SRR}_trimmed.fastq.gz"
    FILTERED="${OUTDIR}/${SRR}_filtered.fastq.gz"    # final quality/length-filtered reads

    # progress marker for the pre-filter QC block
    echo "=== NanoPlot/NanoQC/NanoStat (pre) on ${SRR} ==="
    # length/quality plots of the raw reads (baseline)
    NanoPlot --fastq "${FASTQ}" --outdir "${RESULTS_DIR}/nanoplot_pre_out/${SRR}" --prefix "${SRR}_pre_" \
        --threads "${THREADS}" --loglength --plots dot --title "${SRR} - Pre-filtering QC" \
        > "${RESULTS_DIR}/nanoplot_pre_${SRR}.log" 2>&1
    # per-base quality across read positions on raw reads
    nanoQC -o "${RESULTS_DIR}/nanoqc_out/${SRR}" "${FASTQ}" > "${RESULTS_DIR}/nanoqc_${SRR}.log" 2>&1
    # text summary stats of the raw reads
    NanoStat --fastq "${FASTQ}" --outdir "${RESULTS_DIR}/nanostat_out" --name "${SRR}_pre_stats.txt" \
        --threads "${THREADS}" > "${RESULTS_DIR}/nanostat_pre_${SRR}.log" 2>&1

    echo "=== Porechop_ABI + NanoFilt on ${SRR} ==="  # progress marker for the trim+filter block
    # adapter-trim (porechop_abi, or porechop, or raw) then NanoFilt quality/length-filter, all
    # timed as one step
    measure_and_run "${TIMELOG}" -- bash -c '
        set -uo pipefail
        fastq="$1" trimmed="$2" filtered="$3" min_quality="$4" min_length="$5" threads="$6"
        if porechop_abi --input "${fastq}" --output "${trimmed}" --threads "${threads}"; then
            :
        elif porechop --input "${fastq}" --output "${trimmed}" --threads "${threads}"; then
            :
        else
            echo "WARNING: adapter trimming unavailable, using raw reads for filtering" >&2
            trimmed="${fastq}"
        fi
        gunzip -c "${trimmed}" | NanoFilt --quality "${min_quality}" --length "${min_length}" | gzip > "${filtered}"
    ' _ "${FASTQ}" "${TRIMMED}" "${FILTERED}" "${MIN_QUALITY}" "${MIN_LENGTH}" "${THREADS}" \
        > "${LOG}" 2>&1
    # capture the pipeline's exit status before $? is overwritten
    EXIT_CODE=$?
    # set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file
    parse_time_metrics "${TIMELOG}"
    # drop the intermediate trimmed file; only the filtered output is kept
    rm -f "${TRIMMED}"

    VALID=0                                          # assume invalid until proven otherwise
    METRIC="n/a"                                     # default key metric
    # filtered output must be non-empty and a valid gzip
    if [ -s "${FILTERED}" ] && gzip -t "${FILTERED}" 2>/dev/null; then
        # read count before filtering (lines/4)
        RAW_READS=$(zcat "${FASTQ}" | awk 'END{print NR/4}')
        FILT_READS=$(zcat "${FILTERED}" | awk 'END{print NR/4}')  # read count after filtering
        # at least one read must survive for a valid result
        if [ "${FILT_READS}" -gt 0 ]; then
            VALID=1                                  # mark the output valid
            # percent of reads retained
            RETAINED=$(awk -v f="${FILT_READS}" -v r="${RAW_READS}" 'BEGIN{printf "%.1f", 100*f/r}')
            # human-readable key metric
            METRIC="${FILT_READS}/${RAW_READS} reads retained (${RETAINED}%)"
        fi

        # progress marker for the post-filter QC block
        echo "=== NanoPlot/NanoStat (post) on ${SRR} ==="
        # length/quality plots of the filtered reads (for before/after)
        NanoPlot --fastq "${FILTERED}" --outdir "${RESULTS_DIR}/nanoplot_post_out/${SRR}" --prefix "${SRR}_post_" \
            --threads "${THREADS}" --loglength --plots dot --title "${SRR} - Post-filtering QC" \
            > "${RESULTS_DIR}/nanoplot_post_${SRR}.log" 2>&1
        # text summary stats of the filtered reads
        NanoStat --fastq "${FILTERED}" --outdir "${RESULTS_DIR}/nanostat_out" --name "${SRR}_post_stats.txt" \
            --threads "${THREADS}" > "${RESULTS_DIR}/nanostat_post_${SRR}.log" 2>&1
    fi

    # write this sample's row to summary.tsv
    append_summary_row "download_qc_oxfordnano" "porechop_abi+nanofilt" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"
done

echo "=== MultiQC ==="                               # progress marker for the aggregation step
# aggregate all per-sample QC into one report; --force overwrites a prior run
multiqc "${RESULTS_DIR}" --outdir "${RESULTS_DIR}/multiqc_out" --filename "nanopore_qc_report" \
    --title "PRJNA765218 - Nanopore QC Summary (subset)" --force > "${RESULTS_DIR}/multiqc.log" 2>&1

# final confirmation pointing the user at the results table
echo "Done. See ${SUMMARY_TSV}"
