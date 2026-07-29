#!/bin/bash
# Download-QC step of the tool-comparison harness: runs FastQC, fastp, and
# Trimmomatic directly, then compares fastp vs Trimmomatic on the Illumina
# sample subset. Both trimmers are held to the same quality/length bar (see
# TRIM_* below) so the comparison reflects the algorithms, not mismatched
# settings -- matches scripts/pipelines/illumina/01_download_qc_trim.sh's
# production thresholds.
# Usage: ./download_qc_illumina.sh   (run from this directory, with HIV_U3analysis activated)
# -u errors on unset vars, pipefail fails a pipe if any stage fails; no -e so one tool failing
# doesn't kill the whole comparison
set -uo pipefail

# absolute path of this script's own dir, so paths work regardless of launch location
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
# repo top-level (three dirs up), the base for every other path
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
# load shared helpers: measure_and_run, parse_time_metrics, append_summary_row, subset_accessions
source "${REPO_ROOT}/scripts/common/lib_compare.sh"

# input: downloaded paired FASTQs live here
RAW_DIR="${REPO_ROOT}/data/raw/illumina"
# output: all QC results + summary go here
RESULTS_DIR="${REPO_ROOT}/results/download_qc/illumina"
# the single TSV every tool appends a timing/validity row to
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
# Kraken2 database used for host/bacterial read removal
KRAKEN2_DB="${REPO_ROOT}/data/reference/kraken2_standard_16gb_db"
# create the results dir (and parents) if it doesn't exist
mkdir -p "${RESULTS_DIR}"
# seed the manual notes file from the template only on the first run
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

# Phred quality score Q is defined as Q = -10*log10(P_error) (Ewing & Green
# 1998) -- Q15 corresponds to a 1-in-32 (~97%) per-base call accuracy, Q20
# to 1-in-100 (99%). Q15 is used below (sliding-window trim AND whole-read
# average) as the quality bar, matching Trimmomatic's own textbook
# SLIDINGWINDOW:4:15 example -- both trimmers are held to the same
# threshold so the comparison reflects the algorithms, not mismatched
# settings.

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
# fastp's --cut_right/--cut_right_window_size/--cut_right_mean_quality below
# reproduce this exact same window+threshold for a fair comparison.
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
# OVERALL average quality to be >=Q15 after trimming. A read can pass a 4bp
# sliding window check while still being poor quality on average; this is
# the final whole-read quality gate applied to both trimmers.
# whole-read mean-quality floor (Q15), applied to both trimmers
TRIM_AVGQUAL=15

# locate the TruSeq3 adapter FASTA in the conda env (path varies by version)
ADAPTER_FILE=$(compgen -G "${CONDA_PREFIX}/share/trimmomatic*/adapters/TruSeq3-PE-2.fa" | head -1)
# thread count; honour an externally-set THREADS, else default to 4
THREADS="${THREADS:-4}"
export THREADS                                       # export so child scripts/tools inherit it

# loop over just the Illumina accessions chosen for this comparison
for SRR in $(subset_accessions illumina "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do
    R1="${RAW_DIR}/${SRR}_1.fastq.gz"                # forward-mate raw FASTQ for this accession
    R2="${RAW_DIR}/${SRR}_2.fastq.gz"                # reverse-mate raw FASTQ for this accession
    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ]; then    # if either mate is missing/empty...
        echo "WARNING: raw FASTQs missing for ${SRR}, skipping." >&2  # ...warn on stderr...
        continue                                     # ...and skip to the next accession
    fi

    # --- Pre-trim FastQC: one shared baseline per sample, not tied to
    # either trimmer (there's only one "before", regardless of which
    # trimmer runs next). ---
    # shared output dir for the pre-trim baseline QC
    FASTQC_PRE_DIR="${RESULTS_DIR}/fastqc_pre_out"
    mkdir -p "${FASTQC_PRE_DIR}"                     # ensure it exists
    # skip if the report already exists (makes reruns idempotent)
    if [ ! -s "${FASTQC_PRE_DIR}/${SRR}_1_fastqc.html" ]; then
        echo "=== FastQC (pre-trim) on ${SRR} ==="   # progress marker in the log
        # QC both raw mates; redirect stdout+stderr to the log
        fastqc "${R1}" "${R2}" --outdir "${FASTQC_PRE_DIR}" --threads "${THREADS}" --quiet \
            > "${RESULTS_DIR}/fastqc_pre_${SRR}.log" 2>&1
    fi

    # run each trimmer on the same input for a head-to-head comparison
    for TOOL in fastp trimmomatic; do
        OUTDIR="${RESULTS_DIR}/${TOOL}_out"          # per-tool output dir (e.g. fastp_out/)
        # file where measure_and_run records wallclock/RSS
        TIMELOG="${RESULTS_DIR}/${TOOL}_${SRR}.time"
        LOG="${RESULTS_DIR}/${TOOL}_${SRR}.log"      # captured stdout+stderr of the tool
        mkdir -p "${OUTDIR}"                         # create the per-tool output dir if needed

        echo "=== ${TOOL} on ${SRR} ==="            # progress marker in the log
        # fastp branch (single command with equivalent thresholds)
        if [ "${TOOL}" = "fastp" ]; then
            # fastp's cut_right is off by default -- enable it so it does a
            # real sliding-window trim equivalent to Trimmomatic's
            # SLIDINGWINDOW, rather than being compared on easier defaults.
            # trim + dedup with the shared Q15/len50 bar; emit JSON/HTML reports
            measure_and_run "${TIMELOG}" -- \
                fastp \
                    -i "${R1}" -I "${R2}" \
                    -o "${OUTDIR}/${SRR}_1.trimmed.fastq.gz" -O "${OUTDIR}/${SRR}_2.trimmed.fastq.gz" \
                    --cut_right --cut_right_window_size 4 --cut_right_mean_quality "${TRIM_AVGQUAL}" \
                    --average_qual "${TRIM_AVGQUAL}" \
                    --length_required "${TRIM_MINLEN}" \
                    --dedup \
                    --json "${OUTDIR}/${SRR}_fastp.json" --html "${OUTDIR}/${SRR}_fastp.html" \
                    --thread "${THREADS}" \
                > "${LOG}" 2>&1
        else                                         # trimmomatic branch (builds a step list first)
            # accumulate Trimmomatic's ordered trimming steps here
            TRIM_STEPS=""
            # only add adapter clipping if the adapter FASTA was found
            if [ -n "${ADAPTER_FILE}" ] && [ -f "${ADAPTER_FILE}" ]; then
                # TruSeq3 PE adapter-clip step (recommended params)
                TRIM_STEPS="ILLUMINACLIP:${ADAPTER_FILE}:2:30:10:2:True "
            fi
            # append leading/trailing edge trims
            TRIM_STEPS+="LEADING:${TRIM_LEADING} TRAILING:${TRIM_TRAILING} "
            # append the 4bp/Q15 sliding-window trim
            TRIM_STEPS+="SLIDINGWINDOW:${TRIM_SLIDINGWINDOW} "
            # append the whole-read Q15 average-quality gate
            TRIM_STEPS+="AVGQUAL:${TRIM_AVGQUAL} "
            # append the 50bp minimum-length filter (must come last)
            TRIM_STEPS+="MINLEN:${TRIM_MINLEN}"

            # run Trimmomatic PE with the built step list; -summary feeds the metric below
            measure_and_run "${TIMELOG}" -- \
                trimmomatic PE \
                    -Xmx48g \
                    -threads "${THREADS}" \
                    -phred33 \
                    -summary "${OUTDIR}/${SRR}_trimmomatic_summary.txt" \
                    "${R1}" "${R2}" \
                    "${OUTDIR}/${SRR}_1.trimmed.fastq.gz" "${OUTDIR}/${SRR}_1.unpaired.fastq.gz" \
                    "${OUTDIR}/${SRR}_2.trimmed.fastq.gz" "${OUTDIR}/${SRR}_2.unpaired.fastq.gz" \
                    ${TRIM_STEPS} \
                > "${LOG}" 2>&1
        fi
        # capture the tool's exit status before $? is overwritten
        EXIT_CODE=$?
        # set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file
        parse_time_metrics "${TIMELOG}"

        # Validity check: R1/R2 output pair counts match and are non-zero;
        # pull the relevant surviving-reads metric from each tool's own report.
        R1_OUT="${OUTDIR}/${SRR}_1.trimmed.fastq.gz" # this tool's trimmed forward mate
        R2_OUT="${OUTDIR}/${SRR}_2.trimmed.fastq.gz" # this tool's trimmed reverse mate
        VALID=0                                      # assume invalid until proven otherwise
        METRIC="n/a"                                 # default key metric
        # both mates non-empty and valid gzip
        if [ -s "${R1_OUT}" ] && [ -s "${R2_OUT}" ] && gzip -t "${R1_OUT}" 2>/dev/null && gzip -t "${R2_OUT}" 2>/dev/null; then
            R1_READS=$(zcat "${R1_OUT}" | wc -l)     # line count of forward mate (reads = lines/4)
            R2_READS=$(zcat "${R2_OUT}" | wc -l)     # line count of reverse mate
            # mates must be equal and non-empty for a valid pair
            if [ "${R1_READS}" = "${R2_READS}" ] && [ "${R1_READS}" -gt 0 ]; then
                VALID=1                              # mark this tool's output valid
            fi
            # fastp reports passed-reads in its JSON
            if [ "${TOOL}" = "fastp" ] && [ -s "${OUTDIR}/${SRR}_fastp.json" ]; then
                # pull the first passed_filter_reads number
                PASSED=$(grep -o '"passed_filter_reads"[^,}]*' "${OUTDIR}/${SRR}_fastp.json" | head -1 | grep -o '[0-9]*$')
                # record it as the key metric if found
                [ -n "${PASSED}" ] && METRIC="${PASSED} reads passed"
            # Trimmomatic reports surviving % in its summary
            elif [ "${TOOL}" = "trimmomatic" ] && [ -f "${OUTDIR}/${SRR}_trimmomatic_summary.txt" ]; then
                # extract the both-surviving percentage
                METRIC=$(grep "Both Surviving Read Percent" "${OUTDIR}/${SRR}_trimmomatic_summary.txt" | awk '{print $NF"% surviving"}')
            fi

            # --- Post-trim FastQC on this tool's own output ---
            # per-tool post-trim QC output dir
            FASTQC_POST_DIR="${RESULTS_DIR}/fastqc_post_${TOOL}_out"
            mkdir -p "${FASTQC_POST_DIR}"            # ensure it exists
            # QC this tool's trimmed reads for a before/after comparison
            fastqc "${R1_OUT}" "${R2_OUT}" --outdir "${FASTQC_POST_DIR}" --threads "${THREADS}" --quiet \
                > "${RESULTS_DIR}/fastqc_post_${TOOL}_${SRR}.log" 2>&1
        fi

        # write this trimmer's row to summary.tsv
        append_summary_row "download_qc_illumina" "${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"

        # --- Deduplication + Kraken2 host/bacterial contamination filtering,
        # applied downstream of whichever trimmer just ran (fastp's own
        # --dedup above already deduplicates; Trimmomatic has no native
        # dedup capability, so a fastp dedup-only pass runs first for that
        # arm). This mirrors the review's fastp -> Kraken2 -> SHIVER chain
        # for both trimmer arms, so the comparison stays fair. ---
        # only run downstream filtering when this trimmer produced usable reads
        if [ "${VALID}" -eq 1 ]; then
            # default: feed the trimmer output straight in (fastp already deduped)
            DEDUP_R1="${R1_OUT}" DEDUP_R2="${R2_OUT}"
            # Trimmomatic has no dedup, so add a fastp dedup-only pass first
            if [ "${TOOL}" = "trimmomatic" ]; then
                # deduped forward mate for the Trimmomatic arm
                DEDUP_R1="${OUTDIR}/${SRR}_1.dedup.fastq.gz"
                # deduped reverse mate for the Trimmomatic arm
                DEDUP_R2="${OUTDIR}/${SRR}_2.dedup.fastq.gz"
                # run the fastp dedup-only wrapper on Trimmomatic's output
                bash "${REPO_ROOT}/scripts/utils/fastp_dedup.sh" \
                    "${R1_OUT}" "${R2_OUT}" "${DEDUP_R1}" "${DEDUP_R2}" \
                    "${OUTDIR}/${SRR}_dedup_fastp" \
                    > "${RESULTS_DIR}/dedup_${TOOL}_${SRR}.log" 2>&1
            fi

            # per-trimmer Kraken2 output dir
            KRAKEN_OUTDIR="${RESULTS_DIR}/kraken2_${TOOL}_out"
            # timing file for this Kraken2 run
            KRAKEN_TIMELOG="${RESULTS_DIR}/kraken2_${TOOL}_${SRR}.time"
            # paired-end Kraken2 wrapper: drop host/bacterial reads
            measure_and_run "${KRAKEN_TIMELOG}" -- \
                bash "${REPO_ROOT}/scripts/utils/kraken2_filter_reads.sh" \
                    "${DEDUP_R1}" "${DEDUP_R2}" "${KRAKEN_OUTDIR}" "${SRR}" "${KRAKEN2_DB}" \
                > "${RESULTS_DIR}/kraken2_${TOOL}_${SRR}.log" 2>&1
            KRAKEN_EXIT=$?                           # capture Kraken2's exit status
            # refresh WALLCLOCK_SEC / PEAK_RSS_MB from the Kraken2 timing file
            parse_time_metrics "${KRAKEN_TIMELOG}"

            KRAKEN_VALID=0                           # default to invalid until checked
            KRAKEN_METRIC="n/a"                      # default key metric
            # host-removed forward mate
            KRAKEN_R1_OUT="${KRAKEN_OUTDIR}/${SRR}_1.kraken_filtered.fastq.gz"
            # host-removed reverse mate
            KRAKEN_R2_OUT="${KRAKEN_OUTDIR}/${SRR}_2.kraken_filtered.fastq.gz"
            # both filtered mates present and non-empty
            if [ -s "${KRAKEN_R1_OUT}" ] && [ -s "${KRAKEN_R2_OUT}" ]; then
                KRAKEN_VALID=1                       # mark the Kraken2 output valid
                # pull the retained-reads summary line from the log
                KRAKEN_METRIC=$(grep -o '^Kraken2 filtering.*retained[^.]*' "${RESULTS_DIR}/kraken2_${TOOL}_${SRR}.log" | tail -1)
                # fall back to pointing at the kreport if that line is absent
                [ -z "${KRAKEN_METRIC}" ] && KRAKEN_METRIC="filtered, see ${KRAKEN_OUTDIR}/${SRR}.kreport"
            fi
            # record the Kraken2 row
            append_summary_row "download_qc_illumina" "kraken2_after_${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${KRAKEN_EXIT}" "${KRAKEN_VALID}" "${KRAKEN_METRIC}"
        fi
    done
done

# final confirmation pointing the user at the results table
echo "Done. See ${SUMMARY_TSV}"
