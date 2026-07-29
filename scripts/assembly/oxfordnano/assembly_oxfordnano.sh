#!/bin/bash
# For the assembly step of  Nanopore sequences, minimap2 and bcft00ls were used.
# reference-mapping via minimap2 + bcftools consensus)
# Usage: ./assembly_oxfordnano.sh
# -u errors on unset vars, pipefail fails a pipe if any stage fails
set -uo pipefail

# absolute path of this script's own dir, so paths work regardless of launch location
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
# repo top-level (three dirs up), the base for every other path below
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
# load shared helpers: measure_and_run, parse_time_metrics, append_summary_row, subset_accessions
source "${REPO_ROOT}/scripts/common/lib_compare.sh"

# preferred input: adapter-trimmed + length/quality-filtered reads
FILTERED_DIR="${REPO_ROOT}/results/download_qc/oxfordnano/porechop_nanofilt_out"
# fallback input: raw Nanopore FASTQs
RAW_DIR="${REPO_ROOT}/data/raw/oxnano"
# HXB2 reference genome used for mapping/consensus
REF_FASTA="${REPO_ROOT}/data/reference/K03455.1.fasta"
# output: all assembly results + summary go here
RESULTS_DIR="${REPO_ROOT}/results/assembly/oxfordnano"
# the single TSV every run appends a timing/validity row to
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
# create the results dir (and parents) if it doesn't exist
mkdir -p "${RESULTS_DIR}"
# seed the manual notes file from the template only on the first run
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

# thread count; honour an externally-set THREADS, otherwise default to 4
THREADS="${THREADS:-4}"
# export so the exported assembler function inherits it
export THREADS

# reference-mapping assembler (minimap2 + bcftools consensus)
run_minimap2_consensus() {
    # positional args: sample id, reads, reference, output dir
    local srr="$1" fastq="$2" ref_fasta="$3" outdir="$4"
    mkdir -p "${outdir}"                             # ensure this sample's output dir exists
    # derived output paths
    local bam="${outdir}/${srr}.sorted.bam" vcf="${outdir}/${srr}.vcf.gz" consensus="${outdir}/${srr}_consensus.fasta"

    # map ONT reads to HXB2 (map-ont preset), convert to BAM, sort; bail if any stage fails
    minimap2 -ax map-ont -t "${THREADS:-4}" "${ref_fasta}" "${fastq}" 2>"${outdir}/${srr}_minimap2.log" | \
        samtools view -b - | samtools sort -o "${bam}" || exit 1
    samtools index "${bam}" || exit 1                # index the sorted BAM (required by mpileup)

    # pile up bases and call variants against the reference
    bcftools mpileup -Ou -f "${ref_fasta}" "${bam}" 2>"${outdir}/${srr}_bcftools.log" | \
        bcftools call -c -Oz -o "${vcf}" 2>>"${outdir}/${srr}_bcftools.log" || exit 1
    # index the VCF so bcftools consensus can read it
    tabix -p vcf "${vcf}" || exit 1
    # apply called variants onto the reference to make the consensus
    cat "${ref_fasta}" | bcftools consensus "${vcf}" > "${consensus}" 2>>"${outdir}/${srr}_bcftools.log" || exit 1
    sed -i "1s/.*/>${srr}/" "${consensus}"           # rename the FASTA header to the sample id
}
# export so measure_and_run's child bash can call it
export -f run_minimap2_consensus

# loop over just the Nanopore accessions chosen for this comparison
for SRR in $(subset_accessions nanopore "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do
    # preferred input: this sample's filtered reads
    FASTQ="${FILTERED_DIR}/${SRR}_filtered.fastq.gz"
    if [ ! -s "${FASTQ}" ]; then                     # if the filtered reads are missing...
        # ...warn...
        echo "NOTE: no filtered reads for ${SRR}, run scripts/download_qc/oxfordnano/download_qc_oxfordnano.sh first. Falling back to raw reads." >&2
        FASTQ="${RAW_DIR}/${SRR}.fastq.gz"           # ...fall back to the raw reads
    fi
    if [ ! -s "${FASTQ}" ]; then                     # if even raw reads are absent...
        echo "WARNING: no reads at all for ${SRR}, skipping." >&2  # ...warn...
        continue                                     # ...and skip this sample entirely
    fi

    OUTDIR="${RESULTS_DIR}/minimap2_out"             # output dir for the minimap2 consensus
    # file where measure_and_run records wallclock/RSS
    TIMELOG="${RESULTS_DIR}/minimap2_${SRR}.time"
    LOG="${RESULTS_DIR}/minimap2_${SRR}.log"         # captured stdout+stderr of the run
    # the consensus FASTA the run is expected to produce
    CONSENSUS="${OUTDIR}/${SRR}_consensus.fasta"

    echo "=== minimap2 consensus on ${SRR} ==="      # progress marker in the log
    # time+run the minimap2 consensus function
    measure_and_run "${TIMELOG}" -- \
        bash -c 'run_minimap2_consensus "$@"' _ "${SRR}" "${FASTQ}" "${REF_FASTA}" "${OUTDIR}" > "${LOG}" 2>&1
    # capture the run's exit status before $? is overwritten
    EXIT_CODE=$?
    # set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file
    parse_time_metrics "${TIMELOG}"

    VALID=0                                           # assume invalid until proven otherwise
    METRIC="n/a"                                      # human-readable key metric, filled in below
    if [ -s "${CONSENSUS}" ]; then                    # only evaluate if a consensus was produced
        # consensus length in bp (col 5 of seqkit stats)
        LEN=$(seqkit stats -T "${CONSENSUS}" 2>/dev/null | tail -1 | cut -f5)
        # percent of N (ambiguous) bases
        N_PCT=$(seqkit fx2tab -n -g -B N "${CONSENSUS}" 2>/dev/null | awk -F'\t' '{print $NF}' | tail -1)
        # valid if length is near a full HIV genome (~9kb)
        if [ -n "${LEN}" ] && [ "${LEN}" -ge 8000 ] && [ "${LEN}" -le 10000 ]; then
            VALID=1                                   # mark valid
        fi
        METRIC="${LEN}bp, ${N_PCT:-?}% N"             # record length and N% as the key metric
    fi

    # write this run's row to summary.tsv
    append_summary_row "assembly_oxfordnano" "minimap2" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"
done

# final confirmation pointing the user at the results table
echo "Done. See ${SUMMARY_TSV}"
