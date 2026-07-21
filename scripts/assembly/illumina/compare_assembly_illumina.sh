#!/bin/bash
# Assembly step of the tool-comparison harness: BWA+bcftools-consensus vs
# SPAdes vs SHIVER on the Illumina subset -- one script for all three
# assemblers and their comparison.
# Uses each sample's fastp-trimmed reads from download_qc if available, 
# else falls back to raw reads with a warning.
# Usage: ./compare_assembly_illumina.sh
set -uo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
source "${REPO_ROOT}/scripts/common/lib_compare.sh"

RAW_DIR="${REPO_ROOT}/data/raw/illumina"
FASTP_DIR="${REPO_ROOT}/results/download_qc/illumina/fastp_out"
REF_FASTA="${REPO_ROOT}/data/reference/K03455.1.fasta"
RESULTS_DIR="${REPO_ROOT}/results/assembly/illumina"
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
mkdir -p "${RESULTS_DIR}"
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

# bwa mem requires the reference pre-indexed (bwa index) -- unlike
# minimap2/SPAdes, it errors out immediately without it, so we index it once here

if [ ! -s "${REF_FASTA}.bwt" ]; then
    bwa index "${REF_FASTA}" > "${RESULTS_DIR}/bwa_index.log" 2>&1
fi

THREADS="${THREADS:-4}"
export THREADS

SHIVER_SETUP_DIR="${STAGE_DIR}/shiver_setup"
SHIVER_BIN="${REPO_ROOT}/scripts/tools/shiver/bin"
SHIVER_DATA_DIR="${REPO_ROOT}/scripts/tools/shiver/data/example_input"

# --- One function per assembler. Each is `export -f`'d below so it can be
# invoked as `bash -c '<func> "$@"' _ args...` -- measure_and_run wraps its
# command in /usr/bin/time, which (like any external-process wrapper) can
# only time a real executable, not a shell function living in this script's
# own process, so this is how each function still gets its own wall-clock
# and peak-RSS measurement. ---

# Reference-mapping assembly: the approach currently used by the production
# pipeline (illumina_u3analysis.sh Step 7) -- bwa mem to HXB2 + bcftools
# consensus. The baseline SPAdes/SHIVER are compared against.
run_bwa_consensus() {
    local srr="$1" r1="$2" r2="$3" ref_fasta="$4" outdir="$5"
    mkdir -p "${outdir}"
    local bam="${outdir}/${srr}.sorted.bam" vcf="${outdir}/${srr}.vcf.gz" consensus="${outdir}/${srr}_consensus.fasta"

    bwa mem -t "${THREADS:-4}" "${ref_fasta}" "${r1}" "${r2}" 2>"${outdir}/${srr}_bwa.log" | \
        samtools view -b - | samtools sort -o "${bam}" || exit 1
    samtools index "${bam}" || exit 1

    bcftools mpileup -Ou -f "${ref_fasta}" "${bam}" 2>"${outdir}/${srr}_bcftools.log" | \
        bcftools call -c -Oz -o "${vcf}" 2>>"${outdir}/${srr}_bcftools.log" || exit 1
    tabix -p vcf "${vcf}" || exit 1
    cat "${ref_fasta}" | bcftools consensus "${vcf}" > "${consensus}" 2>>"${outdir}/${srr}_bcftools.log" || exit 1
    sed -i "1s/.*/>${srr}/" "${consensus}"
}
export -f run_bwa_consensus

# De novo assembly with SPAdes, then pick the single best-matching contig
# against HXB2 (by blast) as the "assembly" output for comparison.
run_spades() {
    local srr="$1" r1="$2" r2="$3" ref_fasta="$4" outdir="$5"
    local spades_dir="${outdir}/${srr}_spades"
    mkdir -p "${outdir}"

    spades.py --careful -1 "${r1}" -2 "${r2}" -o "${spades_dir}" -t "${THREADS:-4}" \
        > "${outdir}/${srr}_spades.log" 2>&1
    local contigs="${spades_dir}/contigs.fasta"
    [ -s "${contigs}" ] || { echo "ERROR: SPAdes produced no contigs for ${srr}" >&2; exit 1; }

    # Pick the contig with the best blast hit against HXB2 as the
    # near-full-genome candidate (SPAdes commonly returns many
    # short/host-contaminant contigs).
    if command -v makeblastdb >/dev/null 2>&1 && command -v blastn >/dev/null 2>&1; then
        makeblastdb -in "${ref_fasta}" -dbtype nucl -out "${outdir}/${srr}_hxb2db" >/dev/null 2>&1
        local best_id
        best_id=$(blastn -query "${contigs}" -db "${outdir}/${srr}_hxb2db" -outfmt "6 qseqid length bitscore" 2>/dev/null \
            | sort -k3,3 -rn | head -1 | cut -f1)
        if [ -n "${best_id}" ]; then
            seqkit grep -n -p "${best_id}" "${contigs}" > "${outdir}/${srr}_consensus.fasta"
            sed -i "1s/.*/>${srr}/" "${outdir}/${srr}_consensus.fasta"
            return 0
        fi
    fi

    # Fallback (no blast available): just take the longest contig.
    seqkit sort -l -r "${contigs}" 2>/dev/null | seqkit head -n 1 > "${outdir}/${srr}_consensus.fasta"
    sed -i "1s/.*/>${srr}/" "${outdir}/${srr}_consensus.fasta"
}
export -f run_spades


# SHIVER assembly uses a curated multi-subtype reference to avoid reference bias on
# divergent A1/D/recombinant genomes).
# This handles the fact that Ugandan HIV samples (subtype A1, D, or A/D
# recombinants) differ a lot from HXB2. Instead of forcing every sample onto
# the same reference, it builds a custom reference for each sample from a
# panel of many subtypes so the mapping isn't skewed toward HXB2.
# Three steps run in order:
#   1. SPAdes assembles the reads into contigs
#   2. shiver_align_contigs.sh cleans up the contigs and aligns them to
#      the reference panel
#   3. shiver_map_reads.sh maps the reads to the sample-specific reference
#
# Quirk: the two SHIVER scripts don't take an output directory. They just
# dump files into whatever folder you're standing in. To keep each sample's
# output in its own place, we wrap the whole thing in ( ... ) and 'cd' into
# the sample folder inside it. The parentheses run everything in a subshell,
# so the 'cd' only applies there -- the main script stays put.

run_shiver() {
    local srr="$1" r1="$2" r2="$3" outdir="$4" setup_dir="$5" shiver_bin="$6" shiver_data_dir="$7"
    local ref_alignment="${setup_dir}/HIV1_COM_ref_alignment.fasta"
    if [ ! -s "${ref_alignment}" ]; then
        echo "ERROR: ${ref_alignment} not found. See shiver_setup/SOURCE.md for the one-time manual LANL download step this requires." >&2
        exit 1
    fi

    local adapters="${shiver_data_dir}/adapters_Illumina.fasta"
    local primers="${shiver_data_dir}/primers_GallEtAl2012.fasta"
    local init_dir="${outdir}/shiver_init"
    mkdir -p "${outdir}"

    # Resolve R1/R2 to absolute paths before the subshell cd's, since
    # callers may pass relative paths.
    r1="$(cd "$(dirname "${r1}")" && pwd)/$(basename "${r1}")"
    r2="$(cd "$(dirname "${r2}")" && pwd)/$(basename "${r2}")"

    # shiver_init.sh only needs to run once (its OutDir must not pre-exist)
    if [ ! -d "${init_dir}" ]; then
        "${shiver_bin}/shiver_init.sh" "${init_dir}" "${setup_dir}/config.sh" \
            "${ref_alignment}" "${adapters}" "${primers}" \
            > "${outdir}/shiver_init.log" 2>&1
        if [ $? -ne 0 ]; then
            echo "ERROR: shiver_init.sh failed, see ${outdir}/shiver_init.log" >&2
            exit 1
        fi
    fi

    (
        cd "${outdir}" || exit 1

        # Step 1 of the chain: assemble contigs with SPAdes (cost logged
        # separately so it isn't double-counted against SHIVER's own
        # runtime/RSS when SPAdes is ALSO being benchmarked as its own
        # competing assembler in this same step).
        local contigs_dir="${srr}_spades_for_shiver"
        spades.py --careful -1 "${r1}" -2 "${r2}" -o "${contigs_dir}" -t "${THREADS:-4}" \
            > "${srr}_spades_for_shiver.log" 2>&1
        local contigs="${contigs_dir}/contigs.fasta"
        [ -s "${contigs}" ] || { echo "ERROR: SPAdes (for SHIVER) produced no contigs for ${srr}" >&2; exit 1; }

        # Step 2: align contigs to the reference set, excluding contamination.
        "${shiver_bin}/shiver_align_contigs.sh" "${init_dir}" "${setup_dir}/config.sh" \
            "${contigs}" "${srr}" > "${srr}_align_contigs.log" 2>&1
        if [ $? -ne 0 ]; then
            echo "ERROR: shiver_align_contigs.sh failed, see ${outdir}/${srr}_align_contigs.log" >&2
            exit 1
        fi
        local blast_file="${srr}.blast" cut_wrefs="${srr}_cut_wRefs.fasta"
        [ -s "${blast_file}" ] || { echo "ERROR: no blast hits for ${srr} -- contigs may not be HIV" >&2; exit 1; }

        # Step 3: map reads using the aligned contigs to build the final assembly.
        "${shiver_bin}/shiver_map_reads.sh" "${init_dir}" "${setup_dir}/config.sh" \
            "${contigs}" "${srr}" "${blast_file}" "${cut_wrefs}" "${r1}" "${r2}" \
            > "${srr}_map_reads.log" 2>&1
        if [ $? -ne 0 ]; then
            echo "ERROR: shiver_map_reads.sh failed, see ${outdir}/${srr}_map_reads.log" >&2
            exit 1
        fi

        # SHIVER names its consensus output "<SID>_remap_consensus_MinCov_*.fasta";
        # take the first (default coverage threshold) match.
        local consensus_src
        consensus_src=$(compgen -G "${srr}_remap_consensus_MinCov_*.fasta" | head -1)
        if [ -z "${consensus_src}" ] || [ ! -s "${consensus_src}" ]; then
            echo "ERROR: SHIVER did not produce a consensus FASTA for ${srr}" >&2
            exit 1
        fi
        cp "${consensus_src}" "${srr}_consensus.fasta"
        sed -i "1s/.*/>${srr}/" "${srr}_consensus.fasta"
    )
}
export -f run_shiver

for SRR in $(subset_accessions illumina "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do
    R1="${FASTP_DIR}/${SRR}_1.trimmed.fastq.gz"
    R2="${FASTP_DIR}/${SRR}_2.trimmed.fastq.gz"
    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ]; then
        echo "NOTE: no fastp-trimmed reads for ${SRR}, run scripts/download_qc/illumina/download_qc_illumina.sh first for the intended input. Falling back to raw reads." >&2
        R1="${RAW_DIR}/${SRR}_1.fastq.gz"
        R2="${RAW_DIR}/${SRR}_2.fastq.gz"
    fi
    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ]; then
        echo "WARNING: no reads at all for ${SRR}, skipping." >&2
        continue
    fi

    for TOOL in bwa_consensus spades shiver; do
        OUTDIR="${RESULTS_DIR}/${TOOL}_out"
        TIMELOG="${RESULTS_DIR}/${TOOL}_${SRR}.time"
        LOG="${RESULTS_DIR}/${TOOL}_${SRR}.log"
        CONSENSUS="${OUTDIR}/${SRR}_consensus.fasta"

        echo "=== ${TOOL} on ${SRR} ==="
        case "${TOOL}" in
            bwa_consensus)
                measure_and_run "${TIMELOG}" -- \
                    bash -c 'run_bwa_consensus "$@"' _ "${SRR}" "${R1}" "${R2}" "${REF_FASTA}" "${OUTDIR}" > "${LOG}" 2>&1 ;;
            spades)
                measure_and_run "${TIMELOG}" -- \
                    bash -c 'run_spades "$@"' _ "${SRR}" "${R1}" "${R2}" "${REF_FASTA}" "${OUTDIR}" > "${LOG}" 2>&1 ;;
            shiver)
                measure_and_run "${TIMELOG}" -- \
                    bash -c 'run_shiver "$@"' _ "${SRR}" "${R1}" "${R2}" "${OUTDIR}" "${SHIVER_SETUP_DIR}" "${SHIVER_BIN}" "${SHIVER_DATA_DIR}" > "${LOG}" 2>&1 ;;
        esac
        EXIT_CODE=$?
        parse_time_metrics "${TIMELOG}"

        VALID=0
        METRIC="n/a"
        if [ -s "${CONSENSUS}" ]; then
            LEN=$(seqkit stats -T "${CONSENSUS}" 2>/dev/null | tail -1 | cut -f5)
            N_PCT=$(seqkit fx2tab -n -g -B N "${CONSENSUS}" 2>/dev/null | awk -F'\t' '{print $NF}' | tail -1)
            if [ -n "${LEN}" ] && [ "${LEN}" -ge 8000 ] && [ "${LEN}" -le 10000 ]; then
                VALID=1
            fi
            METRIC="${LEN}bp, ${N_PCT:-?}% N"
        fi

        append_summary_row "assembly_illumina" "${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"
    done
done

echo "Done. See ${SUMMARY_TSV}"
