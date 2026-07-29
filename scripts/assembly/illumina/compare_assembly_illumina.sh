#!/bin/bash
# Assembly step of the tool-comparison harness: BWA+bcftools-consensus vs
# SPAdes vs SHIVER on the Illumina subset -- one script for all three
# assemblers and their comparison.
# Uses each sample's fastp-trimmed reads from download_qc if available,
# else falls back to raw reads with a warning.
# Usage: ./compare_assembly_illumina.sh
# -u errors on unset vars, pipefail fails a pipe if any stage fails; no -e so one assembler failing
# doesn't kill the comparison
set -uo pipefail

# absolute path of this script's own dir, so paths work regardless of launch location
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
# repo top-level (three dirs up), the base for every other path below
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
# load shared helpers: measure_and_run, parse_time_metrics, append_summary_row, subset_accessions
source "${REPO_ROOT}/scripts/common/lib_compare.sh"

# fallback input: raw downloaded Illumina FASTQs
RAW_DIR="${REPO_ROOT}/data/raw/illumina"
# preferred input: fastp-trimmed reads from the QC step
FASTP_DIR="${REPO_ROOT}/results/download_qc/illumina/fastp_out"
# HXB2 reference genome used for mapping/consensus
REF_FASTA="${REPO_ROOT}/data/reference/K03455.1.fasta"
# output: all assembly results + summary go here
RESULTS_DIR="${REPO_ROOT}/results/assembly/illumina"
# the single TSV every assembler appends a timing/validity row to
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
# create the results dir (and parents) if it doesn't exist
mkdir -p "${RESULTS_DIR}"
# seed the manual notes file from the template only on the first run
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

# bwa mem requires the reference pre-indexed (bwa index) -- unlike
# minimap2/SPAdes, it errors out immediately without it, so we index it once here

# skip indexing if the .bwt index already exists (idempotent reruns)
if [ ! -s "${REF_FASTA}.bwt" ]; then
    # build the bwa index once, capturing output to a log
    bwa index "${REF_FASTA}" > "${RESULTS_DIR}/bwa_index.log" 2>&1
fi

# thread count; honour an externally-set THREADS, otherwise default to 4
THREADS="${THREADS:-4}"
# export so the exported assembler functions inherit it
export THREADS

# dir holding SHIVER's config.sh and reference alignment
SHIVER_SETUP_DIR="${STAGE_DIR}/shiver_setup"
# dir with the SHIVER executable scripts
SHIVER_BIN="${REPO_ROOT}/scripts/tools/shiver/bin"
# dir with SHIVER's bundled adapters/primers FASTAs
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
# baseline reference-mapping assembler (bwa + bcftools consensus)
run_bwa_consensus() {
    # positional args: sample id, paired reads, reference, output dir
    local srr="$1" r1="$2" r2="$3" ref_fasta="$4" outdir="$5"
    mkdir -p "${outdir}"                             # ensure this sample's output dir exists
    # derived output paths
    local bam="${outdir}/${srr}.sorted.bam" vcf="${outdir}/${srr}.vcf.gz" consensus="${outdir}/${srr}_consensus.fasta"

    # map reads to HXB2, convert to BAM, sort; bail if any stage fails
    bwa mem -t "${THREADS:-4}" "${ref_fasta}" "${r1}" "${r2}" 2>"${outdir}/${srr}_bwa.log" | \
        samtools view -b - | samtools sort -o "${bam}" || exit 1
    samtools index "${bam}" || exit 1                # index the sorted BAM (required by mpileup)

    # pile up bases and call variants against the reference
    bcftools mpileup -Ou -f "${ref_fasta}" "${bam}" 2>"${outdir}/${srr}_bcftools.log" | \
        bcftools call -c -Oz -o "${vcf}" 2>>"${outdir}/${srr}_bcftools.log" || exit 1
    # index the VCF so bcftools consensus can read it
    tabix -p vcf "${vcf}" || exit 1
    # apply the called variants onto the reference to make the consensus
    cat "${ref_fasta}" | bcftools consensus "${vcf}" > "${consensus}" 2>>"${outdir}/${srr}_bcftools.log" || exit 1
    sed -i "1s/.*/>${srr}/" "${consensus}"           # rename the FASTA header to the sample id
}
# export so measure_and_run's child bash can call it
export -f run_bwa_consensus

# De novo assembly with SPAdes, then pick the single best-matching contig
# against HXB2 (by blast) as the "assembly" output for comparison.
# de novo assembler (SPAdes) with best-contig selection
run_spades() {
    # positional args: sample id, paired reads, reference, output dir
    local srr="$1" r1="$2" r2="$3" ref_fasta="$4" outdir="$5"
    local spades_dir="${outdir}/${srr}_spades"       # per-sample SPAdes working dir
    mkdir -p "${outdir}"                             # ensure output dir exists

    # assemble reads into contigs; --careful reduces mismatches/indels
    spades.py --careful -1 "${r1}" -2 "${r2}" -o "${spades_dir}" -t "${THREADS:-4}" \
        > "${outdir}/${srr}_spades.log" 2>&1
    local contigs="${spades_dir}/contigs.fasta"      # SPAdes' contig output file
    # fail if no contigs were produced
    [ -s "${contigs}" ] || { echo "ERROR: SPAdes produced no contigs for ${srr}" >&2; exit 1; }

    # Pick the contig with the best blast hit against HXB2 as the
    # near-full-genome candidate (SPAdes commonly returns many
    # short/host-contaminant contigs).
    # only do blast selection if the blast tools exist
    if command -v makeblastdb >/dev/null 2>&1 && command -v blastn >/dev/null 2>&1; then
        # build a blast DB from HXB2
        makeblastdb -in "${ref_fasta}" -dbtype nucl -out "${outdir}/${srr}_hxb2db" >/dev/null 2>&1
        local best_id                                # will hold the id of the best-matching contig
        # blast contigs vs HXB2, sort by bitscore, take the top hit's contig id
        best_id=$(blastn -query "${contigs}" -db "${outdir}/${srr}_hxb2db" -outfmt "6 qseqid length bitscore" 2>/dev/null \
            | sort -k3,3 -rn | head -1 | cut -f1)
        if [ -n "${best_id}" ]; then                 # if a best contig was found...
            # extract that one contig as the assembly
            seqkit grep -n -p "${best_id}" "${contigs}" > "${outdir}/${srr}_consensus.fasta"
            # rename the FASTA header to the sample id
            sed -i "1s/.*/>${srr}/" "${outdir}/${srr}_consensus.fasta"
            return 0                                 # done -- skip the length-based fallback
        fi
    fi

    # Fallback (no blast available): just take the longest contig.
    # sort contigs by length desc, keep the longest
    seqkit sort -l -r "${contigs}" 2>/dev/null | seqkit head -n 1 > "${outdir}/${srr}_consensus.fasta"
    # rename the FASTA header to the sample id
    sed -i "1s/.*/>${srr}/" "${outdir}/${srr}_consensus.fasta"
}
# export so measure_and_run's child bash can call it
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

# SHIVER 3-step assembler with per-sample custom reference
run_shiver() {
    # positional args passed in from the loop
    local srr="$1" r1="$2" r2="$3" outdir="$4" setup_dir="$5" shiver_bin="$6" shiver_data_dir="$7"
    # LANL multi-subtype reference alignment panel
    local ref_alignment="${setup_dir}/HIV1_COM_ref_alignment.fasta"
    # the panel is a manual one-time download, so verify it exists
    if [ ! -s "${ref_alignment}" ]; then
        # tell the user how to obtain it
        echo "ERROR: ${ref_alignment} not found. See shiver_setup/SOURCE.md for the one-time manual LANL download step this requires." >&2
        exit 1                                       # cannot proceed without the reference panel
    fi

    # adapter sequences SHIVER trims from reads
    local adapters="${shiver_data_dir}/adapters_Illumina.fasta"
    # PCR primer sequences SHIVER trims from reads
    local primers="${shiver_data_dir}/primers_GallEtAl2012.fasta"
    local init_dir="${outdir}/shiver_init"           # dir for SHIVER's one-time init output
    mkdir -p "${outdir}"                             # ensure this sample's output dir exists

    # Resolve R1/R2 to absolute paths before the subshell cd's, since
    # callers may pass relative paths.
    # canonicalize R1 to an absolute path (survives the later cd)
    r1="$(cd "$(dirname "${r1}")" && pwd)/$(basename "${r1}")"
    # canonicalize R2 to an absolute path (survives the later cd)
    r2="$(cd "$(dirname "${r2}")" && pwd)/$(basename "${r2}")"

    # shiver_init.sh only needs to run once (its OutDir must not pre-exist)
    # only initialise if it hasn't been done yet (dir must not pre-exist)
    if [ ! -d "${init_dir}" ]; then
        # build SHIVER's initialisation dir from the reference/adapters/primers
        "${shiver_bin}/shiver_init.sh" "${init_dir}" "${setup_dir}/config.sh" \
            "${ref_alignment}" "${adapters}" "${primers}" \
            > "${outdir}/shiver_init.log" 2>&1
        if [ $? -ne 0 ]; then                        # if init failed...
            # ...point at the log...
            echo "ERROR: shiver_init.sh failed, see ${outdir}/shiver_init.log" >&2
            exit 1                                   # ...and abort
        fi
    fi

    # subshell so the cd below is local and doesn't move the main script
    (
        # SHIVER writes to cwd, so enter this sample's output dir
        cd "${outdir}" || exit 1

        # Step 1 of the chain: assemble contigs with SPAdes (cost logged
        # separately so it isn't double-counted against SHIVER's own
        # runtime/RSS when SPAdes is ALSO being benchmarked as its own
        # competing assembler in this same step).
        # separate SPAdes dir so its cost is logged apart from SHIVER's
        local contigs_dir="${srr}_spades_for_shiver"
        spades.py --careful -1 "${r1}" -2 "${r2}" -o "${contigs_dir}" -t "${THREADS:-4}" \
            > "${srr}_spades_for_shiver.log" 2>&1     # assemble contigs to feed into SHIVER
        local contigs="${contigs_dir}/contigs.fasta" # the contigs SPAdes produced
        # abort if no contigs
        [ -s "${contigs}" ] || { echo "ERROR: SPAdes (for SHIVER) produced no contigs for ${srr}" >&2; exit 1; }

        # Step 2: align contigs to the reference set, excluding contamination.
        # clean/align contigs to the reference panel
        "${shiver_bin}/shiver_align_contigs.sh" "${init_dir}" "${setup_dir}/config.sh" \
            "${contigs}" "${srr}" > "${srr}_align_contigs.log" 2>&1
        if [ $? -ne 0 ]; then                        # if alignment failed...
            # ...point at the log...
            echo "ERROR: shiver_align_contigs.sh failed, see ${outdir}/${srr}_align_contigs.log" >&2
            exit 1                                   # ...and abort
        fi
        # outputs from the align step needed by the map step
        local blast_file="${srr}.blast" cut_wrefs="${srr}_cut_wRefs.fasta"
        # no hits means contigs aren't HIV
        [ -s "${blast_file}" ] || { echo "ERROR: no blast hits for ${srr} -- contigs may not be HIV" >&2; exit 1; }

        # Step 3: map reads using the aligned contigs to build the final assembly.
        # map reads to the sample-specific reference to build the consensus
        "${shiver_bin}/shiver_map_reads.sh" "${init_dir}" "${setup_dir}/config.sh" \
            "${contigs}" "${srr}" "${blast_file}" "${cut_wrefs}" "${r1}" "${r2}" \
            > "${srr}_map_reads.log" 2>&1
        if [ $? -ne 0 ]; then                        # if mapping failed...
            # ...point at the log...
            echo "ERROR: shiver_map_reads.sh failed, see ${outdir}/${srr}_map_reads.log" >&2
            exit 1                                   # ...and abort
        fi

        # SHIVER names its consensus output "<SID>_remap_consensus_MinCov_*.fasta";
        # take the first (default coverage threshold) match.
        local consensus_src                          # will hold the SHIVER consensus filename
        # glob for SHIVER's consensus output, take the first match
        consensus_src=$(compgen -G "${srr}_remap_consensus_MinCov_*.fasta" | head -1)
        # if no usable consensus was produced...
        if [ -z "${consensus_src}" ] || [ ! -s "${consensus_src}" ]; then
            echo "ERROR: SHIVER did not produce a consensus FASTA for ${srr}" >&2  # ...report it...
            exit 1                                   # ...and abort
        fi
        # copy to the standard consensus name the harness expects
        cp "${consensus_src}" "${srr}_consensus.fasta"
        sed -i "1s/.*/>${srr}/" "${srr}_consensus.fasta"  # rename the FASTA header to the sample id
    )
}
# export so measure_and_run's child bash can call it
export -f run_shiver

# preferred input: fastp + Kraken2 host-removed reads
KRAKEN2_FASTP_DIR="${REPO_ROOT}/results/download_qc/illumina/kraken2_fastp_out"

# loop over just the Illumina accessions chosen for this comparison
for SRR in $(subset_accessions illumina "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do
    # preferred R1: fastp-filtered + host-removed
    R1="${KRAKEN2_FASTP_DIR}/${SRR}_1.kraken_filtered.fastq.gz"
    # preferred R2: fastp-filtered + host-removed
    R2="${KRAKEN2_FASTP_DIR}/${SRR}_2.kraken_filtered.fastq.gz"
    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ]; then    # if the Kraken2-cleaned reads are missing...
        # ...warn...
        echo "NOTE: no Kraken2-filtered reads for ${SRR}, falling back to plain fastp-trimmed reads (no contamination filtering)." >&2
        R1="${FASTP_DIR}/${SRR}_1.trimmed.fastq.gz"  # ...fall back to plain fastp-trimmed R1...
        R2="${FASTP_DIR}/${SRR}_2.trimmed.fastq.gz"  # ...and R2
    fi
    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ]; then    # if the fastp reads are also missing...
        # ...warn...
        echo "NOTE: no fastp-trimmed reads for ${SRR}, run scripts/download_qc/illumina/download_qc_illumina.sh first for the intended input. Falling back to raw reads." >&2
        R1="${RAW_DIR}/${SRR}_1.fastq.gz"            # ...fall back to raw R1...
        R2="${RAW_DIR}/${SRR}_2.fastq.gz"            # ...and R2
    fi
    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ]; then    # if even raw reads are absent...
        echo "WARNING: no reads at all for ${SRR}, skipping." >&2  # ...warn...
        continue                                     # ...and skip this sample entirely
    fi

    # run each assembler (override the default set via ASSEMBLY_TOOLS)
    for TOOL in ${ASSEMBLY_TOOLS:-bwa_consensus spades shiver}; do
        OUTDIR="${RESULTS_DIR}/${TOOL}_out"          # per-assembler output dir
        # file where measure_and_run records wallclock/RSS
        TIMELOG="${RESULTS_DIR}/${TOOL}_${SRR}.time"
        LOG="${RESULTS_DIR}/${TOOL}_${SRR}.log"      # captured stdout+stderr of the assembler
        # the consensus FASTA each assembler is expected to produce
        CONSENSUS="${OUTDIR}/${SRR}_consensus.fasta"

        echo "=== ${TOOL} on ${SRR} ==="             # progress marker in the log
        # dispatch to the matching exported assembler function
        case "${TOOL}" in
            bwa_consensus)
                # time+run the bwa consensus function
                measure_and_run "${TIMELOG}" -- \
                    bash -c 'run_bwa_consensus "$@"' _ "${SRR}" "${R1}" "${R2}" "${REF_FASTA}" "${OUTDIR}" > "${LOG}" 2>&1 ;;
            spades)
                # time+run the SPAdes function
                measure_and_run "${TIMELOG}" -- \
                    bash -c 'run_spades "$@"' _ "${SRR}" "${R1}" "${R2}" "${REF_FASTA}" "${OUTDIR}" > "${LOG}" 2>&1 ;;
            shiver)
                # time+run the SHIVER function
                measure_and_run "${TIMELOG}" -- \
                    bash -c 'run_shiver "$@"' _ "${SRR}" "${R1}" "${R2}" "${OUTDIR}" "${SHIVER_SETUP_DIR}" "${SHIVER_BIN}" "${SHIVER_DATA_DIR}" > "${LOG}" 2>&1 ;;
        esac
        # capture the assembler's exit status before $? is overwritten
        EXIT_CODE=$?
        # set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file
        parse_time_metrics "${TIMELOG}"

        VALID=0                                      # assume invalid until proven otherwise
        METRIC="n/a"                                 # human-readable key metric, filled in below
        if [ -s "${CONSENSUS}" ]; then               # only evaluate if a consensus was produced
            # consensus length in bp (col 5 of seqkit stats)
            LEN=$(seqkit stats -T "${CONSENSUS}" 2>/dev/null | tail -1 | cut -f5)
            # percent of N (ambiguous) bases
            N_PCT=$(seqkit fx2tab -n -g -B N "${CONSENSUS}" 2>/dev/null | awk -F'\t' '{print $NF}' | tail -1)
            # valid if length is near a full HIV genome (~9kb)
            if [ -n "${LEN}" ] && [ "${LEN}" -ge 8000 ] && [ "${LEN}" -le 10000 ]; then
                VALID=1                              # mark valid
            fi
            METRIC="${LEN}bp, ${N_PCT:-?}% N"        # record length and N% as the key metric
        fi

        # write this assembler's row to summary.tsv
        append_summary_row "assembly_illumina" "${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"
    done
done

# final confirmation pointing the user at the results table
echo "Done. See ${SUMMARY_TSV}"
