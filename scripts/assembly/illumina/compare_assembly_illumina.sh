#!/bin/bash
# Assembly step of the tool-comparison harness: BWA+bcftools-consensus vs
# SPAdes vs SHIVER on the Illumina subset -- one script for all three
# assemblers and their comparison.
# Uses each sample's fastp-trimmed reads from download_qc if available,
# else falls back to raw reads with a warning.
# Usage: ./compare_assembly_illumina.sh
set -uo pipefail                                     # -u errors on unset vars, pipefail fails a pipe if any stage fails; no -e so one assembler failing doesn't kill the comparison

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"           # absolute path of this script's own dir, so paths work regardless of launch location
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"     # repo top-level (three dirs up), the base for every other path below
source "${REPO_ROOT}/scripts/common/lib_compare.sh"  # load shared helpers: measure_and_run, parse_time_metrics, append_summary_row, subset_accessions

RAW_DIR="${REPO_ROOT}/data/raw/illumina"                            # fallback input: raw downloaded Illumina FASTQs
FASTP_DIR="${REPO_ROOT}/results/download_qc/illumina/fastp_out"     # preferred input: fastp-trimmed reads from the QC step
REF_FASTA="${REPO_ROOT}/data/reference/K03455.1.fasta"             # HXB2 reference genome used for mapping/consensus
RESULTS_DIR="${REPO_ROOT}/results/assembly/illumina"               # output: all assembly results + summary go here
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"                           # the single TSV every assembler appends a timing/validity row to
mkdir -p "${RESULTS_DIR}"                                          # create the results dir (and parents) if it doesn't exist
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"  # seed the manual notes file from the template only on the first run

# bwa mem requires the reference pre-indexed (bwa index) -- unlike
# minimap2/SPAdes, it errors out immediately without it, so we index it once here

if [ ! -s "${REF_FASTA}.bwt" ]; then                 # skip indexing if the .bwt index already exists (idempotent reruns)
    bwa index "${REF_FASTA}" > "${RESULTS_DIR}/bwa_index.log" 2>&1  # build the bwa index once, capturing output to a log
fi

THREADS="${THREADS:-4}"                              # thread count; honour an externally-set THREADS, otherwise default to 4
export THREADS                                       # export so the exported assembler functions inherit it

SHIVER_SETUP_DIR="${STAGE_DIR}/shiver_setup"                              # dir holding SHIVER's config.sh and reference alignment
SHIVER_BIN="${REPO_ROOT}/scripts/tools/shiver/bin"                       # dir with the SHIVER executable scripts
SHIVER_DATA_DIR="${REPO_ROOT}/scripts/tools/shiver/data/example_input"   # dir with SHIVER's bundled adapters/primers FASTAs

# --- One function per assembler. Each is `export -f`'d below so it can be
# invoked as `bash -c '<func> "$@"' _ args...` -- measure_and_run wraps its
# command in /usr/bin/time, which (like any external-process wrapper) can
# only time a real executable, not a shell function living in this script's
# own process, so this is how each function still gets its own wall-clock
# and peak-RSS measurement. ---

# Reference-mapping assembly: the approach currently used by the production
# pipeline (illumina_u3analysis.sh Step 7) -- bwa mem to HXB2 + bcftools
# consensus. The baseline SPAdes/SHIVER are compared against.
run_bwa_consensus() {                                # baseline reference-mapping assembler (bwa + bcftools consensus)
    local srr="$1" r1="$2" r2="$3" ref_fasta="$4" outdir="$5"  # positional args: sample id, paired reads, reference, output dir
    mkdir -p "${outdir}"                             # ensure this sample's output dir exists
    local bam="${outdir}/${srr}.sorted.bam" vcf="${outdir}/${srr}.vcf.gz" consensus="${outdir}/${srr}_consensus.fasta"  # derived output paths

    bwa mem -t "${THREADS:-4}" "${ref_fasta}" "${r1}" "${r2}" 2>"${outdir}/${srr}_bwa.log" | \
        samtools view -b - | samtools sort -o "${bam}" || exit 1  # map reads to HXB2, convert to BAM, sort; bail if any stage fails
    samtools index "${bam}" || exit 1                # index the sorted BAM (required by mpileup)

    bcftools mpileup -Ou -f "${ref_fasta}" "${bam}" 2>"${outdir}/${srr}_bcftools.log" | \
        bcftools call -c -Oz -o "${vcf}" 2>>"${outdir}/${srr}_bcftools.log" || exit 1  # pile up bases and call variants against the reference
    tabix -p vcf "${vcf}" || exit 1                  # index the VCF so bcftools consensus can read it
    cat "${ref_fasta}" | bcftools consensus "${vcf}" > "${consensus}" 2>>"${outdir}/${srr}_bcftools.log" || exit 1  # apply the called variants onto the reference to make the consensus
    sed -i "1s/.*/>${srr}/" "${consensus}"           # rename the FASTA header to the sample id
}
export -f run_bwa_consensus                          # export so measure_and_run's child bash can call it

# De novo assembly with SPAdes, then pick the single best-matching contig
# against HXB2 (by blast) as the "assembly" output for comparison.
run_spades() {                                       # de novo assembler (SPAdes) with best-contig selection
    local srr="$1" r1="$2" r2="$3" ref_fasta="$4" outdir="$5"  # positional args: sample id, paired reads, reference, output dir
    local spades_dir="${outdir}/${srr}_spades"       # per-sample SPAdes working dir
    mkdir -p "${outdir}"                             # ensure output dir exists

    spades.py --careful -1 "${r1}" -2 "${r2}" -o "${spades_dir}" -t "${THREADS:-4}" \
        > "${outdir}/${srr}_spades.log" 2>&1         # assemble reads into contigs; --careful reduces mismatches/indels
    local contigs="${spades_dir}/contigs.fasta"      # SPAdes' contig output file
    [ -s "${contigs}" ] || { echo "ERROR: SPAdes produced no contigs for ${srr}" >&2; exit 1; }  # fail if no contigs were produced

    # Pick the contig with the best blast hit against HXB2 as the
    # near-full-genome candidate (SPAdes commonly returns many
    # short/host-contaminant contigs).
    if command -v makeblastdb >/dev/null 2>&1 && command -v blastn >/dev/null 2>&1; then  # only do blast selection if the blast tools exist
        makeblastdb -in "${ref_fasta}" -dbtype nucl -out "${outdir}/${srr}_hxb2db" >/dev/null 2>&1  # build a blast DB from HXB2
        local best_id                                # will hold the id of the best-matching contig
        best_id=$(blastn -query "${contigs}" -db "${outdir}/${srr}_hxb2db" -outfmt "6 qseqid length bitscore" 2>/dev/null \
            | sort -k3,3 -rn | head -1 | cut -f1)    # blast contigs vs HXB2, sort by bitscore, take the top hit's contig id
        if [ -n "${best_id}" ]; then                 # if a best contig was found...
            seqkit grep -n -p "${best_id}" "${contigs}" > "${outdir}/${srr}_consensus.fasta"  # extract that one contig as the assembly
            sed -i "1s/.*/>${srr}/" "${outdir}/${srr}_consensus.fasta"  # rename the FASTA header to the sample id
            return 0                                 # done -- skip the length-based fallback
        fi
    fi

    # Fallback (no blast available): just take the longest contig.
    seqkit sort -l -r "${contigs}" 2>/dev/null | seqkit head -n 1 > "${outdir}/${srr}_consensus.fasta"  # sort contigs by length desc, keep the longest
    sed -i "1s/.*/>${srr}/" "${outdir}/${srr}_consensus.fasta"  # rename the FASTA header to the sample id
}
export -f run_spades                                 # export so measure_and_run's child bash can call it


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

run_shiver() {                                       # SHIVER 3-step assembler with per-sample custom reference
    local srr="$1" r1="$2" r2="$3" outdir="$4" setup_dir="$5" shiver_bin="$6" shiver_data_dir="$7"  # positional args passed in from the loop
    local ref_alignment="${setup_dir}/HIV1_COM_ref_alignment.fasta"  # LANL multi-subtype reference alignment panel
    if [ ! -s "${ref_alignment}" ]; then             # the panel is a manual one-time download, so verify it exists
        echo "ERROR: ${ref_alignment} not found. See shiver_setup/SOURCE.md for the one-time manual LANL download step this requires." >&2  # tell the user how to obtain it
        exit 1                                       # cannot proceed without the reference panel
    fi

    local adapters="${shiver_data_dir}/adapters_Illumina.fasta"      # adapter sequences SHIVER trims from reads
    local primers="${shiver_data_dir}/primers_GallEtAl2012.fasta"    # PCR primer sequences SHIVER trims from reads
    local init_dir="${outdir}/shiver_init"           # dir for SHIVER's one-time init output
    mkdir -p "${outdir}"                             # ensure this sample's output dir exists

    # Resolve R1/R2 to absolute paths before the subshell cd's, since
    # callers may pass relative paths.
    r1="$(cd "$(dirname "${r1}")" && pwd)/$(basename "${r1}")"  # canonicalize R1 to an absolute path (survives the later cd)
    r2="$(cd "$(dirname "${r2}")" && pwd)/$(basename "${r2}")"  # canonicalize R2 to an absolute path (survives the later cd)

    # shiver_init.sh only needs to run once (its OutDir must not pre-exist)
    if [ ! -d "${init_dir}" ]; then                  # only initialise if it hasn't been done yet (dir must not pre-exist)
        "${shiver_bin}/shiver_init.sh" "${init_dir}" "${setup_dir}/config.sh" \
            "${ref_alignment}" "${adapters}" "${primers}" \
            > "${outdir}/shiver_init.log" 2>&1        # build SHIVER's initialisation dir from the reference/adapters/primers
        if [ $? -ne 0 ]; then                        # if init failed...
            echo "ERROR: shiver_init.sh failed, see ${outdir}/shiver_init.log" >&2  # ...point at the log...
            exit 1                                   # ...and abort
        fi
    fi

    (                                                # subshell so the cd below is local and doesn't move the main script
        cd "${outdir}" || exit 1                     # SHIVER writes to cwd, so enter this sample's output dir

        # Step 1 of the chain: assemble contigs with SPAdes (cost logged
        # separately so it isn't double-counted against SHIVER's own
        # runtime/RSS when SPAdes is ALSO being benchmarked as its own
        # competing assembler in this same step).
        local contigs_dir="${srr}_spades_for_shiver" # separate SPAdes dir so its cost is logged apart from SHIVER's
        spades.py --careful -1 "${r1}" -2 "${r2}" -o "${contigs_dir}" -t "${THREADS:-4}" \
            > "${srr}_spades_for_shiver.log" 2>&1     # assemble contigs to feed into SHIVER
        local contigs="${contigs_dir}/contigs.fasta" # the contigs SPAdes produced
        [ -s "${contigs}" ] || { echo "ERROR: SPAdes (for SHIVER) produced no contigs for ${srr}" >&2; exit 1; }  # abort if no contigs

        # Step 2: align contigs to the reference set, excluding contamination.
        "${shiver_bin}/shiver_align_contigs.sh" "${init_dir}" "${setup_dir}/config.sh" \
            "${contigs}" "${srr}" > "${srr}_align_contigs.log" 2>&1  # clean/align contigs to the reference panel
        if [ $? -ne 0 ]; then                        # if alignment failed...
            echo "ERROR: shiver_align_contigs.sh failed, see ${outdir}/${srr}_align_contigs.log" >&2  # ...point at the log...
            exit 1                                   # ...and abort
        fi
        local blast_file="${srr}.blast" cut_wrefs="${srr}_cut_wRefs.fasta"  # outputs from the align step needed by the map step
        [ -s "${blast_file}" ] || { echo "ERROR: no blast hits for ${srr} -- contigs may not be HIV" >&2; exit 1; }  # no hits means contigs aren't HIV

        # Step 3: map reads using the aligned contigs to build the final assembly.
        "${shiver_bin}/shiver_map_reads.sh" "${init_dir}" "${setup_dir}/config.sh" \
            "${contigs}" "${srr}" "${blast_file}" "${cut_wrefs}" "${r1}" "${r2}" \
            > "${srr}_map_reads.log" 2>&1             # map reads to the sample-specific reference to build the consensus
        if [ $? -ne 0 ]; then                        # if mapping failed...
            echo "ERROR: shiver_map_reads.sh failed, see ${outdir}/${srr}_map_reads.log" >&2  # ...point at the log...
            exit 1                                   # ...and abort
        fi

        # SHIVER names its consensus output "<SID>_remap_consensus_MinCov_*.fasta";
        # take the first (default coverage threshold) match.
        local consensus_src                          # will hold the SHIVER consensus filename
        consensus_src=$(compgen -G "${srr}_remap_consensus_MinCov_*.fasta" | head -1)  # glob for SHIVER's consensus output, take the first match
        if [ -z "${consensus_src}" ] || [ ! -s "${consensus_src}" ]; then  # if no usable consensus was produced...
            echo "ERROR: SHIVER did not produce a consensus FASTA for ${srr}" >&2  # ...report it...
            exit 1                                   # ...and abort
        fi
        cp "${consensus_src}" "${srr}_consensus.fasta"  # copy to the standard consensus name the harness expects
        sed -i "1s/.*/>${srr}/" "${srr}_consensus.fasta"  # rename the FASTA header to the sample id
    )
}
export -f run_shiver                                 # export so measure_and_run's child bash can call it

KRAKEN2_FASTP_DIR="${REPO_ROOT}/results/download_qc/illumina/kraken2_fastp_out"  # preferred input: fastp + Kraken2 host-removed reads

for SRR in $(subset_accessions illumina "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do  # loop over just the Illumina accessions chosen for this comparison
    R1="${KRAKEN2_FASTP_DIR}/${SRR}_1.kraken_filtered.fastq.gz"  # preferred R1: fastp-filtered + host-removed
    R2="${KRAKEN2_FASTP_DIR}/${SRR}_2.kraken_filtered.fastq.gz"  # preferred R2: fastp-filtered + host-removed
    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ]; then    # if the Kraken2-cleaned reads are missing...
        echo "NOTE: no Kraken2-filtered reads for ${SRR}, falling back to plain fastp-trimmed reads (no contamination filtering)." >&2  # ...warn...
        R1="${FASTP_DIR}/${SRR}_1.trimmed.fastq.gz"  # ...fall back to plain fastp-trimmed R1...
        R2="${FASTP_DIR}/${SRR}_2.trimmed.fastq.gz"  # ...and R2
    fi
    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ]; then    # if the fastp reads are also missing...
        echo "NOTE: no fastp-trimmed reads for ${SRR}, run scripts/download_qc/illumina/download_qc_illumina.sh first for the intended input. Falling back to raw reads." >&2  # ...warn...
        R1="${RAW_DIR}/${SRR}_1.fastq.gz"            # ...fall back to raw R1...
        R2="${RAW_DIR}/${SRR}_2.fastq.gz"            # ...and R2
    fi
    if [ ! -s "${R1}" ] || [ ! -s "${R2}" ]; then    # if even raw reads are absent...
        echo "WARNING: no reads at all for ${SRR}, skipping." >&2  # ...warn...
        continue                                     # ...and skip this sample entirely
    fi

    for TOOL in ${ASSEMBLY_TOOLS:-bwa_consensus spades shiver}; do  # run each assembler (override the default set via ASSEMBLY_TOOLS)
        OUTDIR="${RESULTS_DIR}/${TOOL}_out"          # per-assembler output dir
        TIMELOG="${RESULTS_DIR}/${TOOL}_${SRR}.time" # file where measure_and_run records wallclock/RSS
        LOG="${RESULTS_DIR}/${TOOL}_${SRR}.log"      # captured stdout+stderr of the assembler
        CONSENSUS="${OUTDIR}/${SRR}_consensus.fasta" # the consensus FASTA each assembler is expected to produce

        echo "=== ${TOOL} on ${SRR} ==="             # progress marker in the log
        case "${TOOL}" in                            # dispatch to the matching exported assembler function
            bwa_consensus)
                measure_and_run "${TIMELOG}" -- \
                    bash -c 'run_bwa_consensus "$@"' _ "${SRR}" "${R1}" "${R2}" "${REF_FASTA}" "${OUTDIR}" > "${LOG}" 2>&1 ;;  # time+run the bwa consensus function
            spades)
                measure_and_run "${TIMELOG}" -- \
                    bash -c 'run_spades "$@"' _ "${SRR}" "${R1}" "${R2}" "${REF_FASTA}" "${OUTDIR}" > "${LOG}" 2>&1 ;;  # time+run the SPAdes function
            shiver)
                measure_and_run "${TIMELOG}" -- \
                    bash -c 'run_shiver "$@"' _ "${SRR}" "${R1}" "${R2}" "${OUTDIR}" "${SHIVER_SETUP_DIR}" "${SHIVER_BIN}" "${SHIVER_DATA_DIR}" > "${LOG}" 2>&1 ;;  # time+run the SHIVER function
        esac
        EXIT_CODE=$?                                 # capture the assembler's exit status before $? is overwritten
        parse_time_metrics "${TIMELOG}"              # set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file

        VALID=0                                      # assume invalid until proven otherwise
        METRIC="n/a"                                 # human-readable key metric, filled in below
        if [ -s "${CONSENSUS}" ]; then               # only evaluate if a consensus was produced
            LEN=$(seqkit stats -T "${CONSENSUS}" 2>/dev/null | tail -1 | cut -f5)  # consensus length in bp (col 5 of seqkit stats)
            N_PCT=$(seqkit fx2tab -n -g -B N "${CONSENSUS}" 2>/dev/null | awk -F'\t' '{print $NF}' | tail -1)  # percent of N (ambiguous) bases
            if [ -n "${LEN}" ] && [ "${LEN}" -ge 8000 ] && [ "${LEN}" -le 10000 ]; then  # valid if length is near a full HIV genome (~9kb)
                VALID=1                              # mark valid
            fi
            METRIC="${LEN}bp, ${N_PCT:-?}% N"        # record length and N% as the key metric
        fi

        append_summary_row "assembly_illumina" "${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"  # write this assembler's row to summary.tsv
    done
done

echo "Done. See ${SUMMARY_TSV}"                      # final confirmation pointing the user at the results table
