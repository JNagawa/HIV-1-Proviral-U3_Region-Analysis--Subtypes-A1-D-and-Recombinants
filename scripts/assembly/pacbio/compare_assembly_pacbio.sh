#!/bin/bash
# Assembly step of the PacBio (HIV-SMRTcap) harness: reference-guided
# minimap2->HXB2 + bcftools consensus  vs  de novo hifiasm, on the
# per-sample proviral-core reads from the extraction step. This mirrors the
# Illumina assembly comparison (bwa_consensus vs SPAdes): a reference-guided
# consensus that is robust at low coverage, against a reference-free de novo
# assembler that can capture divergent/structural variation the reference
# would mask. The tools review names minimap2 as the long-read mapper but no
# long-read de novo assembler, so hifiasm (the PacBio-HiFi standard) is the
# de novo counterpart chosen here.
# Usage: ./compare_assembly_pacbio.sh   (via sbatch scripts/utils/run_comparison_step.slurm.sh)
# -u errors on unset vars, pipefail fails a pipe if any stage fails; no -e so one assembler failing
# doesn't kill the comparison
set -uo pipefail

# absolute path of this script's own dir, so paths work regardless of launch location
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
# repo top-level (three dirs up), the base for every other path below
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
# load shared helpers: measure_and_run, parse_time_metrics, append_summary_row, subset_accessions
source "${REPO_ROOT}/scripts/common/lib_compare.sh"

# HXB2 reference genome used for mapping/consensus
REF_FASTA="${REPO_ROOT}/data/reference/K03455.1.fasta"
# input: per-sample proviral-core reads from the extraction step
PROVIRUS_DIR="${REPO_ROOT}/results/proviral_extraction/pacbio"
# output: all assembly results + summary go here
RESULTS_DIR="${REPO_ROOT}/results/assembly/pacbio"
# the single TSV every assembler appends a timing/validity row to
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
# create the results dir (and parents) if it doesn't exist
mkdir -p "${RESULTS_DIR}"
# seed the manual notes file from the template only on the first run
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

# abort early if the reference is missing
if [ ! -s "${REF_FASTA}" ]; then echo "ERROR: HXB2 reference ${REF_FASTA} not found." >&2; exit 1; fi
# thread count; honour an externally-set THREADS, otherwise default to 4
THREADS="${THREADS:-4}"
# export so the exported assembler functions inherit it
export THREADS

# Reference-guided: map proviral-core reads to HXB2 (HiFi preset) and call a
# consensus. Same shape as the Illumina bwa_consensus function.
# reference-guided assembler (minimap2 + bcftools consensus)
run_minimap2_consensus() {
    # positional args: sample id, reads, reference, output dir
    local srr="$1" reads="$2" ref="$3" outdir="$4"
    mkdir -p "${outdir}"                             # ensure this sample's output dir exists
    # derived output paths
    local bam="${outdir}/${srr}.sorted.bam" vcf="${outdir}/${srr}.vcf.gz" consensus="${outdir}/${srr}_consensus.fasta"
    # map HiFi reads to HXB2 (map-hifi preset), convert to BAM, sort; bail on failure
    minimap2 -a -x map-hifi -t "${THREADS:-4}" "${ref}" "${reads}" 2>"${outdir}/${srr}_minimap2.log" \
        | samtools view -b - | samtools sort -o "${bam}" || exit 1
    samtools index "${bam}" || exit 1                # index the sorted BAM (required by mpileup)
    # call variants haploid (--ploidy 1) since a provirus is a single genome
    bcftools mpileup -Ou -f "${ref}" "${bam}" 2>"${outdir}/${srr}_bcftools.log" \
        | bcftools call -c --ploidy 1 -Oz -o "${vcf}" 2>>"${outdir}/${srr}_bcftools.log" || exit 1
    # index the VCF so bcftools consensus can read it
    tabix -p vcf "${vcf}" || exit 1
    # apply called variants onto the reference to make the consensus
    cat "${ref}" | bcftools consensus "${vcf}" > "${consensus}" 2>>"${outdir}/${srr}_bcftools.log" || exit 1
    sed -i "1s/.*/>${srr}/" "${consensus}"           # rename the FASTA header to the sample id
}
# export so measure_and_run's child bash can call it
export -f run_minimap2_consensus

# De novo: hifiasm assembles the proviral-core reads; take the primary contig
# best matching HXB2 (by BLAST) as the assembly, else the longest -- same
# best-contig selection as the Illumina SPAdes function.
# de novo assembler (hifiasm) with best-contig selection
run_hifiasm() {
    # positional args: sample id, reads, reference, output dir
    local srr="$1" reads="$2" ref="$3" outdir="$4"
    mkdir -p "${outdir}"                             # ensure this sample's output dir exists
    # output-file prefix hifiasm names everything from
    local prefix="${outdir}/${srr}"
    # de novo assemble the HiFi reads; bail on failure
    hifiasm -o "${prefix}" -t "${THREADS:-4}" "${reads}" > "${outdir}/${srr}_hifiasm.log" 2>&1 || exit 1
    # expected primary-contig graph (newer hifiasm naming)
    local gfa="${prefix}.bp.p_ctg.gfa"
    # fall back to the older primary-contig graph name
    [ -s "${gfa}" ] || gfa="${prefix}.p_ctg.gfa"
    # abort if no graph was produced
    [ -s "${gfa}" ] || { echo "ERROR: hifiasm produced no primary-contig GFA for ${srr}" >&2; exit 1; }
    local contigs="${outdir}/${srr}_contigs.fasta"   # FASTA the contigs get extracted into
    # convert GFA segment (S) lines into FASTA records
    awk '/^S/{print ">"$2"\n"$3}' "${gfa}" > "${contigs}"
    # abort if the GFA had no segments
    [ -s "${contigs}" ] || { echo "ERROR: no contigs parsed from ${gfa}" >&2; exit 1; }

    # only do blast selection if the blast tools exist
    if command -v makeblastdb >/dev/null 2>&1 && command -v blastn >/dev/null 2>&1; then
        # build a blast DB from HXB2
        makeblastdb -in "${ref}" -dbtype nucl -out "${outdir}/${srr}_hxb2db" >/dev/null 2>&1
        local best                                   # will hold the id of the best-matching contig
        # blast contigs vs HXB2, sort by bitscore, take the top hit's contig id
        best=$(blastn -query "${contigs}" -db "${outdir}/${srr}_hxb2db" -outfmt "6 qseqid length bitscore" 2>/dev/null \
            | sort -k3,3 -rn | head -1 | cut -f1)
        if [ -n "${best}" ]; then                    # if a best contig was found...
            # extract that one contig as the assembly
            seqkit grep -n -p "${best}" "${contigs}" > "${outdir}/${srr}_consensus.fasta"
            # rename the FASTA header to the sample id
            sed -i "1s/.*/>${srr}/" "${outdir}/${srr}_consensus.fasta"
            return 0                                 # done -- skip the length-based fallback
        fi
    fi
    # fallback: sort contigs by length desc, keep the longest
    seqkit sort -l -r "${contigs}" 2>/dev/null | seqkit head -n 1 > "${outdir}/${srr}_consensus.fasta"
    # rename the FASTA header to the sample id
    sed -i "1s/.*/>${srr}/" "${outdir}/${srr}_consensus.fasta"
}
# export so measure_and_run's child bash can call it
export -f run_hifiasm

# loop over just the PacBio accessions chosen for this comparison
for SRR in $(subset_accessions pacbio "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do
    READS="${PROVIRUS_DIR}/${SRR}.provirus.fasta"    # this sample's extracted proviral-core reads
    if [ ! -s "${READS}" ]; then                     # if the proviral cores are missing...
        # ...warn...
        echo "WARNING: no proviral cores for ${SRR} (run extract_provirus_pacbio.sh first), skipping." >&2
        continue                                     # ...and skip this sample entirely
    fi

    # run each assembler (override the default set via ASSEMBLY_TOOLS)
    for TOOL in ${ASSEMBLY_TOOLS:-minimap2_consensus hifiasm}; do
        # hifiasm is optional...
        if [ "${TOOL}" = "hifiasm" ] && ! command -v hifiasm >/dev/null 2>&1; then
            # ...note its absence...
            echo "NOTE: hifiasm not installed, skipping (see HIV_U3analysis_env.yml)." >&2
            continue                                 # ...and skip it if the binary isn't on PATH
        fi
        OUTDIR="${RESULTS_DIR}/${TOOL}_out"          # per-assembler output dir
        # file where measure_and_run records wallclock/RSS
        TIMELOG="${RESULTS_DIR}/${TOOL}_${SRR}.time"
        LOG="${RESULTS_DIR}/${TOOL}_${SRR}.log"      # captured stdout+stderr of the assembler
        # the consensus FASTA each assembler is expected to produce
        CONSENSUS="${OUTDIR}/${SRR}_consensus.fasta"

        echo "=== ${TOOL} on ${SRR} ==="             # progress marker in the log
        # time+run the matching run_<tool> function (name built from ${TOOL})
        measure_and_run "${TIMELOG}" -- \
            bash -c 'run_'"${TOOL}"' "$@"' _ "${SRR}" "${READS}" "${REF_FASTA}" "${OUTDIR}" > "${LOG}" 2>&1
        # capture the assembler's exit status before $? is overwritten
        EXIT_CODE=$?
        # set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file
        parse_time_metrics "${TIMELOG}"

        VALID=0; METRIC="n/a"                        # assume invalid until proven otherwise
        if [ -s "${CONSENSUS}" ]; then               # only evaluate if a consensus was produced
            # consensus length in bp (col 5 of seqkit stats)
            LEN=$(seqkit stats -T "${CONSENSUS}" 2>/dev/null | tail -1 | cut -f5)
            # percent of N (ambiguous) bases
            N_PCT=$(seqkit fx2tab -n -g -B N "${CONSENSUS}" 2>/dev/null | awk -F'\t' '{print $NF}' | tail -1)
            # valid if length is near a full HIV genome (~9kb)
            if [ -n "${LEN}" ] && [ "${LEN}" -ge 8000 ] && [ "${LEN}" -le 10000 ]; then VALID=1; fi
            METRIC="${LEN}bp, ${N_PCT:-?}% N"        # record length and N% as the key metric
        fi
        # write this assembler's row to summary.tsv
        append_summary_row "assembly_pacbio" "${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"
    done
done

# final confirmation pointing the user at the results table
echo "Done. See ${SUMMARY_TSV}"
