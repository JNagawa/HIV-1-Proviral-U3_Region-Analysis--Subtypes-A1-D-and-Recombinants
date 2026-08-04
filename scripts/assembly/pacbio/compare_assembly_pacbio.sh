#!/bin/bash
# Assembly step of the PacBio (HIV-SMRTcap) harness, on the per-sample
# proviral-core reads from the extraction step. Three arms:
#   minimap2_consensus  reference-guided against HXB2 (subtype B), the fixed
#                       baseline and the coordinate system everything
#                       downstream is used to
#   minimap2_bestref    the identical pipeline against the closest match from
#                       LANL's Group M compendium, so the cost of forcing
#                       subtype-A1/D Rakai reads onto a subtype-B backbone is
#                       measurable rather than assumed
#   hifiasm             de novo, for the reference-free comparison
# This mirrors the Illumina assembly comparison (bwa_consensus vs SPAdes): a
# reference-guided consensus that is robust at low coverage, against a
# reference-free de novo assembler that can capture divergent/structural
# variation the reference would mask. The tools review names minimap2 as the
# long-read mapper but no long-read de novo assembler, so hifiasm (the
# PacBio-HiFi standard) is the de novo counterpart chosen here.
#
# Note on the de novo arm: these libraries are amplicon, not shotgun -- reads
# stack into 1-3 discrete intervals instead of tiling the genome. OLC assembly
# needs dovetail overlaps, so hifiasm can only ever assemble the samples that
# have two overlapping amplicons (211_0 of the current four). That is recorded
# as a result, not treated as a tool failure. See run_hifiasm below.
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

# hifiasm's defaults assume a gigabase-scale diploid genome and produce an empty
# overlap graph (0 contigs) on these proviral read sets, which are 19-75 reads
# covering a 9.7kb genome. The overrides match the actual target:
#   -f0        disable the bloom-filter k-mer pre-filter -- it is sized for large
#              genomes and discards the low-count k-mers that are all we have
#   -l0        no purge-dups; there are no duplicate haplotigs to purge here
#   --n-hap 1  a provirus is a single haplotype, not the default diploid pair
#   --hg-size  pin the genome size to HIV's ~9.7kb so hifiasm stops inferring
#              coverage as if this were a gigabase genome
#   -n 0       do NOT drop tip unitigs of <=3 reads (the default -n 3). On a
#              30-read library almost every unitig is a "tip" by that measure,
#              so the default deletes the whole graph. This is the single
#              override that lets 211_0 assemble at all.
#   --ctg-n 0  same reasoning, at the contig rather than unitig stage.
# Measured 2026-08-04 across a 7-config x 4-sample sweep: --hom-cov tuning
# (1..8), read dedup and contained-read removal all change nothing, and -a 0
# actively breaks the one sample that works, so none of them are used here.
HIFIASM_OPTS="${HIFIASM_OPTS:--f0 -l0 --n-hap 1 --hg-size 10k -n 0 --ctg-n 0}"
# export so the exported run_hifiasm inherits it
export HIFIASM_OPTS

# These libraries are subtype A1/D (Rakai cohort) but HXB2 is subtype B, so the
# reads sit at only ~90% identity to it. Every position bcftools does not call
# is left as the reference base, which means a subtype-B sequence is silently
# substituted for un-called sample sequence. Picking the closest reference from
# LANL's Group M compendium (180 seqs, all subtypes + CRFs -- already vendored
# for the Illumina SHIVER arm) shrinks that gap before masking closes the rest.
REF_PANEL_ALN="${REPO_ROOT}/scripts/assembly/illumina/shiver_setup/HIV1_COM_ref_alignment.fasta"
# the panel as ungapped mapping targets, derived from the alignment above
REF_PANEL="${RESULTS_DIR}/refpanel/panel.fasta"
# HXB2 coordinates of the two U3 regions a candidate reference has to contain
U3_5P_START=1 ; U3_5P_END=455
U3_3P_START=9086 ; U3_3P_END=9550
if [ -s "${REF_PANEL_ALN}" ]; then
    mkdir -p "${RESULTS_DIR}/refpanel"
    # rebuild only when missing or older than the alignment it derives from
    if [ ! -s "${REF_PANEL}" ] || [ "${REF_PANEL_ALN}" -nt "${REF_PANEL}" ]; then
        echo "Building reference panel from ${REF_PANEL_ALN}..."
        # LANL compendium genomes conventionally run gag->nef and omit the LTRs:
        # only 22 of the 180 carry both. Mapping to an LTR-less reference throws
        # away exactly the reads this project exists to analyse -- measured
        # 2026-08-04, the unfiltered panel put 124_4 on a reference spanning HXB2
        # 649-9594, which silently discarded its 24 3'LTR reads (40 mapped -> 16).
        # So a candidate only qualifies if it actually carries both U3 regions.
        # Presence is read off the alignment directly: HXB2's ungapped positions
        # give the column range for each U3, and a candidate must be non-gap
        # across at least 90% of both.
        seqkit seq -w 0 "${REF_PANEL_ALN}" \
            | awk -v u5s="${U3_5P_START}" -v u5e="${U3_5P_END}" \
                  -v u3s="${U3_3P_START}" -v u3e="${U3_3P_END}" -v minf=0.9 '
                /^>/ {id=substr($0,2); order[++n]=id; next}
                {seq[id]=$0}
                END{
                    for (i=1;i<=n;i++) if (order[i] ~ /K03455/) hx=seq[order[i]]
                    if (hx=="") exit 1
                    p=0
                    for (c=1;c<=length(hx);c++) if (substr(hx,c,1)!="-") col[++p]=c
                    for (i=1;i<=n;i++) {
                        id=order[i]; s=seq[id]
                        g5=0; for (c=col[u5s]; c<=col[u5e]; c++) if (substr(s,c,1)!="-") g5++
                        g3=0; for (c=col[u3s]; c<=col[u3e]; c++) if (substr(s,c,1)!="-") g3++
                        if (g5/(u5e-u5s+1) >= minf && g3/(u3e-u3s+1) >= minf) print id
                    }
                }' > "${RESULTS_DIR}/refpanel/ltr_complete.ids"
        seqkit grep -f "${RESULTS_DIR}/refpanel/ltr_complete.ids" "${REF_PANEL_ALN}" 2>/dev/null \
            | seqkit seq -g -w 0 > "${REF_PANEL}"
        echo "  panel: $(grep -c '^>' "${REF_PANEL}" 2>/dev/null || echo 0) LTR-complete of $(grep -c '^>' "${REF_PANEL_ALN}") references"
    fi
    # an empty filter result would silently send every sample to the fallback
    if [ ! -s "${REF_PANEL}" ]; then
        echo "NOTE: no LTR-complete references in the panel; falling back to HXB2." >&2
        REF_PANEL=""
    fi
else
    # no panel vendored -- every sample falls back to HXB2, as before
    echo "NOTE: ${REF_PANEL_ALN} not found; falling back to HXB2 for all samples." >&2
    REF_PANEL=""
fi
# export so the exported assembler functions inherit it
export REF_PANEL

# Picks the closest reference for one sample and prints its FASTA path.
# Maps the sample's reads against the whole panel and scores each candidate by
# total matching bases (PAF column 10) summed over all its alignments, which
# rewards both identity and how much of the read set a candidate explains.
# Progress goes to stderr so stdout stays a clean single path for the caller.
select_reference() {
    # positional args: sample id, reads, fallback reference, output dir
    local srr="$1" reads="$2" default_ref="$3" outdir="$4"
    # this sample's chosen reference, copied out of the panel
    local chosen="${outdir}/${srr}_reference.fasta"
    # human-readable record of which reference won, for the summary table
    local report="${outdir}/${srr}_reference.txt"
    # Whichever branch below wins, the chosen FASTA is rewritten in place. Its
    # .fai must be dropped with it: samtools/bcftools trust an existing index
    # over the file, so a stale .fai left by a previous run (different sample,
    # different panel) makes every later mpileup fail with "sequence not found"
    # and silently emit an empty VCF -- a consensus with no variants applied.
    finalise_reference() {
        rm -f "${chosen}.fai"
        samtools faidx "${chosen}" 2>/dev/null
        printf '%s\n' "${chosen}"
    }
    # with no panel available there is nothing to choose between
    if [ -z "${REF_PANEL:-}" ] || [ ! -s "${REF_PANEL}" ]; then
        cp "${default_ref}" "${chosen}"
        printf 'HXB2 (no panel)\n' > "${report}"
        finalise_reference
        return 0
    fi
    local best
    # --secondary=no keeps one alignment per read so the sum is not inflated
    best=$(minimap2 -x map-hifi -t "${THREADS:-4}" --secondary=no "${REF_PANEL}" "${reads}" 2>/dev/null \
        | awk -F'\t' '{m[$6]+=$10} END{for (r in m) print m[r]"\t"r}' \
        | sort -k1,1 -rn | head -1 | cut -f2)
    # no read mapped anywhere in the panel -- fall back rather than fail
    if [ -z "${best}" ]; then
        cp "${default_ref}" "${chosen}"
        printf 'HXB2 (no panel hit)\n' > "${report}"
        finalise_reference
        return 0
    fi
    seqkit grep -n -p "${best}" "${REF_PANEL}" > "${chosen}"
    # panel IDs are dot-delimited as <subtype>.<country>.<year>.<name>.<accession>
    printf '%s\n' "${best%%.*}:${best##*.}" > "${report}"
    echo "  ${srr}: closest reference ${best}" >&2
    finalise_reference
}
# export so measure_and_run's child bash can call it
export -f select_reference

# Reference-guided: map proviral-core reads (HiFi preset) to whichever reference
# the caller chose and call a consensus. Same shape as the Illumina
# bwa_consensus function. Split out from the two wrappers below so the HXB2 and
# subtype-matched arms differ ONLY in their reference -- everything else about
# the two runs is identical, which is what makes them comparable.
consensus_against_ref() {
    # positional args: sample id, reads, reference to map against, output dir
    local srr="$1" reads="$2" ref="$3" outdir="$4"
    # derived output paths
    local bam="${outdir}/${srr}.sorted.bam" vcf="${outdir}/${srr}.vcf.gz" consensus="${outdir}/${srr}_consensus.fasta"
    # map reads to the reference (map-hifi preset), convert to BAM, sort; bail on failure
    minimap2 -a -x map-hifi -t "${THREADS:-4}" "${ref}" "${reads}" 2>"${outdir}/${srr}_minimap2.log" \
        | samtools view -b - | samtools sort -o "${bam}" || exit 1
    samtools index "${bam}" || exit 1                # index the sorted BAM (required by mpileup)
    # call variants haploid (--ploidy 1) since a provirus is a single genome
    bcftools mpileup -Ou -f "${ref}" "${bam}" 2>"${outdir}/${srr}_bcftools.log" \
        | bcftools call -c --ploidy 1 -Oz -o "${vcf}" 2>>"${outdir}/${srr}_bcftools.log" || exit 1
    # index the VCF so bcftools consensus can read it
    tabix -p vcf "${vcf}" || exit 1
    # A VCF with no records while reads did map means mpileup rejected the
    # reference (usually a name mismatch against the BAM) and the consensus
    # would silently come out as the bare reference. Fail loudly instead.
    if [ "$(bcftools view -H "${vcf}" 2>/dev/null | wc -l)" -eq 0 ] \
       && [ "$(samtools view -c -F 0x904 "${bam}" 2>/dev/null)" -gt 0 ]; then
        echo "ERROR: ${srr} has mapped reads but an empty VCF -- see ${outdir}/${srr}_bcftools.log" >&2
        exit 1
    fi

    # bcftools consensus starts FROM the reference and only edits positions it
    # called a variant at, so any position with no read coverage silently keeps
    # the reference base -- a subtype-B sequence presented as sample data, with
    # no N to flag it. Measured 2026-08-04: that left 203_3's consensus 88%
    # HXB2, 128_5's 44%, 211_0's 39%, every one of them reported as "0.00% N".
    # Build a BED of the zero-coverage stretches so -m can mask them instead.
    local mask="${outdir}/${srr}.uncovered.bed"
    local refid reflen
    # the reference's sequence name, which must match the BED/VCF CHROM field
    refid=$(seqkit fx2tab -n -i "${ref}" | head -1 | cut -f1)
    reflen=$(seqkit fx2tab -n -l -i "${ref}" | head -1 | cut -f2)
    # Positions with at least one read. samtools depth -a is asked for all
    # positions, but emits nothing at all when a BAM has no alignments, so the
    # uncovered intervals are derived by complementing this list against the
    # reference length rather than by reading zero-depth rows back.
    samtools depth -a "${bam}" 2>/dev/null | awk -F'\t' '$3>0{print $2}' > "${outdir}/${srr}.covered.pos"
    awk -v L="${reflen}" -v CHR="${refid}" -v OFS='\t' '
        {cov[$1]=1}
        END{
            inrun=0
            for (i=1; i<=L; i++) {
                if (!(i in cov)) { if (!inrun) {s=i; inrun=1} }
                else if (inrun) { print CHR, s-1, i-1; inrun=0 }
            }
            if (inrun) print CHR, s-1, L
        }' "${outdir}/${srr}.covered.pos" > "${mask}"

    # apply called variants onto the reference, N-masking everything unsequenced
    bcftools consensus -m "${mask}" -f "${ref}" "${vcf}" \
        > "${consensus}" 2>>"${outdir}/${srr}_bcftools.log" || exit 1
    sed -i "1s/.*/>${srr}/" "${consensus}"           # rename the FASTA header to the sample id
}
# export so the two wrappers below inherit it
export -f consensus_against_ref

# Arm 1: always HXB2, the subtype-B standard coordinate reference. Kept as-is so
# there is a fixed baseline to compare the subtype-matched arm against.
run_minimap2_consensus() {
    # positional args: sample id, reads, reference, output dir
    local srr="$1" reads="$2" ref="$3" outdir="$4"
    mkdir -p "${outdir}"                             # ensure this sample's output dir exists
    # a stale consensus from an earlier run would otherwise be scored as this run's
    rm -f "${outdir}/${srr}_consensus.fasta" "${outdir}/${srr}_metric.txt"
    # record the reference for the summary column, uniform with the bestref arm
    printf 'HXB2\n' > "${outdir}/${srr}_reference.txt"
    consensus_against_ref "${srr}" "${reads}" "${ref}" "${outdir}"
}
# export so measure_and_run's child bash can call it
export -f run_minimap2_consensus

# Arm 2: same pipeline, but mapped to the closest Group M reference for this
# sample instead of HXB2. Isolates how much of the consensus is an artefact of
# forcing subtype-A1/D reads onto a subtype-B backbone.
run_minimap2_bestref() {
    # positional args: sample id, reads, fallback reference, output dir
    local srr="$1" reads="$2" ref="$3" outdir="$4"
    mkdir -p "${outdir}"                             # ensure this sample's output dir exists
    # a stale consensus from an earlier run would otherwise be scored as this run's
    rm -f "${outdir}/${srr}_consensus.fasta" "${outdir}/${srr}_metric.txt"
    local chosen
    # pick this sample's closest panel reference (falls back to HXB2 internally)
    chosen=$(select_reference "${srr}" "${reads}" "${ref}" "${outdir}") || exit 1
    consensus_against_ref "${srr}" "${reads}" "${chosen}" "${outdir}"
}
# export so measure_and_run's child bash can call it
export -f run_minimap2_bestref

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
    # hifiasm caches error-corrected reads and overlaps in <prefix>.*.bin and
    # silently reuses them whenever the prefix matches -- even if the input reads
    # have changed since. That made the Jul-29 re-run assemble Jul-22 overlaps and
    # emit zero contigs. Drop the cache so a re-run always uses the current reads.
    rm -f "${prefix}.ec.bin" "${prefix}.ovlp.source.bin" "${prefix}.ovlp.reverse.bin"
    # a stale consensus from an earlier run would otherwise be scored as this run's
    rm -f "${outdir}/${srr}_consensus.fasta" "${outdir}/${srr}_metric.txt"
    # de novo assemble the HiFi reads; bail on failure
    # HIFIASM_OPTS is deliberately unquoted so it word-splits into separate flags
    hifiasm ${HIFIASM_OPTS} -o "${prefix}" -t "${THREADS:-4}" "${reads}" \
        > "${outdir}/${srr}_hifiasm.log" 2>&1 || exit 1
    # expected primary-contig graph (newer hifiasm naming)
    local gfa="${prefix}.bp.p_ctg.gfa"
    # fall back to the older primary-contig graph name
    [ -s "${gfa}" ] || gfa="${prefix}.p_ctg.gfa"
    # An empty graph here is a RESULT, not a crash. These are amplicon libraries:
    # reads stack into 1-3 discrete intervals rather than tiling the genome, and
    # OLC assembly needs dovetail overlaps. 128_5 and 203_3 are single amplicons
    # where every read is contained in every other, so no overlap graph exists to
    # build at any parameter setting. Record it as an outcome so the summary can
    # tell it apart from hifiasm actually failing, and let the run exit 0.
    if [ ! -s "${gfa}" ]; then
        printf '0 contigs (empty overlap graph; amplicon input)\n' > "${outdir}/${srr}_metric.txt"
        echo "NOTE: hifiasm produced no contigs for ${srr} -- expected on single-amplicon input." >&2
        return 0
    fi
    local contigs="${outdir}/${srr}_contigs.fasta"   # FASTA the contigs get extracted into
    # convert GFA segment (S) lines into FASTA records
    awk '/^S/{print ">"$2"\n"$3}' "${gfa}" > "${contigs}"
    # a GFA with no S lines is the same outcome as no GFA at all, not a failure
    if [ ! -s "${contigs}" ]; then
        printf '0 contigs (GFA had no segments)\n' > "${outdir}/${srr}_metric.txt"
        echo "NOTE: no contigs parsed from ${gfa} for ${srr}." >&2
        return 0
    fi

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

    # Run each assembler (override the default set via ASSEMBLY_TOOLS). The two
    # minimap2 arms are the same pipeline over different references -- HXB2 as a
    # fixed subtype-B baseline, and the closest Group M panel entry -- so their
    # summary rows can be read against each other directly.
    for TOOL in ${ASSEMBLY_TOOLS:-minimap2_consensus minimap2_bestref hifiasm}; do
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
        # Time+run the matching run_<tool> function (name built from ${TOOL}).
        # -o pipefail matters: this child shell does NOT inherit the set at the
        # top of this script, so without it a failing `bcftools mpileup` piped
        # into `bcftools call` reports the exit status of `call` alone and the
        # `|| exit 1` never fires -- which is how an empty VCF once produced a
        # variant-free consensus that still looked like a clean run.
        measure_and_run "${TIMELOG}" -- \
            bash -o pipefail -c 'run_'"${TOOL}"' "$@"' _ "${SRR}" "${READS}" "${REF_FASTA}" "${OUTDIR}" > "${LOG}" 2>&1
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
            # Bases actually supported by reads. Total length is no longer a
            # useful validity test on its own: now that unsequenced positions are
            # N-masked rather than filled with reference, a consensus can be a
            # full 9719bp and still be 88% N, so the called count is what says
            # whether a genome was really recovered.
            CALLED=$(awk -v l="${LEN:-0}" -v n="${N_PCT:-100}" 'BEGIN{printf "%d", l*(100-n)/100}')
            # valid if enough of a full HIV genome (~9kb) was actually called
            if [ -n "${CALLED}" ] && [ "${CALLED}" -ge 8000 ]; then VALID=1; fi
            # which reference this arm mapped against, for side-by-side reading
            REFID=$(head -1 "${OUTDIR}/${SRR}_reference.txt" 2>/dev/null)
            METRIC="${LEN}bp, ${CALLED}bp called, ${N_PCT:-?}% N${REFID:+, ref=${REFID}}"
        elif [ -s "${OUTDIR}/${SRR}_metric.txt" ]; then
            # the tool ran to completion but produced no assembly, and said why
            METRIC=$(head -1 "${OUTDIR}/${SRR}_metric.txt")
        fi
        # write this assembler's row to summary.tsv
        append_summary_row "assembly_pacbio" "${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"
    done
done

# final confirmation pointing the user at the results table
echo "Done. See ${SUMMARY_TSV}"
