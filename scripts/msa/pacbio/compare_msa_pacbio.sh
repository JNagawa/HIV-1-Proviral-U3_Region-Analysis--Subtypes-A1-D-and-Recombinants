#!/bin/bash
# PacBio arm: compares MAFFT (L-INS-i) vs MUSCLE vs Clustal Omega, aligning
# HXB2 + the per-sample proviral consensus sequences from the PacBio assembly
# stage (minimap2_consensus_out, the reference-guided output that doesn't
# depend on hifiasm succeeding). Reuses the Illumina arm's run_<tool>.sh
# scripts unchanged -- alignment is platform-agnostic, only the input consensus
# sequences differ. The review prefers MAFFT L-INS-i for the final coordinate
# alignment; MUSCLE/Clustal Omega are the speed/accuracy comparators.
# Usage: ./compare_msa_pacbio.sh   (via sbatch scripts/utils/run_comparison_step.slurm.sh)
# -u errors on unset vars, pipefail fails a pipe if any stage fails; no -e so one aligner failing
# doesn't kill the comparison
set -uo pipefail

# absolute path of this script's own dir, so paths work regardless of launch location
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
# repo top-level (three dirs up), the base for every other path below
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
# load shared helpers: measure_and_run, parse_time_metrics, append_summary_row
source "${REPO_ROOT}/scripts/common/lib_compare.sh"
RUN_DIR="${REPO_ROOT}/scripts/msa/illumina"   # reuse the shared run_<tool>.sh scripts

# HXB2 (K03455.1) reference sequence, the alignment anchor
REF_FASTA="${REPO_ROOT}/data/reference/K03455.1.fasta"
# Which assembly arm feeds this stage. The PacBio assembly step produces two
# reference-guided consensus sets -- minimap2_consensus (HXB2 baseline) and
# minimap2_bestref (closest LTR-complete Group M reference) -- and both are
# carried through downstream so the effect of the reference choice stays
# measurable all the way to the motif hits. Each arm writes to its own
# subdirectory rather than sharing one summary.tsv with an extra column,
# because append_summary_row's 9-column schema is shared by every stage on all
# three platforms and widening it would desync the existing summaries.
ASSEMBLY_ARM="${ASSEMBLY_ARM:-minimap2_consensus}"
# input: per-sample proviral consensus FASTAs from the PacBio assembly stage
ASSEMBLY_DIR="${REPO_ROOT}/results/assembly/pacbio/${ASSEMBLY_ARM}_out"
# output: alignments + summary go here, scoped to this arm
RESULTS_DIR="${REPO_ROOT}/results/msa/pacbio/${ASSEMBLY_ARM}"
# the single TSV every aligner appends a timing/validity row to
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
# create the results dir (and parents) if it doesn't exist
mkdir -p "${RESULTS_DIR}"
# seed the manual notes file from the template only on the first run
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

# single multi-FASTA (reference + all consensuses) fed to every aligner
COMBINED="${RESULTS_DIR}/combined_input.fasta"
# concatenate HXB2 + every sample consensus; silence errors if some are missing
cat "${REF_FASTA}" "${ASSEMBLY_DIR}"/*_consensus.fasta > "${COMBINED}" 2>/dev/null
# count sequences (FASTA headers) in the combined input
N_SEQS=$(grep -c "^>" "${COMBINED}" 2>/dev/null || echo 0)
# an MSA needs at least two sequences to be meaningful
if [ "${N_SEQS}" -lt 2 ]; then
    # tell the user what to run first
    echo "ERROR: fewer than 2 sequences available (run scripts/assembly/pacbio/compare_assembly_pacbio.sh first). Found ${N_SEQS}." >&2
    exit 1                                           # bail out; nothing to align
fi
# progress marker in the log
echo "Aligning ${N_SEQS} sequences (HXB2 + PacBio subset consensus sequences)."

# thread count; honour an externally-set THREADS, otherwise default to 4
THREADS="${THREADS:-4}"
# export so the run_<tool>.sh child scripts inherit it
export THREADS

# run each aligner on the same input for a head-to-head comparison
for TOOL in mafft muscle clustalo; do
    OUT="${RESULTS_DIR}/${TOOL}_aligned.fasta"       # this aligner's output alignment
    # file where measure_and_run records wallclock/RSS
    TIMELOG="${RESULTS_DIR}/${TOOL}.time"
    LOG="${RESULTS_DIR}/${TOOL}.log"                 # captured stdout+stderr of the aligner

    echo "=== ${TOOL} ==="                           # progress marker in the log
    # time+run the reused Illumina per-tool wrapper on the combined FASTA
    measure_and_run "${TIMELOG}" -- "${RUN_DIR}/run_${TOOL}.sh" "${COMBINED}" "${OUT}" "${LOG}"
    # capture the aligner's exit status before $? is overwritten
    EXIT_CODE=$?
    # set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file
    parse_time_metrics "${TIMELOG}"

    # assume invalid with a default metric until proven otherwise
    VALID=0; METRIC="n/a"
    # only inspect the alignment if it was actually produced
    if [ -s "${OUT}" ]; then
        # number of sequences present in the output alignment
        OUT_SEQS=$(grep -c "^>" "${OUT}")
        # alignment width (column count) from seqkit's max-length field
        COLS=$(seqkit stats -T "${OUT}" 2>/dev/null | tail -1 | cut -f7 | tr -d ',')
        # total gap characters across all aligned sequences
        GAP_CHARS=$(grep -v "^>" "${OUT}" | tr -cd '-' | wc -c)
        # total residue+gap characters, the denominator for gap %
        TOTAL_CHARS=$(grep -v "^>" "${OUT}" | tr -d '\n' | wc -c)
        # gap percentage, a rough alignment-quality proxy
        GAP_PCT=$(awk -v g="${GAP_CHARS}" -v t="${TOTAL_CHARS}" 'BEGIN{if(t>0) printf "%.1f", 100*g/t; else print "0"}')
        # valid only if the aligner kept every input sequence
        [ "${OUT_SEQS}" = "${N_SEQS}" ] && VALID=1
        # record the key metrics for the summary
        METRIC="${OUT_SEQS}/${N_SEQS} seqs retained, ${COLS:-?} columns, ${GAP_PCT}% gap"
    fi
    # write this aligner's row to summary.tsv
    append_summary_row "msa_pacbio" "${TOOL}" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"
done

# final confirmation pointing the user at the results table
echo "Done. See ${SUMMARY_TSV}"
