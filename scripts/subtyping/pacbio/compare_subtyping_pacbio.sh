#!/bin/bash
# PacBio arm: compares jpHMM (per-sample recombination/subtype calls) vs
# IQ-TREE2 (ML phylogeny confirming subtype clustering) on the PacBio subset.
# Reuses the Illumina arm's run_jphmm.sh / run_iqtree.sh unchanged (subtyping
# is platform-agnostic) -- jpHMM scans each per-sample proviral consensus from
# minimap2_consensus_out, IQ-TREE2 builds one ML tree from the PacBio msa
# stage's MAFFT alignment. COMET and REGA v3 are web-only (run separately and
# fold into the notes by hand). This is where the A1/D/recombinant subtype of
# each Rakai cohort sample is actually determined (it is not in the SRA
# metadata); the 92UG/93UG anchors give known-subtype reference points.
# Usage: ./compare_subtyping_pacbio.sh   (via sbatch scripts/utils/run_comparison_step.slurm.sh)
# -u errors on unset vars, pipefail catches mid-pipe failures; no -e so one tool failing doesn't
# abort the comparison
set -uo pipefail

# absolute path of this script's own dir, so paths resolve regardless of launch location
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
# repo top-level (three dirs up), the base for every other path below
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
# load shared helpers: measure_and_run, parse_time_metrics, append_summary_row, subset_accessions
source "${REPO_ROOT}/scripts/common/lib_compare.sh"
RUN_DIR="${REPO_ROOT}/scripts/subtyping/illumina"   # reuse shared run_<tool>.sh

# input: per-sample proviral consensus FASTAs from the PacBio assembly stage
# which assembly arm feeds this stage (see compare_msa_pacbio.sh for the rationale)
ASSEMBLY_ARM="${ASSEMBLY_ARM:-minimap2_consensus}"
ASSEMBLY_DIR="${REPO_ROOT}/results/assembly/pacbio/${ASSEMBLY_ARM}_out"
# input: whole-subset MAFFT alignment from the PacBio msa stage (feeds IQ-TREE2)
ALIGNMENT="${REPO_ROOT}/results/msa/pacbio/${ASSEMBLY_ARM}/mafft_aligned.fasta"
# output: all subtyping results + summary go here
RESULTS_DIR="${REPO_ROOT}/results/subtyping/pacbio/${ASSEMBLY_ARM}"
# the single TSV every tool appends a timing/validity row to
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
# create the results dir (and parents) if it doesn't exist
mkdir -p "${RESULTS_DIR}"
# seed the manual notes file from the template only on first run
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

# thread count; honour an externally-set THREADS, otherwise default to 4
THREADS="${THREADS:-4}"
# export so child scripts (run_jphmm.sh, run_iqtree.sh) inherit it
export THREADS

# progress marker for the per-sample subtyping pass
echo "=== jpHMM (per sample) ==="
# loop over just the PacBio accessions chosen for this comparison
for SRR in $(subset_accessions pacbio "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do
    # the per-sample proviral consensus jpHMM will subtype
    QUERY="${ASSEMBLY_DIR}/${SRR}_consensus.fasta"
    # skip samples with no consensus yet
    [ -s "${QUERY}" ] || { echo "NOTE: no consensus for ${SRR}, skipping (run scripts/assembly/pacbio/compare_assembly_pacbio.sh first)." >&2; continue; }

    OUTDIR="${RESULTS_DIR}/jphmm_out/${SRR}"         # per-sample jpHMM output dir
    # file where measure_and_run records wallclock/RSS
    TIMELOG="${RESULTS_DIR}/jphmm_${SRR}.time"
    LOG="${RESULTS_DIR}/jphmm_${SRR}.log"            # captured stdout+stderr of the jpHMM run

    # time and run the shared jpHMM wrapper on this sample's consensus
    measure_and_run "${TIMELOG}" -- "${RUN_DIR}/run_jphmm.sh" "${QUERY}" "${OUTDIR}" > "${LOG}" 2>&1
    # capture the wrapper's exit status before $? is overwritten
    EXIT_CODE=$?
    # set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file
    parse_time_metrics "${TIMELOG}"

    # assume invalid with a default metric until proven otherwise
    VALID=0; METRIC="n/a"
    # find jpHMM's recombination output file (its key subtype result)
    RECOMB_FILE=$(compgen -G "${OUTDIR}/*recombination*" | head -1)
    # valid if that file exists and is non-empty; point summary at it
    if [ -n "${RECOMB_FILE}" ] && [ -s "${RECOMB_FILE}" ]; then VALID=1; METRIC="see ${RECOMB_FILE}"; fi
    # write this sample's jpHMM row to summary.tsv
    append_summary_row "subtyping_pacbio" "jphmm" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"
done

# progress marker for the phylogeny pass
echo "=== IQ-TREE2 (whole-subset ML tree, from the PacBio msa stage's MAFFT alignment) ==="
# only build a tree if the msa stage produced an alignment
if [ -s "${ALIGNMENT}" ]; then
    TIMELOG="${RESULTS_DIR}/iqtree.time"             # timing file for the IQ-TREE2 run
    LOG="${RESULTS_DIR}/iqtree.log"                  # log for the IQ-TREE2 run
    # time and run the shared IQ-TREE2 wrapper on the whole-subset alignment
    measure_and_run "${TIMELOG}" -- "${RUN_DIR}/run_iqtree.sh" "${ALIGNMENT}" "${RESULTS_DIR}/iqtree_out" subset > "${LOG}" 2>&1
    EXIT_CODE=$?                                     # capture IQ-TREE2's exit status
    # refresh WALLCLOCK_SEC / PEAK_RSS_MB from the timing file
    parse_time_metrics "${TIMELOG}"

    # assume invalid with a default metric until proven otherwise
    VALID=0; METRIC="n/a"
    TREEFILE="${RESULTS_DIR}/iqtree_out/subset.treefile"  # the ML tree IQ-TREE2 writes on success
    # valid only if a non-empty treefile was produced
    if [ -s "${TREEFILE}" ]; then
        VALID=1                                      # mark this run valid
        # summarise: tree written plus a bootstrap-reference count from the report
        METRIC="tree written, $(grep -o "bp\." "${RESULTS_DIR}/iqtree_out/subset.iqtree" 2>/dev/null | wc -l) bootstrap refs"
    fi
    # write the IQ-TREE2 row to summary.tsv
    append_summary_row "subtyping_pacbio" "iqtree2" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"
else
    # tell the user to build the alignment first
    echo "NOTE: no alignment at ${ALIGNMENT}, run scripts/msa/pacbio/compare_msa_pacbio.sh first." >&2
fi

# final confirmation pointing the user at the results table
echo "Done. See ${SUMMARY_TSV}"
