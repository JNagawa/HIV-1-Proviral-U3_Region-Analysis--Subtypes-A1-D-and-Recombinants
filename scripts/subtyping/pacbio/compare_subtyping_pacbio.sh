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
set -uo pipefail                                     # -u errors on unset vars, pipefail catches mid-pipe failures; no -e so one tool failing doesn't abort the comparison

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"           # absolute path of this script's own dir, so paths resolve regardless of launch location
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"     # repo top-level (three dirs up), the base for every other path below
source "${REPO_ROOT}/scripts/common/lib_compare.sh"  # load shared helpers: measure_and_run, parse_time_metrics, append_summary_row, subset_accessions
RUN_DIR="${REPO_ROOT}/scripts/subtyping/illumina"   # reuse shared run_<tool>.sh

ASSEMBLY_DIR="${REPO_ROOT}/results/assembly/pacbio/minimap2_consensus_out"  # input: per-sample proviral consensus FASTAs from the PacBio assembly stage
ALIGNMENT="${REPO_ROOT}/results/msa/pacbio/mafft_aligned.fasta"             # input: whole-subset MAFFT alignment from the PacBio msa stage (feeds IQ-TREE2)
RESULTS_DIR="${REPO_ROOT}/results/subtyping/pacbio"                         # output: all subtyping results + summary go here
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"             # the single TSV every tool appends a timing/validity row to
mkdir -p "${RESULTS_DIR}"                            # create the results dir (and parents) if it doesn't exist
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"  # seed the manual notes file from the template only on first run

THREADS="${THREADS:-4}"                              # thread count; honour an externally-set THREADS, otherwise default to 4
export THREADS                                       # export so child scripts (run_jphmm.sh, run_iqtree.sh) inherit it

echo "=== jpHMM (per sample) ==="                    # progress marker for the per-sample subtyping pass
for SRR in $(subset_accessions pacbio "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do  # loop over just the PacBio accessions chosen for this comparison
    QUERY="${ASSEMBLY_DIR}/${SRR}_consensus.fasta"   # the per-sample proviral consensus jpHMM will subtype
    [ -s "${QUERY}" ] || { echo "NOTE: no consensus for ${SRR}, skipping (run scripts/assembly/pacbio/compare_assembly_pacbio.sh first)." >&2; continue; }  # skip samples with no consensus yet

    OUTDIR="${RESULTS_DIR}/jphmm_out/${SRR}"         # per-sample jpHMM output dir
    TIMELOG="${RESULTS_DIR}/jphmm_${SRR}.time"       # file where measure_and_run records wallclock/RSS
    LOG="${RESULTS_DIR}/jphmm_${SRR}.log"            # captured stdout+stderr of the jpHMM run

    measure_and_run "${TIMELOG}" -- "${RUN_DIR}/run_jphmm.sh" "${QUERY}" "${OUTDIR}" > "${LOG}" 2>&1  # time and run the shared jpHMM wrapper on this sample's consensus
    EXIT_CODE=$?                                     # capture the wrapper's exit status before $? is overwritten
    parse_time_metrics "${TIMELOG}"                  # set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file

    VALID=0; METRIC="n/a"                            # assume invalid with a default metric until proven otherwise
    RECOMB_FILE=$(compgen -G "${OUTDIR}/*recombination*" | head -1)  # find jpHMM's recombination output file (its key subtype result)
    if [ -n "${RECOMB_FILE}" ] && [ -s "${RECOMB_FILE}" ]; then VALID=1; METRIC="see ${RECOMB_FILE}"; fi  # valid if that file exists and is non-empty; point summary at it
    append_summary_row "subtyping_pacbio" "jphmm" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"  # write this sample's jpHMM row to summary.tsv
done

echo "=== IQ-TREE2 (whole-subset ML tree, from the PacBio msa stage's MAFFT alignment) ==="  # progress marker for the phylogeny pass
if [ -s "${ALIGNMENT}" ]; then                       # only build a tree if the msa stage produced an alignment
    TIMELOG="${RESULTS_DIR}/iqtree.time"             # timing file for the IQ-TREE2 run
    LOG="${RESULTS_DIR}/iqtree.log"                  # log for the IQ-TREE2 run
    measure_and_run "${TIMELOG}" -- "${RUN_DIR}/run_iqtree.sh" "${ALIGNMENT}" "${RESULTS_DIR}/iqtree_out" subset > "${LOG}" 2>&1  # time and run the shared IQ-TREE2 wrapper on the whole-subset alignment
    EXIT_CODE=$?                                     # capture IQ-TREE2's exit status
    parse_time_metrics "${TIMELOG}"                  # refresh WALLCLOCK_SEC / PEAK_RSS_MB from the timing file

    VALID=0; METRIC="n/a"                            # assume invalid with a default metric until proven otherwise
    TREEFILE="${RESULTS_DIR}/iqtree_out/subset.treefile"  # the ML tree IQ-TREE2 writes on success
    if [ -s "${TREEFILE}" ]; then                    # valid only if a non-empty treefile was produced
        VALID=1                                      # mark this run valid
        METRIC="tree written, $(grep -o "bp\." "${RESULTS_DIR}/iqtree_out/subset.iqtree" 2>/dev/null | wc -l) bootstrap refs"  # summarise: tree written plus a bootstrap-reference count from the report
    fi
    append_summary_row "subtyping_pacbio" "iqtree2" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"  # write the IQ-TREE2 row to summary.tsv
else
    echo "NOTE: no alignment at ${ALIGNMENT}, run scripts/msa/pacbio/compare_msa_pacbio.sh first." >&2  # tell the user to build the alignment first
fi

echo "Done. See ${SUMMARY_TSV}"                      # final confirmation pointing the user at the results table
