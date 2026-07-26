#!/bin/bash
# Compares jpHMM (per-sample recombination/subtype calls) vs IQ-TREE2
# (ML phylogeny to confirm subtype clustering) on the Illumina subset.
# COMET and REGA v3 are excluded (web-only, run those separately).
# Usage: ./compare_subtyping_oxfordnano.sh
set -uo pipefail                                     # -u errors on unset vars, pipefail catches mid-pipe failures; no -e so one tool failing doesn't abort the comparison

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"           # absolute path of this script's own dir, so paths resolve regardless of launch location
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"     # repo top-level (three dirs up), the base for every other path below
source "${REPO_ROOT}/scripts/common/lib_compare.sh"  # load shared helpers: measure_and_run, parse_time_metrics, append_summary_row, subset_accessions

ASSEMBLY_DIR="${REPO_ROOT}/results/assembly/oxfordnano/minimap2_out"  # input: per-sample consensus FASTAs from the Nanopore assembly stage
ALIGNMENT="${REPO_ROOT}/results/msa/oxfordnano/mafft_aligned.fasta"   # input: whole-subset MAFFT alignment from the msa stage (feeds IQ-TREE2)
RESULTS_DIR="${REPO_ROOT}/results/subtyping/oxfordnano"              # output: all subtyping results + summary go here
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"             # the single TSV every tool appends a timing/validity row to
mkdir -p "${RESULTS_DIR}"                            # create the results dir (and parents) if it doesn't exist
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"  # seed the manual notes file from the template only on first run

THREADS="${THREADS:-4}"                              # thread count; honour an externally-set THREADS, otherwise default to 4
export THREADS                                       # export so child scripts (run_jphmm.sh, run_iqtree.sh) inherit it

echo "=== jpHMM (per sample) ==="                    # progress marker for the per-sample subtyping pass
for SRR in $(subset_accessions nanopore "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do  # loop over just the Nanopore accessions chosen for this comparison
    QUERY="${ASSEMBLY_DIR}/${SRR}_consensus.fasta"   # the per-sample consensus jpHMM will subtype
    [ -s "${QUERY}" ] || { echo "NOTE: no consensus for ${SRR}, skipping (run scripts/assembly/oxfordnano/assembly_oxfordnano.sh first)." >&2; continue; }  # skip samples with no consensus yet

    OUTDIR="${RESULTS_DIR}/jphmm_out/${SRR}"         # per-sample jpHMM output dir
    TIMELOG="${RESULTS_DIR}/jphmm_${SRR}.time"       # file where measure_and_run records wallclock/RSS
    LOG="${RESULTS_DIR}/jphmm_${SRR}.log"            # captured stdout+stderr of the jpHMM run

    measure_and_run "${TIMELOG}" -- "${STAGE_DIR}/run_jphmm.sh" "${QUERY}" "${OUTDIR}" > "${LOG}" 2>&1  # time and run the jpHMM wrapper on this sample's consensus
    EXIT_CODE=$?                                     # capture the wrapper's exit status before $? is overwritten
    parse_time_metrics "${TIMELOG}"                  # set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file

    VALID=0                                          # assume invalid until proven otherwise
    METRIC="n/a"                                     # default human-readable metric
    RECOMB_FILE=$(compgen -G "${OUTDIR}/*recombination*" | head -1)  # find jpHMM's recombination output file (its key subtype result)
    if [ -n "${RECOMB_FILE}" ] && [ -s "${RECOMB_FILE}" ]; then  # valid only if that file exists and is non-empty
        VALID=1                                      # mark this run valid
        METRIC="see ${RECOMB_FILE}"                  # point the summary at the recombination result file
    fi
    append_summary_row "subtyping_oxfordnano" "jphmm" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"  # write this sample's jpHMM row to summary.tsv
done

echo "=== IQ-TREE2 (whole-subset ML tree, from the msa stage's MAFFT alignment) ==="  # progress marker for the phylogeny pass
if [ -s "${ALIGNMENT}" ]; then                       # only build a tree if the msa stage produced an alignment
    TIMELOG="${RESULTS_DIR}/iqtree.time"             # timing file for the IQ-TREE2 run
    LOG="${RESULTS_DIR}/iqtree.log"                  # log for the IQ-TREE2 run
    measure_and_run "${TIMELOG}" -- "${STAGE_DIR}/run_iqtree.sh" "${ALIGNMENT}" "${RESULTS_DIR}/iqtree_out" subset > "${LOG}" 2>&1  # time and run the IQ-TREE2 wrapper on the whole-subset alignment
    EXIT_CODE=$?                                     # capture IQ-TREE2's exit status
    parse_time_metrics "${TIMELOG}"                  # refresh WALLCLOCK_SEC / PEAK_RSS_MB from the timing file

    VALID=0                                          # assume invalid until proven otherwise
    METRIC="n/a"                                     # default human-readable metric
    TREEFILE="${RESULTS_DIR}/iqtree_out/subset.treefile"  # the ML tree IQ-TREE2 writes on success
    if [ -s "${TREEFILE}" ]; then                    # valid only if a non-empty treefile was produced
        VALID=1                                      # mark this run valid
        METRIC="tree written, $(grep -o "bp\." "${RESULTS_DIR}/iqtree_out/subset.iqtree" 2>/dev/null | wc -l) bootstrap refs"  # summarise: tree written plus a bootstrap-reference count from the report
    fi
    append_summary_row "subtyping_oxfordnano" "iqtree2" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"  # write the IQ-TREE2 row to summary.tsv
else
    echo "NOTE: no alignment at ${ALIGNMENT}, run scripts/msa/oxfordnano/compare_msa_oxfordnano.sh first." >&2  # tell the user to build the alignment first
fi

echo "Done. See ${SUMMARY_TSV}"                      # final confirmation pointing the user at the results table
