#!/bin/bash
# Compares MAFFT vs MUSCLE vs Clustal Omega, aligning HXB2 + whatever
# per-sample consensus sequences are available from the assembly stage (preferring
# bwa_consensus_out since it doesn't depend on SPAdes/SHIVER succeeding).
# Usage: ./compare_msa_illumina.sh
set -uo pipefail                                     # -u errors on unset vars, pipefail fails a pipe if any stage fails; no -e so one aligner failing doesn't kill the comparison

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"           # absolute path of this script's own dir, so paths work regardless of launch location
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"     # repo top-level (three dirs up), the base for every other path below
source "${REPO_ROOT}/scripts/common/lib_compare.sh"  # load shared helpers: measure_and_run, parse_time_metrics, append_summary_row

REF_FASTA="${REPO_ROOT}/data/reference/K03455.1.fasta"                    # HXB2 (K03455.1) reference sequence, the alignment anchor
ASSEMBLY_DIR="${REPO_ROOT}/results/assembly/illumina/bwa_consensus_out"   # input: per-sample consensus FASTAs from the assembly stage
RESULTS_DIR="${REPO_ROOT}/results/msa/illumina"                           # output: alignments + summary go here
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"                                  # the single TSV every aligner appends a timing/validity row to
mkdir -p "${RESULTS_DIR}"                                                 # create the results dir (and parents) if it doesn't exist
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"  # seed the manual notes file from the template only on the first run

COMBINED="${RESULTS_DIR}/combined_input.fasta"       # single multi-FASTA (reference + all consensuses) fed to every aligner
cat "${REF_FASTA}" "${ASSEMBLY_DIR}"/*_consensus.fasta > "${COMBINED}" 2>/dev/null  # concatenate HXB2 + every sample consensus; silence errors if some are missing
N_SEQS=$(grep -c "^>" "${COMBINED}" 2>/dev/null || echo 0)  # count sequences (FASTA headers) in the combined input
if [ "${N_SEQS}" -lt 2 ]; then                       # an MSA needs at least two sequences to be meaningful
    echo "ERROR: fewer than 2 sequences available (run scripts/assembly/illumina/compare_assembly_illumina.sh (bwa_consensus) first). Found ${N_SEQS}." >&2  # tell the user what to run first
    exit 1                                           # bail out; nothing to align
fi
echo "Aligning ${N_SEQS} sequences (HXB2 + subset consensus sequences)."  # progress marker in the log

THREADS="${THREADS:-4}"                              # thread count; honour an externally-set THREADS, otherwise default to 4
export THREADS                                       # export so the run_<tool>.sh child scripts inherit it

for TOOL in mafft muscle clustalo; do                # run each aligner on the same input for a head-to-head comparison
    OUT="${RESULTS_DIR}/${TOOL}_aligned.fasta"       # this aligner's output alignment
    TIMELOG="${RESULTS_DIR}/${TOOL}.time"            # file where measure_and_run records wallclock/RSS
    LOG="${RESULTS_DIR}/${TOOL}.log"                 # captured stdout+stderr of the aligner

    echo "=== ${TOOL} ==="                           # progress marker in the log
    measure_and_run "${TIMELOG}" -- "${STAGE_DIR}/run_${TOOL}.sh" "${COMBINED}" "${OUT}" "${LOG}"  # time+run the per-tool wrapper on the combined FASTA
    EXIT_CODE=$?                                      # capture the aligner's exit status before $? is overwritten
    parse_time_metrics "${TIMELOG}"                  # set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file

    VALID=0                                          # assume invalid until proven otherwise
    METRIC="n/a"                                     # default human-readable metric
    if [ -s "${OUT}" ]; then                         # only inspect the alignment if it was actually produced
        OUT_SEQS=$(grep -c "^>" "${OUT}")            # number of sequences present in the output alignment
        COLS=$(seqkit stats -T "${OUT}" 2>/dev/null | tail -1 | cut -f7 | tr -d ',')  # alignment width (column count) from seqkit's max-length field
        GAP_CHARS=$(grep -v "^>" "${OUT}" | tr -cd '-' | wc -c)      # total gap characters across all aligned sequences
        TOTAL_CHARS=$(grep -v "^>" "${OUT}" | tr -d '\n' | wc -c)    # total residue+gap characters, the denominator for gap %
        GAP_PCT=$(awk -v g="${GAP_CHARS}" -v t="${TOTAL_CHARS}" 'BEGIN{if(t>0) printf "%.1f", 100*g/t; else print "0"}')  # gap percentage, a rough alignment-quality proxy
        if [ "${OUT_SEQS}" = "${N_SEQS}" ]; then     # valid only if the aligner kept every input sequence
            VALID=1                                  # mark this run as valid
        fi
        METRIC="${OUT_SEQS}/${N_SEQS} seqs retained, ${COLS:-?} columns, ${GAP_PCT}% gap"  # record the key metrics for the summary
    fi

    append_summary_row "msa_illumina" "${TOOL}" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"  # write this aligner's row to summary.tsv
done

echo "Done. See ${SUMMARY_TSV}"                      # final confirmation pointing the user at the results table
