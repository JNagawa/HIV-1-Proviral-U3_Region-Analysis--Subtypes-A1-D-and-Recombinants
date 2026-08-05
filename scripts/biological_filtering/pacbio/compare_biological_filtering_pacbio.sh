#!/bin/bash
# PacBio arm: compares Poplars (Hypermut 3) vs HIVSeqinR vs HIVIntact on the
# PacBio assembled subset (HXB2 + per-sample proviral consensus from
# minimap2_consensus_out). Reuses the Illumina arm's run_<tool>.sh scripts and
# the shared scripts/tools/{Poplars,HIVSeqinR,HIVIntact} clones -- intactness/
# hypermutation classification is platform-agnostic (operates on FASTA
# assemblies), only the input consensus sequences differ. As on the Illumina
# arm, no gold standard exists for this cohort, so "validity" is just "ran and
# produced a classification"; the review's recommended order is Poplars then
# HIVSeqinR, with HIVIntact as an optional cross-check (subtype-B-biased, so
# read with caution on A1/D).
# Usage: ./compare_biological_filtering_pacbio.sh   (via sbatch scripts/utils/run_comparison_step.slurm.sh)
# -u errors on unset vars, pipefail fails a pipe if any stage fails; no -e so one tool failing
# doesn't kill the comparison
set -uo pipefail

# absolute path of this script's own dir, so paths work regardless of launch location
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
# repo top-level (three dirs up), the base for every other path below
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
# load shared helpers: measure_and_run, parse_time_metrics, append_summary_row
source "${REPO_ROOT}/scripts/common/lib_compare.sh"
RUN_DIR="${REPO_ROOT}/scripts/biological_filtering/illumina"   # reuse shared run_<tool>.sh

# HXB2 reference genome, used as the first record (reference) for the tools
REF_FASTA="${REPO_ROOT}/data/reference/K03455.1.fasta"
# input: per-sample proviral consensus FASTAs from the PacBio assembly stage
# which assembly arm feeds this stage (see compare_msa_pacbio.sh for the rationale)
ASSEMBLY_ARM="${ASSEMBLY_ARM:-minimap2_consensus}"
ASSEMBLY_DIR="${REPO_ROOT}/results/assembly/pacbio/${ASSEMBLY_ARM}_out"
# output: all classification results + summary go here
RESULTS_DIR="${REPO_ROOT}/results/biological_filtering/pacbio/${ASSEMBLY_ARM}"
# the single TSV every tool appends a timing/validity row to
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
# create the results dir (and parents) if it doesn't exist
mkdir -p "${RESULTS_DIR}"
# seed the manual notes file from the template only on first run
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

# single FASTA fed to all three tools (HXB2 first, then consensuses)
INPUT_FASTA="${RESULTS_DIR}/combined_input.fasta"
# concatenate reference + all consensus seqs; HXB2 first so Poplars uses it as reference
cat "${REF_FASTA}" "${ASSEMBLY_DIR}"/*_consensus.fasta > "${INPUT_FASTA}" 2>/dev/null
# count FASTA records to check we have enough input
N_SEQS=$(grep -c "^>" "${INPUT_FASTA}" 2>/dev/null || echo 0)
# need at least the reference plus one consensus to compare
if [ "${N_SEQS}" -lt 2 ]; then
    # tell the user what to run first
    echo "ERROR: fewer than 2 sequences available (run scripts/assembly/pacbio/compare_assembly_pacbio.sh first). Found ${N_SEQS}." >&2
    exit 1                                           # bail; nothing to classify
fi
if [ ! -d "${REPO_ROOT}/scripts/tools/Poplars" ]; then  # the tool clones are a prerequisite
    # point the user at the shared setup script
    echo "NOTE: tools not set up yet -- run scripts/biological_filtering/illumina/setup_tools.sh first." >&2
    exit 1                                           # bail; tools not installed
fi

echo "=== Poplars (Hypermut 3) ==="                  # progress marker in the log
OUT="${RESULTS_DIR}/poplars_out.tsv"                 # Poplars hypermutation results table
# file where measure_and_run records wallclock/RSS
TIMELOG="${RESULTS_DIR}/poplars.time"
# Poplars' hypermut compares each sequence to a consensus column by column and
# asserts every record is the same length ("Sequences are not aligned."). The raw
# consensuses satisfied that only by accident, back when zero indels were being
# called and every one came out at exactly the reference length; now that indels
# are applied they differ (9719-9751bp) and the assertion fires. An ALIGNED FASTA
# is the correct input for a positional comparison anyway, so this stage now
# prefers the MSA stage's MAFFT alignment for this arm and only falls back to the
# unaligned consensuses if it is missing.
POPLARS_IN="${REPO_ROOT}/results/msa/pacbio/${ASSEMBLY_ARM}/mafft_aligned.fasta"
if [ ! -s "${POPLARS_IN}" ]; then
    echo "NOTE: ${POPLARS_IN} missing; falling back to unaligned consensuses, which Poplars may reject." >&2
    POPLARS_IN="${INPUT_FASTA}"
fi
# run the shared Poplars wrapper under timing, capturing output
measure_and_run "${TIMELOG}" -- "${RUN_DIR}/run_poplars.sh" "${POPLARS_IN}" "${OUT}" > "${RESULTS_DIR}/poplars.log" 2>&1
# capture the wrapper's exit status before $? is overwritten
EXIT_CODE=$?
# set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file
parse_time_metrics "${TIMELOG}"
VALID=0; METRIC="n/a"                                 # assume invalid until proven otherwise
# valid if output non-empty; record row count as the key metric
if [ -s "${OUT}" ]; then VALID=1; METRIC="$(wc -l < "${OUT}") result rows"; fi
# write Poplars' row to summary.tsv
append_summary_row "biological_filtering_pacbio" "poplars" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"

echo "=== HIVSeqinR ==="                              # progress marker in the log
OUTDIR="${RESULTS_DIR}/hivseqinr_out"                 # per-tool output dir for HIVSeqinR
TIMELOG="${RESULTS_DIR}/hivseqinr.time"               # timing file for this run
# run the shared HIVSeqinR wrapper under timing, capturing output
measure_and_run "${TIMELOG}" -- "${RUN_DIR}/run_hivseqinr.sh" "${INPUT_FASTA}" "${OUTDIR}" > "${RESULTS_DIR}/hivseqinr_wrapper.log" 2>&1
EXIT_CODE=$?                                          # capture the wrapper's exit status
# refresh WALLCLOCK_SEC / PEAK_RSS_MB from the timing file
parse_time_metrics "${TIMELOG}"
VALID=0; METRIC="n/a"                                 # default to invalid until checked
CSV="${OUTDIR}/Output_MyBigSummary_DF_FINAL.csv"      # HIVSeqinR's final summary table
# valid if CSV exists; subtract 1 for the header row to get sequences classified
if [ -s "${CSV}" ]; then
    VALID=1; METRIC="$(($(wc -l < "${CSV}") - 1)) classified"
else
    # Record WHY it failed instead of a bare "n/a". HIVSeqinR is an RStudio-era
    # script being driven headlessly on input it was not designed for (amplicon
    # segments, no PCR primer flanks), so a reason in the summary is the
    # difference between a diagnosable result and an unexplained blank.
    REASON=$(grep -m1 -A2 "^Error" "${OUTDIR}/hivseqinr.log" 2>/dev/null \
             | tr '\n' ' ' | tr -s ' ' | cut -c1-160)
    METRIC="${REASON:-no summary CSV produced; see hivseqinr.log}"
fi
# write HIVSeqinR's row to summary.tsv
append_summary_row "biological_filtering_pacbio" "hivseqinr" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"

echo "=== HIVIntact ==="                              # progress marker in the log
OUTDIR="${RESULTS_DIR}/hivintact_out"                 # per-tool output dir for HIVIntact
TIMELOG="${RESULTS_DIR}/hivintact.time"               # timing file for this run
# HIVIntact scores a sequence against a subtype-specific reference, so the
# --subtype flag has to match the sequence being scored. Running the whole batch
# as A1 (the previous behaviour) scored HXB2 itself against an A1 reference and
# duly reported it as defective -- FrameshiftInOrf, MajorSpliceDonorSiteMutated
# and two MisplacedORFs on the canonical reference genome, which is the signature
# of a misconfigured run rather than a finding. Measured 2026-08-05.
#
# Each sequence is therefore scored on its own, with the subtype taken from the
# reference selected for it (<sample>_reference.txt, e.g. "D:KU168271" -> D).
#
# The subtype is read from the minimap2_bestref arm regardless of which arm is
# being processed, because that arm's reference selection IS our subtype estimate.
# Reading it from the current arm would make every sample "HXB2" in the HXB2
# baseline arm, which is true of the coordinate frame but not of the biology --
# these are subtype A1/D sequences whichever frame they are expressed in.
# HXB2 itself is scored as HXB2 and acts as the positive control: if it does not
# come back intact, this stage is not to be trusted.
SUBTYPE_SRC_DIR="${REPO_ROOT}/results/assembly/pacbio/minimap2_bestref_out"
mkdir -p "${OUTDIR}"
# HIVIntact only accepts subtypes it ships a reference for
SUPPORTED_SUBTYPES="A1 A2 B C D F1 F2 G H HXB2"
INTACT_MANIFEST="${OUTDIR}/per_sample_subtype.tsv"
printf 'sequence\tsubtype_used\tsource\n' > "${INTACT_MANIFEST}"

: > "${RESULTS_DIR}/hivintact_wrapper.log"            # start a fresh wrapper log
SPLIT_DIR="${OUTDIR}/per_sample_input"
rm -rf "${SPLIT_DIR}"; mkdir -p "${SPLIT_DIR}"
# split the combined FASTA into one file per record so each can carry its own subtype
awk -v d="${SPLIT_DIR}" '
    /^>/ { id = substr($0,2); split(id, t, /[ \t]/); id = t[1]
           gsub(/[^A-Za-z0-9_.]/, "_", id)
           f = d "/" id ".fasta"; print > f; next }
    { print > f }
' "${INPUT_FASTA}"

measure_and_run "${TIMELOG}" -- bash -o pipefail -c '
for FA in "$1"/*.fasta; do
    ID=$(basename "${FA}" .fasta)
    case "${ID}" in
        K03455*) SUB=HXB2; SRC="reference control" ;;
        *)
            REFTXT="$2/${ID}_reference.txt"
            RAW=$(head -1 "${REFTXT}" 2>/dev/null)
            # "D:KU168271" -> D ; "HXB2" -> HXB2 ; anything else falls back
            CAND="${RAW%%:*}"
            case " $3 " in
                *" ${CAND} "*) SUB="${CAND}"; SRC="assembly reference ${RAW}" ;;
                *) SUB=HXB2; SRC="unsupported/absent (${RAW:-none}) -- fell back to HXB2" ;;
            esac
            ;;
    esac
    printf "%s\t%s\t%s\n" "${ID}" "${SUB}" "${SRC}" >> "$4"
    echo "--- HIVIntact ${ID} as subtype ${SUB} (${SRC})"
    "$5/run_hivintact.sh" "${FA}" "$6/${ID}" "${SUB}" >> "$7" 2>&1 || echo "WARNING: HIVIntact failed for ${ID}"
done
' _ "${SPLIT_DIR}" "${SUBTYPE_SRC_DIR}" "${SUPPORTED_SUBTYPES}" "${INTACT_MANIFEST}" "${RUN_DIR}" "${OUTDIR}" "${RESULTS_DIR}/hivintact_wrapper.log"
EXIT_CODE=$?                                          # capture the wrapper's exit status

# Merge the per-sample FASTA outputs back into the batch-level filenames the
# validity check below and the downstream report already expect. Concatenating
# FASTA is safe; the per-sample errors.json files are deliberately NOT merged --
# stitching JSON objects together with text tools is fragile, and the report
# reads each sample's own errors.json instead.
cat "${OUTDIR}"/*/intact.fasta    > "${OUTDIR}/intact.fasta"    2>/dev/null
cat "${OUTDIR}"/*/nonintact.fasta > "${OUTDIR}/nonintact.fasta" 2>/dev/null
# refresh WALLCLOCK_SEC / PEAK_RSS_MB from the timing file
parse_time_metrics "${TIMELOG}"
VALID=0; METRIC="n/a"                                 # default to invalid until checked
# valid if either the intact or non-intact output was produced
if [ -s "${OUTDIR}/intact.fasta" ] || [ -s "${OUTDIR}/nonintact.fasta" ]; then
    VALID=1                                           # mark this run as valid
    # count intact proviruses (0 if file missing)
    N_INTACT=$(grep -c "^>" "${OUTDIR}/intact.fasta" 2>/dev/null); N_INTACT="${N_INTACT:-0}"
    # count non-intact proviruses (0 if file missing)
    N_NONINTACT=$(grep -c "^>" "${OUTDIR}/nonintact.fasta" 2>/dev/null); N_NONINTACT="${N_NONINTACT:-0}"
    # record the intact/non-intact split as the key metric
    METRIC="${N_INTACT} intact, ${N_NONINTACT} non-intact"
fi
# write HIVIntact's row to summary.tsv
append_summary_row "biological_filtering_pacbio" "hivintact" "all_subset" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"

# Generate the basis-for-verdict report. summary.tsv records the intact/non-intact
# split; the report pairs each verdict with per-gene read coverage so a "defective"
# call cannot be confused with a gene that was simply never sequenced. Non-fatal:
# a failed report must not discard classifications that succeeded.
if ! bash "${REPO_ROOT}/scripts/utils/report_pacbio_biological_filtering.sh" "${ASSEMBLY_ARM}"; then
    echo "WARNING: intactness basis report failed; summary.tsv is still valid." >&2
fi

# final confirmation pointing the user at the results table
echo "Done. See ${SUMMARY_TSV}"
