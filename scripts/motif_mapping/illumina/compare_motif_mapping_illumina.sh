#!/bin/bash
# Compares FIMO vs MOODS vs TFBSTools -- PWM motif scanning for the 6 core
# JASPAR TFs (see setup_jaspar.sh) -- and separately gquad vs pqsfinder --
# G-quadruplex prediction, which doesn't use JASPAR PWMs at all, just the
# raw sequence -- against U3 sequences extracted from the msa stage's
# alignment via scripts/utils/extract_u3_by_hxb2_anchor.sh.
#
# Per the README, each tool runs a positive control first: scanning HXB2's
# own U3 alone. For the TFBS tools this is a real positive control (JASPAR
# motifs are known/expected to be present); for gquad/pqsfinder it's a
# weaker sanity check (no literature-confirmed G4 site in HXB2's U3 is
# cited anywhere in this repo) -- both G4 tools finding zero hits there
# isn't itself damning the way it would be for the TFBS tools, but a
# non-zero, non-implausible result is still worth confirming before
# trusting the subset numbers.
# Usage: ./compare_motif_mapping_illumina.sh
# -u errors on unset vars, pipefail catches pipeline failures; no -e so one tool failing doesn't
# abort the comparison
set -uo pipefail

# absolute path of this script's dir, so it runs from anywhere
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
# repo top-level (three dirs up), base for every path below
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
# load shared helpers: measure_and_run, parse_time_metrics, append_summary_row
source "${REPO_ROOT}/scripts/common/lib_compare.sh"

# input: MAFFT alignment from the msa stage
ALIGNMENT="${REPO_ROOT}/results/msa/illumina/mafft_aligned.fasta"
# output: all motif-mapping results + summary go here
RESULTS_DIR="${REPO_ROOT}/results/motif_mapping/illumina"
# the single TSV every tool appends a timing/validity row to
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
# create the results dir (and parents) if needed
mkdir -p "${RESULTS_DIR}"
# seed the manual notes file from the template on first run only
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

# can't map motifs without the upstream alignment...
if [ ! -s "${ALIGNMENT}" ]; then
    # ...tell the user which stage to run first...
    echo "ERROR: no alignment at ${ALIGNMENT} (run scripts/msa/illumina/compare_msa_illumina.sh first)." >&2
    exit 1                                           # ...and abort
fi

# progress marker in the log
echo "=== Extracting U3 (anchored to HXB2's own annotated LTR/R-region) ==="
# ungapped U3 sequences, the actual scanning input for every tool
U3_FASTA="${RESULTS_DIR}/U3_extracted.fasta"
# slice U3 out of each aligned seq using HXB2's annotated LTR/R boundary as the anchor
bash "${REPO_ROOT}/scripts/utils/extract_u3_by_hxb2_anchor.sh" \
    --alignment "${ALIGNMENT}" \
    --hxb2-id K03455.1 \
    --gb-cache "${REPO_ROOT}/data/reference/K03455.1.gb" \
    --out-gapped "${RESULTS_DIR}/U3_aligned.fasta" \
    --out "${U3_FASTA}" \
    --warnings-log "${RESULTS_DIR}/u3_extraction_warnings.log"
# if extraction yielded nothing there's nothing to scan...
if [ ! -s "${U3_FASTA}" ]; then
    # ...report it...
    echo "ERROR: U3 extraction produced no sequences, see extraction output above." >&2
    exit 1                                           # ...and abort
fi

# HXB2's own U3 alone, used as each tool's positive control
HXB2_U3_ONLY="${RESULTS_DIR}/hxb2_u3_only.fasta"
# pull just the HXB2 record (ID starting K03455) out of the U3 set
seqkit grep -n -r -p "^K03455" "${U3_FASTA}" > "${HXB2_U3_ONLY}"
# if HXB2's U3 couldn't be isolated, the positive control is impossible...
if [ ! -s "${HXB2_U3_ONLY}" ]; then
    # ...report it...
    echo "ERROR: could not isolate HXB2's own U3 from ${U3_FASTA} for the positive control." >&2
    exit 1                                           # ...and abort
fi

# run one tool on one input, time it, count hits, and log a summary row
run_tool() {
    # args: tool name, sample label, input FASTA, output path
    local tool="$1" sample="$2" input="$3" out="$4"
    local timelog="${RESULTS_DIR}/${tool}_${sample}.time"  # per-run timing file for measure_and_run
    local log="${RESULTS_DIR}/${tool}_${sample}.log"       # per-run captured stdout+stderr
    local exit_code n_hits valid metric              # locals set per tool below

    # dispatch to the right wrapper and hit-counting rule per tool
    case "${tool}" in
        fimo)
            # run FIMO under timing, capturing output
            measure_and_run "${timelog}" -- "${STAGE_DIR}/run_fimo.sh" "${input}" "${out}" > "${log}" 2>&1
            # capture its exit status before $? is overwritten
            exit_code=$?
            n_hits=0                                  # default to zero hits
            # count data rows in fimo.tsv, excluding comments and the header
            [ -s "${out}/fimo.tsv" ] && n_hits=$(grep -vc "^#\|^motif_id" "${out}/fimo.tsv" 2>/dev/null || echo 0)
            ;;
        moods)
            # run MOODS under timing, capturing output
            measure_and_run "${timelog}" -- "${STAGE_DIR}/run_moods.sh" "${input}" "${out}" > "${log}" 2>&1
            exit_code=$?                              # capture its exit status
            n_hits=0                                  # default to zero hits
            # MOODS output has one hit per line, so line count = hit count
            [ -s "${out}" ] && n_hits=$(wc -l < "${out}")
            ;;
        tfbstools)
            # run TFBSTools under timing, capturing output
            measure_and_run "${timelog}" -- "${STAGE_DIR}/run_tfbstools.sh" "${input}" "${out}" > "${log}" 2>&1
            exit_code=$?                              # capture its exit status
            n_hits=0                                  # default to zero hits
            # count non-comment GFF3 lines as hits
            [ -s "${out}" ] && n_hits=$(grep -vc "^#" "${out}" 2>/dev/null || echo 0)
            ;;
        gquad)
            # run gquad G4 prediction under timing
            measure_and_run "${timelog}" -- "${STAGE_DIR}/run_gquad.sh" "${input}" "${out}" > "${log}" 2>&1
            exit_code=$?                              # capture its exit status
            n_hits=0                                  # default to zero predictions
            # count non-comment GFF3 lines as predictions
            [ -s "${out}" ] && n_hits=$(grep -vc "^#" "${out}" 2>/dev/null || echo 0)
            ;;
        pqsfinder)
            # run pqsfinder G4 prediction under timing
            measure_and_run "${timelog}" -- "${STAGE_DIR}/run_pqsfinder.sh" "${input}" "${out}" > "${log}" 2>&1
            exit_code=$?                              # capture its exit status
            n_hits=0                                  # default to zero predictions
            # count non-comment GFF3 lines as predictions
            [ -s "${out}" ] && n_hits=$(grep -vc "^#" "${out}" 2>/dev/null || echo 0)
            ;;
    esac

    # set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file
    parse_time_metrics "${timelog}"
    # assume the run is invalid until it produces hits
    valid=0
    metric="0 hits"                                   # default human-readable metric
    if [ "${n_hits}" -gt 0 ] 2>/dev/null; then        # if the tool actually found something...
        valid=1                                       # ...mark the run valid...
        metric="${n_hits} hits"                       # ...and record the count
    fi
    # write this run's row to summary.tsv
    append_summary_row "motif_mapping_illumina" "${tool}" "${sample}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${exit_code}" "${valid}" "${metric}"

    # a positive control that finds nothing is suspect
    if [ "${sample}" = "positive_control" ] && [ "${valid}" -eq 0 ]; then
        # the severity differs between TFBS tools and G4 tools
        case "${tool}" in
            gquad|pqsfinder)
                # softer note: no known G4 to require here
                echo "NOTE: ${tool} found zero G-quadruplex predictions on HXB2's own U3 -- unlike the TFBS tools, there's no literature-confirmed G4 site here to treat as a strict positive control, but worth a manual look if this seems implausible." >&2
                ;;
            *)
                # hard warning: JASPAR motifs are expected here
                echo "WARNING: ${tool} found zero motif hits on HXB2's own U3 (positive control failed) -- treat its subset numbers as suspect until this is investigated." >&2
                ;;
        esac
    fi
}

# run every tool through both the positive control and the full subset
for TOOL in fimo moods tfbstools gquad pqsfinder; do
    echo "=== ${TOOL}: positive control (HXB2 own U3) ==="  # progress marker
    # each tool writes its positive-control output to a tool-specific path/format
    case "${TOOL}" in
        fimo) PC_OUT="${RESULTS_DIR}/fimo_positive_control_out" ;;       # FIMO writes a directory
        moods) PC_OUT="${RESULTS_DIR}/moods_positive_control.tsv" ;;     # MOODS writes a TSV
        # TFBSTools writes GFF3
        tfbstools) PC_OUT="${RESULTS_DIR}/tfbstools_positive_control.gff3" ;;
        gquad) PC_OUT="${RESULTS_DIR}/gquad_positive_control.gff3" ;;    # gquad writes GFF3
        # pqsfinder writes GFF3
        pqsfinder) PC_OUT="${RESULTS_DIR}/pqsfinder_positive_control.gff3" ;;
    esac
    # scan HXB2's own U3 first as a sanity check
    run_tool "${TOOL}" "positive_control" "${HXB2_U3_ONLY}" "${PC_OUT}"

    echo "=== ${TOOL}: full subset ==="              # progress marker
    # same per-tool output paths, now for the full U3 subset
    case "${TOOL}" in
        fimo) OUT="${RESULTS_DIR}/fimo_out" ;;        # FIMO writes a directory
        moods) OUT="${RESULTS_DIR}/moods_out.tsv" ;;  # MOODS writes a TSV
        tfbstools) OUT="${RESULTS_DIR}/tfbstools_out.gff3" ;;  # TFBSTools writes GFF3
        gquad) OUT="${RESULTS_DIR}/gquad_out.gff3" ;;  # gquad writes GFF3
        pqsfinder) OUT="${RESULTS_DIR}/pqsfinder_out.gff3" ;;  # pqsfinder writes GFF3
    esac
    # scan the full set of extracted U3 sequences
    run_tool "${TOOL}" "all_subset" "${U3_FASTA}" "${OUT}"
done

# final confirmation pointing the user at the results table
echo "Done. See ${SUMMARY_TSV}"
