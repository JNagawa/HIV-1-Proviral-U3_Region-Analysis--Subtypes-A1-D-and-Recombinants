#!/bin/bash
# PacBio arm: compares FIMO vs MOODS vs TFBSTools (PWM TFBS scanning, 6 core
# JASPAR TFs) and, separately, gquad vs pqsfinder (G-quadruplex prediction),
# against U3 sequences extracted from the PacBio msa stage's MAFFT alignment.
# Reuses the Illumina arm's run_<tool>.sh scripts and the shared U3-extraction
# util unchanged (motif scanning is platform-agnostic; only the aligned input
# differs). Each tool runs a positive control (HXB2's own U3) before the
# subset, same as the Illumina harness.
# Usage: ./compare_motif_mapping_pacbio.sh   (via sbatch scripts/utils/run_comparison_step.slurm.sh)
set -uo pipefail                                     # -u errors on unset vars, pipefail catches pipeline failures; no -e so one tool failing doesn't abort the comparison

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"            # absolute path of this script's dir, so it runs from anywhere
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"      # repo top-level (three dirs up), base for every path below
source "${REPO_ROOT}/scripts/common/lib_compare.sh"  # load shared helpers: measure_and_run, parse_time_metrics, append_summary_row
RUN_DIR="${REPO_ROOT}/scripts/motif_mapping/illumina"   # reuse shared run_<tool>.sh

ALIGNMENT="${REPO_ROOT}/results/msa/pacbio/mafft_aligned.fasta"  # input: MAFFT alignment from the PacBio msa stage
RESULTS_DIR="${REPO_ROOT}/results/motif_mapping/pacbio"          # output: all motif-mapping results + summary go here
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"             # the single TSV every tool appends a timing/validity row to
mkdir -p "${RESULTS_DIR}"                            # create the results dir (and parents) if needed
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"  # seed the manual notes file from the template on first run only

if [ ! -s "${ALIGNMENT}" ]; then                     # can't map motifs without the upstream alignment...
    echo "ERROR: no alignment at ${ALIGNMENT} (run scripts/msa/pacbio/compare_msa_pacbio.sh first)." >&2  # ...tell the user which stage to run first...
    exit 1                                           # ...and abort
fi

echo "=== Extracting U3 (anchored to HXB2's own annotated LTR/R-region) ==="  # progress marker in the log
U3_FASTA="${RESULTS_DIR}/U3_extracted.fasta"         # ungapped U3 sequences, the actual scanning input for every tool
bash "${REPO_ROOT}/scripts/utils/extract_u3_by_hxb2_anchor.sh" \
    --alignment "${ALIGNMENT}" \
    --hxb2-id K03455.1 \
    --gb-cache "${REPO_ROOT}/data/reference/K03455.1.gb" \
    --out-gapped "${RESULTS_DIR}/U3_aligned.fasta" \
    --out "${U3_FASTA}" \
    --warnings-log "${RESULTS_DIR}/u3_extraction_warnings.log"  # slice U3 out of each aligned seq using HXB2's annotated LTR/R boundary as the anchor
if [ ! -s "${U3_FASTA}" ]; then                      # if extraction yielded nothing there's nothing to scan...
    echo "ERROR: U3 extraction produced no sequences, see extraction output above." >&2  # ...report it...
    exit 1                                           # ...and abort
fi

HXB2_U3_ONLY="${RESULTS_DIR}/hxb2_u3_only.fasta"     # HXB2's own U3 alone, used as each tool's positive control
seqkit grep -n -r -p "^K03455" "${U3_FASTA}" > "${HXB2_U3_ONLY}"  # pull just the HXB2 record (ID starting K03455) out of the U3 set
if [ ! -s "${HXB2_U3_ONLY}" ]; then                  # if HXB2's U3 couldn't be isolated, the positive control is impossible...
    echo "ERROR: could not isolate HXB2's own U3 from ${U3_FASTA} for the positive control." >&2  # ...report it...
    exit 1                                           # ...and abort
fi

run_tool() {                                         # run one tool on one input, time it, count hits, and log a summary row
    local tool="$1" sample="$2" input="$3" out="$4"  # args: tool name, sample label, input FASTA, output path
    local timelog="${RESULTS_DIR}/${tool}_${sample}.time"  # per-run timing file for measure_and_run
    local log="${RESULTS_DIR}/${tool}_${sample}.log"       # per-run captured stdout+stderr
    local exit_code n_hits valid metric              # locals set per tool below

    case "${tool}" in                                # dispatch to the right wrapper and hit-counting rule per tool
        fimo)
            measure_and_run "${timelog}" -- "${RUN_DIR}/run_fimo.sh" "${input}" "${out}" > "${log}" 2>&1  # run the shared FIMO wrapper under timing
            exit_code=$?; n_hits=0                    # capture exit status; default to zero hits
            [ -s "${out}/fimo.tsv" ] && n_hits=$(grep -vc "^#\|^motif_id" "${out}/fimo.tsv" 2>/dev/null || echo 0) ;;  # count data rows in fimo.tsv, excluding comments and the header
        moods)
            measure_and_run "${timelog}" -- "${RUN_DIR}/run_moods.sh" "${input}" "${out}" > "${log}" 2>&1  # run the shared MOODS wrapper under timing
            exit_code=$?; n_hits=0                    # capture exit status; default to zero hits
            [ -s "${out}" ] && n_hits=$(wc -l < "${out}") ;;  # MOODS output has one hit per line, so line count = hit count
        tfbstools)
            measure_and_run "${timelog}" -- "${RUN_DIR}/run_tfbstools.sh" "${input}" "${out}" > "${log}" 2>&1  # run the shared TFBSTools wrapper under timing
            exit_code=$?; n_hits=0                    # capture exit status; default to zero hits
            [ -s "${out}" ] && n_hits=$(grep -vc "^#" "${out}" 2>/dev/null || echo 0) ;;  # count non-comment GFF3 lines as hits
        gquad)
            measure_and_run "${timelog}" -- "${RUN_DIR}/run_gquad.sh" "${input}" "${out}" > "${log}" 2>&1  # run the shared gquad G4-prediction wrapper under timing
            exit_code=$?; n_hits=0                    # capture exit status; default to zero predictions
            [ -s "${out}" ] && n_hits=$(grep -vc "^#" "${out}" 2>/dev/null || echo 0) ;;  # count non-comment GFF3 lines as predictions
        pqsfinder)
            measure_and_run "${timelog}" -- "${RUN_DIR}/run_pqsfinder.sh" "${input}" "${out}" > "${log}" 2>&1  # run the shared pqsfinder G4-prediction wrapper under timing
            exit_code=$?; n_hits=0                    # capture exit status; default to zero predictions
            [ -s "${out}" ] && n_hits=$(grep -vc "^#" "${out}" 2>/dev/null || echo 0) ;;  # count non-comment GFF3 lines as predictions
    esac

    parse_time_metrics "${timelog}"                  # set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file
    valid=0; metric="0 hits"                         # assume invalid with zero hits until proven otherwise
    if [ "${n_hits}" -gt 0 ] 2>/dev/null; then valid=1; metric="${n_hits} hits"; fi  # if the tool found something, mark valid and record the count
    append_summary_row "motif_mapping_pacbio" "${tool}" "${sample}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${exit_code}" "${valid}" "${metric}"  # write this run's row to summary.tsv

    if [ "${sample}" = "positive_control" ] && [ "${valid}" -eq 0 ]; then  # a positive control that finds nothing is suspect
        case "${tool}" in                             # severity differs between G4 tools and TFBS tools
            gquad|pqsfinder)
                echo "NOTE: ${tool} found zero G4 predictions on HXB2's own U3 -- no literature-confirmed G4 site here to treat as a strict positive control, but worth a manual look if implausible." >&2 ;;  # softer note: no known G4 to require here
            *)
                echo "WARNING: ${tool} found zero motif hits on HXB2's own U3 (positive control failed) -- treat its subset numbers as suspect until investigated." >&2 ;;  # hard warning: JASPAR motifs are expected here
        esac
    fi
}

for TOOL in fimo moods tfbstools gquad pqsfinder; do  # run every tool through both the positive control and the full subset
    echo "=== ${TOOL}: positive control (HXB2 own U3) ==="  # progress marker
    case "${TOOL}" in                                 # each tool writes its positive-control output to a tool-specific path/format
        fimo) PC_OUT="${RESULTS_DIR}/fimo_positive_control_out" ;;       # FIMO writes a directory
        moods) PC_OUT="${RESULTS_DIR}/moods_positive_control.tsv" ;;     # MOODS writes a TSV
        tfbstools) PC_OUT="${RESULTS_DIR}/tfbstools_positive_control.gff3" ;;  # TFBSTools writes GFF3
        gquad) PC_OUT="${RESULTS_DIR}/gquad_positive_control.gff3" ;;    # gquad writes GFF3
        pqsfinder) PC_OUT="${RESULTS_DIR}/pqsfinder_positive_control.gff3" ;;  # pqsfinder writes GFF3
    esac
    run_tool "${TOOL}" "positive_control" "${HXB2_U3_ONLY}" "${PC_OUT}"  # scan HXB2's own U3 first as a sanity check

    echo "=== ${TOOL}: full subset ==="              # progress marker
    case "${TOOL}" in                                 # same per-tool output paths, now for the full U3 subset
        fimo) OUT="${RESULTS_DIR}/fimo_out" ;;        # FIMO writes a directory
        moods) OUT="${RESULTS_DIR}/moods_out.tsv" ;;  # MOODS writes a TSV
        tfbstools) OUT="${RESULTS_DIR}/tfbstools_out.gff3" ;;  # TFBSTools writes GFF3
        gquad) OUT="${RESULTS_DIR}/gquad_out.gff3" ;;  # gquad writes GFF3
        pqsfinder) OUT="${RESULTS_DIR}/pqsfinder_out.gff3" ;;  # pqsfinder writes GFF3
    esac
    run_tool "${TOOL}" "all_subset" "${U3_FASTA}" "${OUT}"  # scan the full set of extracted U3 sequences
done

echo "Done. See ${SUMMARY_TSV}"                      # final confirmation pointing the user at the results table
