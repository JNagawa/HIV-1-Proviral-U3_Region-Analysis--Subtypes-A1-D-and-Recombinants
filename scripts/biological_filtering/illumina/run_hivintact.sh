#!/bin/bash
# HIVIntact ("proviral" CLI). Real flags confirmed from `proviral intact
# --help` (the tools review/README guessed `-o` for output dir; the actual
# flag is `--working-folder`).
# NOTE: its heuristics are optimised for subtype B (per the tools review).
# This cohort (PRJNA207834) is A1/D/recombinant, not B -- A1 and D are both
# present in scripts/tools/HIVIntact/util/subtype_alignments, so default to
# A1 rather than B; per-sample subtype (e.g. from jpHMM's own call) should
# override this once available, since HIVIntact's heuristics are still
# validated on B and treating A1 as a drop-in match is an approximation.
# Usage: run_hivintact.sh <INPUT_FASTA> <OUTDIR> [SUBTYPE]
# -u errors on unset vars, pipefail fails a pipe if any stage fails
set -uo pipefail
# positional args: input FASTA, output dir, subtype (defaults to A1 for this cohort)
IN="$1" OUTDIR="$2" SUBTYPE="${3:-A1}"

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"           # absolute path of this script's own dir
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"     # repo top-level (three dirs up)
# the dedicated venv setup_tools.sh created for the proviral CLI
HIVINTACT_ENV="${REPO_ROOT}/scripts/tools/HIVIntact/env/bin/activate"
# bail early with guidance if the venv is missing
[ -f "${HIVINTACT_ENV}" ] || { echo "ERROR: HIVIntact venv not found, run setup_tools.sh first." >&2; exit 1; }

mkdir -p "${OUTDIR}"                                 # ensure the output dir exists
LOG="${OUTDIR}/hivintact.log"                        # captured stdout+stderr of the proviral run
# proviral writes intact.fasta/nonintact.fasta/orfs.json/errors.json relative
# to the CWD, not --working-folder (that flag is for something else --
# confirmed by testing: without this cd, output landed in the caller's cwd).
# move into OUTDIR so proviral's CWD-relative outputs land there
cd "${OUTDIR}" || exit 1
# shellcheck disable=SC1090
# activate the proviral venv so the CLI and its deps are on PATH
source "${HIVINTACT_ENV}"
# run the intactness classifier, capturing all output to the log
proviral intact --subtype "${SUBTYPE}" --working-folder "${OUTDIR}" "${IN}" > "${LOG}" 2>&1
deactivate                                            # leave the venv cleanly
