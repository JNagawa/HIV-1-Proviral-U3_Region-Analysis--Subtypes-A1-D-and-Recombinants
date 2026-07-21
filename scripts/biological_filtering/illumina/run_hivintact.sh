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
set -uo pipefail
IN="$1" OUTDIR="$2" SUBTYPE="${3:-A1}"

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
HIVINTACT_ENV="${REPO_ROOT}/scripts/tools/HIVIntact/env/bin/activate"
[ -f "${HIVINTACT_ENV}" ] || { echo "ERROR: HIVIntact venv not found, run setup_tools.sh first." >&2; exit 1; }

mkdir -p "${OUTDIR}"
LOG="${OUTDIR}/hivintact.log"
# proviral writes intact.fasta/nonintact.fasta/orfs.json/errors.json relative
# to the CWD, not --working-folder (that flag is for something else --
# confirmed by testing: without this cd, output landed in the caller's cwd).
cd "${OUTDIR}" || exit 1
# shellcheck disable=SC1090
source "${HIVINTACT_ENV}"
proviral intact --subtype "${SUBTYPE}" --working-folder "${OUTDIR}" "${IN}" > "${LOG}" 2>&1
deactivate
