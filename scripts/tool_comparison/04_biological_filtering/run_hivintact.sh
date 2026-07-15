#!/bin/bash
# HIVIntact ("proviral" CLI, confirmed from its own README).
# NOTE: its heuristics are optimised for subtype B (per the tools review);
# results on A1/D/recombinant sequences should be treated with caution.
# Confirm your target subtype is actually present in
# tools/HIVIntact/util/subtype_alignments before trusting output for it.
# Usage: run_hivintact.sh <INPUT_FASTA> <OUTDIR> [SUBTYPE]
set -uo pipefail
IN="$1" OUTDIR="$2" SUBTYPE="${3:-B}"

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
HIVINTACT_ENV="${STAGE_DIR}/tools/HIVIntact/env/bin/activate"
[ -f "${HIVINTACT_ENV}" ] || { echo "ERROR: HIVIntact venv not found, run setup_tools.sh first." >&2; exit 1; }

mkdir -p "${OUTDIR}"
cd "${OUTDIR}" || exit 1
# shellcheck disable=SC1090
source "${HIVINTACT_ENV}"
proviral intact --subtype "${SUBTYPE}" "${IN}" > "${OUTDIR}/hivintact.log" 2>&1
deactivate
