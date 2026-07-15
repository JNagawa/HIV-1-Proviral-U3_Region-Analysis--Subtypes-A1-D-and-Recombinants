#!/bin/bash
# jpHMM -- recombination/subtype caller. CLI confirmed directly from
# src/main.cpp's getopt() call (not guessed):
#   -s <query_fasta> -v <virus_type> -I <input_dir> -P <priors_dir> -o <output_dir>
# Uses jpHMM's own bundled HIV reference alignment/priors (tools/jpHMM/input,
# tools/jpHMM/priors) -- no separate reference sourcing needed, unlike SHIVER.
# Usage: run_jphmm.sh <QUERY_FASTA> <OUTDIR>
set -uo pipefail
QUERY="$1" OUTDIR="$2"

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
JPHMM_DIR="${STAGE_DIR}/tools/jpHMM"
JPHMM_BIN="${JPHMM_DIR}/src/jpHMM"

[ -x "${JPHMM_BIN}" ] || { echo "ERROR: jpHMM not built, run setup_jphmm.sh first." >&2; exit 1; }

mkdir -p "${OUTDIR}"
"${JPHMM_BIN}" \
    -s "${QUERY}" \
    -v HIV \
    -I "${JPHMM_DIR}/input" \
    -P "${JPHMM_DIR}/priors" \
    -o "${OUTDIR}"
