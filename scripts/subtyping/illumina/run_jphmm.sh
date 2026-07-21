#!/bin/bash
# jpHMM -- recombination/subtype caller. CLI confirmed directly from
# src/main.cpp's getopt() call (not guessed):
#   -s <query_fasta> -v <virus_type> -I <input_dir> -P <priors_dir> -o <output_dir>
# Uses jpHMM's own bundled HIV reference alignment/priors
# (scripts/tools/jpHMM/input, scripts/tools/jpHMM/priors) -- no separate
# reference sourcing needed, unlike SHIVER.
# Usage: run_jphmm.sh <QUERY_FASTA> <OUTDIR>
set -uo pipefail
QUERY="$1" OUTDIR="$2"

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
JPHMM_DIR="${REPO_ROOT}/scripts/tools/jpHMM"
JPHMM_BIN="${JPHMM_DIR}/src/jpHMM"

[ -x "${JPHMM_BIN}" ] || { echo "ERROR: jpHMM not built, run setup_jphmm.sh first." >&2; exit 1; }

# This HPC's Lmod-managed LD_LIBRARY_PATH puts /opt/ohpc/pub/apps/anaconda3/lib
# (an older libstdc++) ahead of the active conda env's own lib dir, so jpHMM
# (compiled with a newer g++) fails to dynamically link against GLIBCXX
# symbols that are actually present in the env's libstdc++ -- just shadowed.
# Prepend CONDA_PREFIX/lib so it's found first, scoped to this invocation only.
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

mkdir -p "${OUTDIR}"
"${JPHMM_BIN}" \
    -s "${QUERY}" \
    -v HIV \
    -I "${JPHMM_DIR}/input" \
    -P "${JPHMM_DIR}/priors" \
    -o "${OUTDIR}"
