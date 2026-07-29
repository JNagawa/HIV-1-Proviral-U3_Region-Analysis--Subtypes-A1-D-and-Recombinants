#!/bin/bash
# jpHMM -- recombination/subtype caller. CLI confirmed directly from
# src/main.cpp's getopt() call (not guessed):
#   -s <query_fasta> -v <virus_type> -I <input_dir> -P <priors_dir> -o <output_dir>
# Uses jpHMM's own bundled HIV reference alignment/priors
# (scripts/tools/jpHMM/input, scripts/tools/jpHMM/priors) -- no separate
# reference sourcing needed, unlike SHIVER.
# Usage: run_jphmm.sh <QUERY_FASTA> <OUTDIR>
# -u errors on unset vars, pipefail catches mid-pipe failures
set -uo pipefail
# args: query consensus FASTA to subtype, and output dir
QUERY="$1" OUTDIR="$2"

# absolute path of this script's own dir, so paths resolve regardless of launch location
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
# repo top-level (three dirs up), the base for the tool paths below
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
# where setup_jphmm.sh installed jpHMM (source, input, priors)
JPHMM_DIR="${REPO_ROOT}/scripts/tools/jpHMM"
JPHMM_BIN="${JPHMM_DIR}/src/jpHMM"                   # the compiled jpHMM binary

# bail early if the binary is missing/unbuilt
[ -x "${JPHMM_BIN}" ] || { echo "ERROR: jpHMM not built, run setup_jphmm.sh first." >&2; exit 1; }

# This HPC's Lmod-managed LD_LIBRARY_PATH puts /opt/ohpc/pub/apps/anaconda3/lib
# (an older libstdc++) ahead of the active conda env's own lib dir, so jpHMM
# (compiled with a newer g++) fails to dynamically link against GLIBCXX
# symbols that are actually present in the env's libstdc++ -- just shadowed.
# Prepend CONDA_PREFIX/lib so it's found first, scoped to this invocation only.
# put the env's newer libstdc++ ahead of the HPC's older one for this run
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# ensure the output dir exists before jpHMM writes into it
mkdir -p "${OUTDIR}"
# run jpHMM: subtype the query against its bundled HIV alignment/priors
"${JPHMM_BIN}" \
    -s "${QUERY}" \
    -v HIV \
    -I "${JPHMM_DIR}/input" \
    -P "${JPHMM_DIR}/priors" \
    -o "${OUTDIR}"
