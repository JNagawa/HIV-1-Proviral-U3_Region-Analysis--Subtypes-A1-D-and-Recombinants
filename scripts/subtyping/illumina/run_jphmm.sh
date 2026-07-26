#!/bin/bash
# jpHMM -- recombination/subtype caller. CLI confirmed directly from
# src/main.cpp's getopt() call (not guessed):
#   -s <query_fasta> -v <virus_type> -I <input_dir> -P <priors_dir> -o <output_dir>
# Uses jpHMM's own bundled HIV reference alignment/priors
# (scripts/tools/jpHMM/input, scripts/tools/jpHMM/priors) -- no separate
# reference sourcing needed, unlike SHIVER.
# Usage: run_jphmm.sh <QUERY_FASTA> <OUTDIR>
set -uo pipefail                                     # -u errors on unset vars, pipefail catches mid-pipe failures
QUERY="$1" OUTDIR="$2"                               # args: query consensus FASTA to subtype, and output dir

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"           # absolute path of this script's own dir, so paths resolve regardless of launch location
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"     # repo top-level (three dirs up), the base for the tool paths below
JPHMM_DIR="${REPO_ROOT}/scripts/tools/jpHMM"         # where setup_jphmm.sh installed jpHMM (source, input, priors)
JPHMM_BIN="${JPHMM_DIR}/src/jpHMM"                   # the compiled jpHMM binary

[ -x "${JPHMM_BIN}" ] || { echo "ERROR: jpHMM not built, run setup_jphmm.sh first." >&2; exit 1; }  # bail early if the binary is missing/unbuilt

# This HPC's Lmod-managed LD_LIBRARY_PATH puts /opt/ohpc/pub/apps/anaconda3/lib
# (an older libstdc++) ahead of the active conda env's own lib dir, so jpHMM
# (compiled with a newer g++) fails to dynamically link against GLIBCXX
# symbols that are actually present in the env's libstdc++ -- just shadowed.
# Prepend CONDA_PREFIX/lib so it's found first, scoped to this invocation only.
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"  # put the env's newer libstdc++ ahead of the HPC's older one for this run

mkdir -p "${OUTDIR}"                                 # ensure the output dir exists before jpHMM writes into it
"${JPHMM_BIN}" \
    -s "${QUERY}" \
    -v HIV \
    -I "${JPHMM_DIR}/input" \
    -P "${JPHMM_DIR}/priors" \
    -o "${OUTDIR}"                                   # run jpHMM: subtype the query against its bundled HIV alignment/priors
