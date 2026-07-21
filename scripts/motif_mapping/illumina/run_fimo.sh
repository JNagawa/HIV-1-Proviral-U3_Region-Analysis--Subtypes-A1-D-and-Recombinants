#!/bin/bash
# FIMO (MEME Suite) -- PWM motif scanning against the 6 core-TF JASPAR set.
# Usage: run_fimo.sh <SEQUENCES_FASTA> <OUTDIR>
set -uo pipefail
SEQS="$1" OUTDIR="$2"
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
fimo --oc "${OUTDIR}" --thresh 1e-4 "${REPO_ROOT}/data/reference/jaspar/core6_pfms.meme" "${SEQS}"
