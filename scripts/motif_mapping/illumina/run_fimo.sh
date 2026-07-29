#!/bin/bash
# FIMO (MEME Suite) -- PWM motif scanning against the 6 core-TF JASPAR set.
# Usage: run_fimo.sh <SEQUENCES_FASTA> <OUTDIR>
# -u errors on unset vars, pipefail surfaces failures in a pipeline
set -uo pipefail
# arg 1 = input FASTA to scan; arg 2 = FIMO output dir
SEQS="$1" OUTDIR="$2"
# absolute path of this script's dir, so it runs from anywhere
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
# repo top-level (three dirs up), base for the JASPAR path
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
# scan seqs with the 6 MEME PWMs; --oc overwrites outdir; p-value threshold 1e-4
fimo --oc "${OUTDIR}" --thresh 1e-4 "${REPO_ROOT}/data/reference/jaspar/core6_pfms.meme" "${SEQS}"
