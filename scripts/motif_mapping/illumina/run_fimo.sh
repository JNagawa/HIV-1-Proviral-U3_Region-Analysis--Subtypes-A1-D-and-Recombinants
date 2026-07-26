#!/bin/bash
# FIMO (MEME Suite) -- PWM motif scanning against the 6 core-TF JASPAR set.
# Usage: run_fimo.sh <SEQUENCES_FASTA> <OUTDIR>
set -uo pipefail                                     # -u errors on unset vars, pipefail surfaces failures in a pipeline
SEQS="$1" OUTDIR="$2"                                 # arg 1 = input FASTA to scan; arg 2 = FIMO output dir
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"            # absolute path of this script's dir, so it runs from anywhere
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"      # repo top-level (three dirs up), base for the JASPAR path
fimo --oc "${OUTDIR}" --thresh 1e-4 "${REPO_ROOT}/data/reference/jaspar/core6_pfms.meme" "${SEQS}"  # scan seqs with the 6 MEME PWMs; --oc overwrites outdir; p-value threshold 1e-4
