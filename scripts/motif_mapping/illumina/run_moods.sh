#!/bin/bash
# MOODS -- PWM motif scanning via its moods-dna.py CLI (confirmed usage from
# the MOODS GitHub wiki, not guessed):
#   moods-dna.py -m matrices*.pfm -s sequences.fa -p pvalue
# Usage: run_moods.sh <SEQUENCES_FASTA> <OUT_TSV>
set -uo pipefail                                     # -u errors on unset vars, pipefail surfaces failures in a pipeline
SEQS="$1" OUT="$2"                                    # arg 1 = input FASTA to scan; arg 2 = output TSV of hits
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"            # absolute path of this script's dir, so it runs from anywhere
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"      # repo top-level (three dirs up), base for the JASPAR path
moods-dna.py -m "${REPO_ROOT}"/data/reference/jaspar/MA*.pfm -s "${SEQS}" -p 0.0001 > "${OUT}"  # scan with the per-matrix .pfm files at p<=1e-4; MOODS wants one matrix per file
