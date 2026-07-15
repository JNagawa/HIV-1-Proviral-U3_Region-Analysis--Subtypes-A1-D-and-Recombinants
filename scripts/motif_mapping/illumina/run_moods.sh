#!/bin/bash
# MOODS -- PWM motif scanning via its moods-dna.py CLI (confirmed usage from
# the MOODS GitHub wiki, not guessed):
#   moods-dna.py -m matrices*.pfm -s sequences.fa -p pvalue
# Usage: run_moods.sh <SEQUENCES_FASTA> <OUT_TSV>
set -uo pipefail
SEQS="$1" OUT="$2"
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../.." && pwd)"
moods-dna.py -m "${REPO_ROOT}"/data/reference/jaspar/MA*.pfm -s "${SEQS}" -p 0.0001 > "${OUT}"
