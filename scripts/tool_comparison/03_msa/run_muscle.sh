#!/bin/bash
# Usage: run_muscle.sh <COMBINED_FASTA> <OUT_FASTA> <LOG>
set -uo pipefail
IN="$1" OUT="$2" LOG="$3"
muscle -align "${IN}" -output "${OUT}" -threads "${THREADS:-4}" > "${LOG}" 2>&1
