#!/bin/bash
# Usage: run_clustalo.sh <COMBINED_FASTA> <OUT_FASTA> <LOG>
set -uo pipefail
IN="$1" OUT="$2" LOG="$3"
clustalo -i "${IN}" -o "${OUT}" --threads "${THREADS:-4}" --force > "${LOG}" 2>&1
