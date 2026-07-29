#!/bin/bash
# Usage: run_clustalo.sh <COMBINED_FASTA> <OUT_FASTA> <LOG>
# -u errors on unset vars, pipefail fails a pipe if any stage fails
set -uo pipefail
# positional args: input FASTA, output alignment, log file
IN="$1" OUT="$2" LOG="$3"
# run Clustal Omega; --force overwrites a stale output; redirect stdout+stderr to the log
clustalo -i "${IN}" -o "${OUT}" --threads "${THREADS:-4}" --force > "${LOG}" 2>&1
