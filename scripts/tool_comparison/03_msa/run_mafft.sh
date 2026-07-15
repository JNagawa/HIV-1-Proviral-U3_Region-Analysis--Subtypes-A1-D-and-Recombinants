#!/bin/bash
# Usage: run_mafft.sh <COMBINED_FASTA> <OUT_FASTA> <LOG>
set -uo pipefail
IN="$1" OUT="$2" LOG="$3"
mafft --localpair --maxiterate 1000 --thread "${THREADS:-4}" "${IN}" > "${OUT}" 2> "${LOG}"
