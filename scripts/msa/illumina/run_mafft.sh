#!/bin/bash
# Usage: run_mafft.sh <COMBINED_FASTA> <OUT_FASTA> <LOG>
# -u errors on unset vars, pipefail fails a pipe if any stage fails
set -uo pipefail
# positional args: input FASTA, output alignment, log file
IN="$1" OUT="$2" LOG="$3"
# MAFFT L-INS-i (accurate local-pair mode); alignment to stdout, diagnostics to the log
mafft --localpair --maxiterate 1000 --thread "${THREADS:-4}" "${IN}" > "${OUT}" 2> "${LOG}"
