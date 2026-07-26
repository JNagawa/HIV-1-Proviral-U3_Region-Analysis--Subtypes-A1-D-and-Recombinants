#!/bin/bash
# Usage: run_mafft.sh <COMBINED_FASTA> <OUT_FASTA> <LOG>
set -uo pipefail                                     # -u errors on unset vars, pipefail fails a pipe if any stage fails
IN="$1" OUT="$2" LOG="$3"                            # positional args: input FASTA, output alignment, log file
mafft --localpair --maxiterate 1000 --thread "${THREADS:-4}" "${IN}" > "${OUT}" 2> "${LOG}"  # MAFFT L-INS-i (accurate local-pair mode); alignment to stdout, diagnostics to the log
