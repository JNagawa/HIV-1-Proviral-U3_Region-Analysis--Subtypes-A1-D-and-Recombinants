#!/bin/bash
# Pinned to classic MUSCLE 3.8 (see HIV_U3analysis_env.yml for why: 5.3's
# bioconda build crashes with SIGILL on this node's CPU). v3's CLI is
# -in/-out, single-threaded (no -threads flag), unlike v5's -align/-output.
# Usage: run_muscle.sh <COMBINED_FASTA> <OUT_FASTA> <LOG>
# -u errors on unset vars, pipefail fails a pipe if any stage fails
set -uo pipefail
# positional args: input FASTA, output alignment, log file
IN="$1" OUT="$2" LOG="$3"
# run MUSCLE 3.8 with its -in/-out CLI; redirect stdout+stderr to the log
muscle -in "${IN}" -out "${OUT}" > "${LOG}" 2>&1
