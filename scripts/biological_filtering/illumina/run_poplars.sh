#!/bin/bash
# Poplars Hypermut 3 -- APOBEC3G/F G-to-A hypermutation screen.
# Real CLI (confirmed from poplars/hypermut.py source, not guessed):
#   python -m poplars hypermut <fasta> [--consensus] [--skip N] [--out FILE]
# With no reference specified, the FIRST sequence in the FASTA is used as
# the reference -- so INPUT_FASTA must have HXB2 as its first record.
# Usage: run_poplars.sh <INPUT_FASTA_HXB2_FIRST> <OUT_TSV>
set -uo pipefail
IN="$1" OUT="$2"
python3 -m poplars hypermut "${IN}" --out "${OUT}"
