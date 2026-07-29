#!/bin/bash
# Poplars Hypermut 3 -- APOBEC3G/F G-to-A hypermutation screen.
# Poplars has no __main__.py / console_scripts entry point despite being pip
# installed -- `python -m poplars hypermut` fails with "'poplars' is a
# package and cannot be directly executed". The real invocation is running
# poplars/hypermut.py directly as a script:
#   python <path-to>/poplars/hypermut.py <fasta> [--consensus] [--skip N] [--out FILE]
# With no reference specified, the FIRST sequence in the FASTA is used as
# the reference -- so INPUT_FASTA must have HXB2 as its first record.
# Usage: run_poplars.sh <INPUT_FASTA_HXB2_FIRST> <OUT_TSV>
# -u errors on unset vars, pipefail fails a pipe if any stage fails
set -uo pipefail
# positional args: input FASTA (HXB2 first) and output TSV path
IN="$1" OUT="$2"
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"           # absolute path of this script's own dir
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"     # repo top-level (three dirs up)
# run hypermut.py directly (no working entry point); first seq is the reference
python3 "${REPO_ROOT}/scripts/tools/Poplars/poplars/hypermut.py" "${IN}" --out "${OUT}"
