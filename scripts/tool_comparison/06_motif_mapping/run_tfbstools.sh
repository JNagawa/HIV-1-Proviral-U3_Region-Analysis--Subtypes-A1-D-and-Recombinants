#!/bin/bash
# Thin wrapper around run_tfbstools.R (kept as a separate .sh so
# compare_06.sh's per-tool invocation pattern matches fimo/moods).
# Usage: run_tfbstools.sh <SEQUENCES_FASTA> <OUT_GFF3>
set -uo pipefail
SEQS="$1" OUT="$2"
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
Rscript "${STAGE_DIR}/run_tfbstools.R" "${SEQS}" "${STAGE_DIR}/jaspar/core6_pfms.jaspar" "${OUT}"
