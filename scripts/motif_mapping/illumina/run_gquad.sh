#!/bin/bash
# Thin wrapper around run_gquad.R (kept separate so
# compare_motif_mapping_illumina.sh's per-tool invocation pattern matches
# fimo/moods/tfbstools/pqsfinder).
# Usage: run_gquad.sh <SEQUENCES_FASTA> <OUT_GFF3>
set -uo pipefail
SEQS="$1" OUT="$2"
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
Rscript "${STAGE_DIR}/run_gquad.R" "${SEQS}" "${OUT}"
