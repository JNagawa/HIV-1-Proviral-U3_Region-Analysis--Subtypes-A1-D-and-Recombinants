#!/bin/bash
# Thin wrapper around run_pqsfinder.R (kept separate so
# compare_motif_mapping_illumina.sh's per-tool invocation pattern matches
# fimo/moods/tfbstools).
# Usage: run_pqsfinder.sh <SEQUENCES_FASTA> <OUT_GFF3>
# -u errors on unset vars, pipefail surfaces failures in a pipeline
set -uo pipefail
# arg 1 = input FASTA to scan; arg 2 = output GFF3 of G4 predictions
SEQS="$1" OUT="$2"
# absolute path of this script's dir, so it can find the R script
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
# delegate the actual pqsfinder G-quadruplex prediction to the R implementation
Rscript "${STAGE_DIR}/run_pqsfinder.R" "${SEQS}" "${OUT}"
