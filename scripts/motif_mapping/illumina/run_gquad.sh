#!/bin/bash
# Thin wrapper around run_gquad.R (kept separate so
# compare_motif_mapping_illumina.sh's per-tool invocation pattern matches
# fimo/moods/tfbstools/pqsfinder).
# Usage: run_gquad.sh <SEQUENCES_FASTA> <OUT_GFF3>
set -uo pipefail                                     # -u errors on unset vars, pipefail surfaces failures in a pipeline
SEQS="$1" OUT="$2"                                    # arg 1 = input FASTA to scan; arg 2 = output GFF3 of G4 predictions
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"            # absolute path of this script's dir, so it can find the R script
Rscript "${STAGE_DIR}/run_gquad.R" "${SEQS}" "${OUT}"  # delegate the actual G-quadruplex prediction to the R implementation
