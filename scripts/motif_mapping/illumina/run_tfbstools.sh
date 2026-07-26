#!/bin/bash
# Thin wrapper around run_tfbstools.R (kept as a separate .sh so
# compare_motif_mapping_illumina.sh's per-tool invocation pattern matches fimo/moods).
# Usage: run_tfbstools.sh <SEQUENCES_FASTA> <OUT_GFF3>
set -uo pipefail                                     # -u errors on unset vars, pipefail surfaces failures in a pipeline
SEQS="$1" OUT="$2"                                    # arg 1 = input FASTA to scan; arg 2 = output GFF3 of hits
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"            # absolute path of this script's dir, so it can find the R script
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"      # repo top-level (three dirs up), base for the JASPAR path
Rscript "${STAGE_DIR}/run_tfbstools.R" "${SEQS}" "${REPO_ROOT}/data/reference/jaspar/core6_pfms.jaspar" "${OUT}"  # delegate to R; pass seqs, the 6-TF JASPAR flat file, and the output path
