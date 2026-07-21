#!/bin/bash
# Thin wrapper around run_tfbstools.R (kept as a separate .sh so
# compare_motif_mapping_illumina.sh's per-tool invocation pattern matches fimo/moods).
# Usage: run_tfbstools.sh <SEQUENCES_FASTA> <OUT_GFF3>
set -uo pipefail
SEQS="$1" OUT="$2"
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
Rscript "${STAGE_DIR}/run_tfbstools.R" "${SEQS}" "${REPO_ROOT}/data/reference/jaspar/core6_pfms.jaspar" "${OUT}"
