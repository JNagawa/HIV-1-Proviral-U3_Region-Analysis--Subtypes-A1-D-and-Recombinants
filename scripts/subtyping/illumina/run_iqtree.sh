#!/bin/bash
# IQ-TREE 2 -- maximum-likelihood phylogenetic confirmation of subtype
# assignment. Not a subtype *caller* by itself the way jpHMM is; used here
# per the review's recommendation to confirm subtype calls with a
# bootstrap-supported ML tree (-m MFP auto-selects the substitution model,
# -B 1000 runs ultrafast bootstrap).
# Usage: run_iqtree.sh <ALIGNMENT_FASTA> <OUTDIR> <PREFIX>
set -uo pipefail
ALN="$1" OUTDIR="$2" PREFIX="$3"
mkdir -p "${OUTDIR}"
iqtree2 -s "${ALN}" -m MFP -B 1000 -T "${THREADS:-4}" --prefix "${OUTDIR}/${PREFIX}" -redo
