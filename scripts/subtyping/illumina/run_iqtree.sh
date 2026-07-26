#!/bin/bash
# IQ-TREE 2 -- maximum-likelihood phylogenetic confirmation of subtype
# assignment. Not a subtype *caller* by itself the way jpHMM is; used here
# per the review's recommendation to confirm subtype calls with a
# bootstrap-supported ML tree (-m MFP auto-selects the substitution model).
# No -B/ultrafast-bootstrap here: IQ-TREE2 refuses UFBoot when the alignment
# has too few effectively-distinct sequences, which the small (4-sample)
# comparison-harness subset triggers. Re-add -B 1000 once this is run
# against the full cohort, where bootstrap support is meaningful.
# Usage: run_iqtree.sh <ALIGNMENT_FASTA> <OUTDIR> <PREFIX>
set -uo pipefail                                     # -u errors on unset vars, pipefail catches mid-pipe failures
ALN="$1" OUTDIR="$2" PREFIX="$3"                     # args: input alignment, output dir, and output filename prefix
mkdir -p "${OUTDIR}"                                 # ensure the output dir exists before IQ-TREE2 writes into it
iqtree2 -s "${ALN}" -m MFP -T "${THREADS:-4}" --prefix "${OUTDIR}/${PREFIX}" -redo  # build ML tree; MFP auto-picks the model, -redo overwrites any prior run
