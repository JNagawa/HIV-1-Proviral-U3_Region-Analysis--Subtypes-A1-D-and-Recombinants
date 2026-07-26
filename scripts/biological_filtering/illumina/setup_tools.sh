#!/bin/bash
# One-time clone + setup for the three biological-filtering candidates. All
# three are GitHub-hosted (not on conda/PyPI), confirmed by direct lookup:
#   - Poplars (Hypermut 3):   https://github.com/PoonLab/Poplars
#   - HIVSeqinR:              https://github.com/guineverelee/HIVSeqinR
#   - HIVIntact ("proviral"): https://github.com/ramics/HIVIntact
#     (the actual pip-installable package inside this repo is named
#     "intactness-pipeline" per its own README; the CLI it installs is
#     called `proviral`)
set -euo pipefail                                    # strict mode: -e aborts on any error since a broken clone/install shouldn't proceed
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"           # absolute path of this script's own dir
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"     # repo top-level (three dirs up)
mkdir -p "${REPO_ROOT}/scripts/tools"                # ensure the shared tools dir exists
cd "${REPO_ROOT}/scripts/tools"                      # clone everything into scripts/tools/

if [ ! -d Poplars ]; then                            # only clone/install Poplars if not already present (idempotent)
    git clone https://github.com/PoonLab/Poplars.git # fetch the Poplars source
    pip install -e ./Poplars                         # editable install so its deps resolve (run via hypermut.py directly)
fi

if [ ! -d HIVSeqinR ]; then                          # only clone HIVSeqinR if missing
    git clone https://github.com/guineverelee/HIVSeqinR.git  # fetch it; it's an R script run in-place, so no install step
fi

if [ ! -d HIVIntact ]; then                          # only set up HIVIntact if missing
    git clone --recurse-submodules https://github.com/ramics/HIVIntact.git  # clone with submodules (subtype alignments live in a submodule)
    python3 -m venv HIVIntact/env                    # dedicated venv so proviral's deps don't clash with the base env
    # Python 3.12's venv module no longer bundles setuptools, but the
    # installed `proviral` CLI imports pkg_resources (from setuptools) --
    # confirmed by the ModuleNotFoundError this throws without it.
    HIVIntact/env/bin/pip install setuptools         # add setuptools back so pkg_resources is importable at runtime
    HIVIntact/env/bin/pip install ./HIVIntact        # install the intactness-pipeline package, providing the `proviral` CLI
fi

echo "Setup complete. IMPORTANT caveats before running compare_biological_filtering_illumina.sh:"  # print post-setup warnings the user must read
echo "1. HIVSeqinR is NOT a CLI tool -- it's an R script meant to be run in"
echo "   RStudio ('highlight all, run all'). It has hardcoded config at the"
echo "   top of R_HIVSeqinR_Combined_ver04.R (blast DB path, primer"
echo "   sequences) that must be edited before it will run, and it REQUIRES"
echo "   input FASTA with no dashes or IUPAC ambiguity codes -- but bcftools"
echo "   consensus output legitimately contains IUPAC codes (R/Y/W/etc) at"
echo "   heterozygous positions. run_hivseqinr.sh strips these to the"
echo "   reference-matching base as a workaround; note this in your methods"
echo "   write-up as a real limitation, not a transparent substitution."
echo "2. HIVIntact requires a --subtype argument from a fixed list in"
echo "   util/subtype_alignments -- confirm A1/D are actually present before"
echo "   trusting its output on this cohort (the tools review already flags"
echo "   its heuristics as subtype-B-optimised)."
