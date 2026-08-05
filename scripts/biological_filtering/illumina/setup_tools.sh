#!/bin/bash
# One-time clone + setup for the three biological-filtering candidates. All
# three are GitHub-hosted (not on conda/PyPI), confirmed by direct lookup:
#   - Poplars (Hypermut 3):   https://github.com/PoonLab/Poplars
#   - HIVSeqinR:              https://github.com/guineverelee/HIVSeqinR
#   - HIVIntact ("proviral"): https://github.com/ramics/HIVIntact
#     (the actual pip-installable package inside this repo is named
#     "intactness-pipeline" per its own README; the CLI it installs is
#     called `proviral`)
# strict mode: -e aborts on any error since a broken clone/install shouldn't proceed
set -euo pipefail
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"           # absolute path of this script's own dir
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"     # repo top-level (three dirs up)
mkdir -p "${REPO_ROOT}/scripts/tools"                # ensure the shared tools dir exists
cd "${REPO_ROOT}/scripts/tools"                      # clone everything into scripts/tools/

# only clone/install Poplars if not already present (idempotent)
if [ ! -d Poplars ]; then
    git clone https://github.com/PoonLab/Poplars.git # fetch the Poplars source
    # editable install so its deps resolve (run via hypermut.py directly)
    pip install -e ./Poplars
fi

if [ ! -d HIVSeqinR ]; then                          # only clone HIVSeqinR if missing
    # fetch it; it's an R script run in-place, so no install step
    git clone https://github.com/guineverelee/HIVSeqinR.git
fi

# The clone alone is NOT enough to run HIVSeqinR, which is why it had never
# executed: the repo ships only release ZIPs, so R_HIVSeqinR_Combined_ver*.R does
# not exist until one is unpacked, the R `muscle` dependency is absent, and the
# script's MyBlastnDir points at the author's own machine. All three are handled
# here so the tool is provisioned reproducibly rather than by hand.
HIVSEQINR_DIR="$(pwd)/HIVSeqinR"
# 1. unpack the highest-numbered release ZIP if the R script is not there yet
if ! compgen -G "${HIVSEQINR_DIR}/R_HIVSeqinR_Combined_ver*.R" > /dev/null; then
    # sort -V so ver2.7.1 beats ver2.6.4 rather than sorting lexically
    LATEST_ZIP=$(ls -1 "${HIVSEQINR_DIR}"/HIVSeqinR_ver*.zip 2>/dev/null | sort -V | tail -1)
    if [ -n "${LATEST_ZIP}" ]; then
        echo "Unpacking $(basename "${LATEST_ZIP}")"
        # -x __MACOSX skips the AppleDouble metadata copies in these archives
        unzip -o -q "${LATEST_ZIP}" -x "__MACOSX/*" -d "${HIVSEQINR_DIR}"
        rm -rf "${HIVSEQINR_DIR}/__MACOSX"
    else
        echo "WARNING: no HIVSeqinR_ver*.zip found to unpack." >&2
    fi
fi

# 2. R dependencies. `muscle` is a hard requirement and absent from the conda env.
#    `pwalign` is needed because Bioconductor moved pairwiseAlignment() out of
#    Biostrings and made it formally defunct there in Biostrings >= 2.77.1; this
#    2019 script calls it unqualified in 9 places and dies with
#    "pairwiseAlignment() has moved from Biostrings to the pwalign package".
Rscript -e 'for (p in c("muscle", "pwalign")) {
              if (!requireNamespace(p, quietly=TRUE)) {
                  if (!requireNamespace("BiocManager", quietly=TRUE))
                      install.packages("BiocManager", repos="https://cloud.r-project.org")
                  BiocManager::install(p, ask=FALSE, update=FALSE)
              }
            }' || echo "WARNING: could not install the R muscle/pwalign packages." >&2

# 3. HIVSeqinR blasts against an HXB2 database it expects to find in MyBlastnDir,
#    built from R_HXB2.fasta renamed to HXB2.fasta (per the script's own header)
BLASTDB_DIR="${HIVSEQINR_DIR}/blastdb"
if [ -s "${HIVSEQINR_DIR}/R_HXB2.fasta" ] && [ ! -s "${BLASTDB_DIR}/HXB2.fasta.nin" ]; then
    mkdir -p "${BLASTDB_DIR}"
    cp "${HIVSEQINR_DIR}/R_HXB2.fasta" "${BLASTDB_DIR}/HXB2.fasta"
    makeblastdb -in "${BLASTDB_DIR}/HXB2.fasta" -parse_seqids -dbtype nucl > /dev/null \
        || echo "WARNING: makeblastdb failed for the HXB2 database." >&2
fi

# 4. point MyBlastnDir at that database instead of the author's /Users/guin/.
#    Patched with sed rather than by hand so a fresh checkout is reproducible;
#    the trailing slash matters because the R script concatenates paths directly.
RSCRIPT_FILE=$(compgen -G "${HIVSEQINR_DIR}/R_HIVSeqinR_Combined_ver*.R" | head -1 || true)
if [ -n "${RSCRIPT_FILE}" ]; then
    sed -i "s|^MyBlastnDir <- .*|MyBlastnDir <- \"${BLASTDB_DIR}/\" #set by setup_tools.sh|" "${RSCRIPT_FILE}"
    # 4b. Bioconductor moved the whole pairwise-alignment API out of Biostrings
    #     into pwalign and made the originals formally defunct in Biostrings
    #     >= 2.77.1. This 2019 script calls several of them unqualified
    #     (pairwiseAlignment, pattern, subject, ...) and each one aborts the run.
    #
    #     Rather than chase individual names, pwalign is attached immediately
    #     AFTER every library(Biostrings) call, so it always sits ahead of
    #     Biostrings in the search path and its live functions mask the defunct
    #     stubs. Ordering matters: attaching it before Biostrings would let a
    #     later library(Biostrings) shadow it again.
    #
    #     The grep guard keeps setup idempotent -- re-running must not insert
    #     the same line repeatedly.
    if ! grep -q "^library(pwalign)" "${RSCRIPT_FILE}"; then
        sed -i 's|^library(Biostrings).*|&\nlibrary(pwalign) #added by setup_tools.sh: pairwise-alignment API moved out of Biostrings|' \
            "${RSCRIPT_FILE}"
    fi
    N_ATTACH=$(grep -c "^library(pwalign)" "${RSCRIPT_FILE}" || true)
    echo "pwalign attached after ${N_ATTACH} library(Biostrings) call(s)"
    # Record what was configured automatically AND what was left at its default,
    # so the marker file is an audit trail rather than a bare "trust me" flag.
    # The primers are the author's own 2nd-round PCR primers; the ones used for
    # the PacBio SMRTcap libraries are not documented in the SRA metadata, so
    # they are left alone and autotrim is expected to find no flanking primer.
    {
        echo "Configured by setup_tools.sh on $(date +%Y-%m-%d)"
        echo "MyBlastnDir     = ${BLASTDB_DIR}/  (HXB2 blast db built here)"
        echo "R script        = $(basename "${RSCRIPT_FILE}")"
        echo "Primer2ndF/R    = LEFT AT AUTHOR DEFAULTS -- SMRTcap primers unknown."
        echo "                  Treat primer autotrim as NOT performed."
    } > "${HIVSEQINR_DIR}/.CONFIGURED"
fi

if [ ! -d HIVIntact ]; then                          # only set up HIVIntact if missing
    # clone with submodules (subtype alignments live in a submodule)
    git clone --recurse-submodules https://github.com/ramics/HIVIntact.git
    # dedicated venv so proviral's deps don't clash with the base env
    python3 -m venv HIVIntact/env
    # Python 3.12's venv module no longer bundles setuptools, but the
    # installed `proviral` CLI imports pkg_resources (from setuptools) --
    # confirmed by the ModuleNotFoundError this throws without it.
    # add setuptools back so pkg_resources is importable at runtime
    HIVIntact/env/bin/pip install setuptools
    # install the intactness-pipeline package, providing the `proviral` CLI
    HIVIntact/env/bin/pip install ./HIVIntact
fi

# print post-setup warnings the user must read
echo "Setup complete. IMPORTANT caveats before running compare_biological_filtering_illumina.sh:"
echo "1. HIVSeqinR is NOT a CLI tool -- it's an R script meant to be run in"
echo "   RStudio ('highlight all, run all'). It has hardcoded config at the"
echo "   top of R_HIVSeqinR_Combined_ver04.R (blast DB path, primer"
echo "   sequences) that must be edited before it will run, and it REQUIRES"
echo "   input FASTA with no dashes or IUPAC ambiguity codes -- but bcftools"
echo "   consensus output legitimately contains IUPAC codes (R/Y/W/etc) at"
echo "   heterozygous positions. run_hivseqinr.sh now SPLITS each record on"
echo "   runs of N and classifies only the sequenced segments; it no longer"
echo "   substitutes N->A, which fabricated bases and mimicked the APOBEC"
echo "   G->A hypermutation signature."
echo "   Provisioning (release ZIP, R muscle, HXB2 blast DB, MyBlastnDir) is"
echo "   now automatic above -- see HIVSeqinR/.CONFIGURED for what was set."
echo "   STILL AT DEFAULTS: the 2nd-round PCR primers (SMRTcap primers are"
echo "   undocumented), so report primer autotrim as not performed."
echo "2. HIVIntact requires a --subtype argument from a fixed list in"
echo "   util/subtype_alignments -- confirm A1/D are actually present before"
echo "   trusting its output on this cohort (the tools review already flags"
echo "   its heuristics as subtype-B-optimised)."
