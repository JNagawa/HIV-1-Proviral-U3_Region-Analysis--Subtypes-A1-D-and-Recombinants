#!/bin/bash
# HIVSeqinR is NOT a normal CLI tool -- confirmed from its own README: it's
# an R script meant to be run inside RStudio ("highlight all, run all"),
# with configuration (blast DB directory, 2nd-round PCR primer sequences,
# amino-acid length cutoffs) hardcoded at the top of the script rather than
# passed as arguments. This wrapper does NOT attempt to auto-patch that
# config, since guessing at variable names risks silently misconfiguring a
# tool that classifies genome intactness. Instead it requires you to have
# manually edited tools/HIVSeqinR/R_HIVSeqinR_Combined_ver*.R yourself and
# confirmed that by creating a CONFIGURED marker file.
#
# HIVSeqinR also requires input FASTA with no dashes or IUPAC ambiguity
# codes -- bcftools consensus output legitimately contains IUPAC codes
# (R/Y/W/etc) at heterozygous positions, so this wrapper resolves each
# ambiguity code to one of its bases (arbitrarily, the first alphabetically)
# before copying sequences in. Record this substitution as a real limitation
# in the methods write-up, not a transparent equivalent to the true sequence.
#
# Usage: run_hivseqinr.sh <INPUT_FASTA> <OUTDIR>
set -uo pipefail
IN="$1" OUTDIR="$2"

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
HIVSEQINR_DIR="${STAGE_DIR}/tools/HIVSeqinR"
CONFIGURED_MARKER="${HIVSEQINR_DIR}/.CONFIGURED"

if [ ! -f "${CONFIGURED_MARKER}" ]; then
    echo "ERROR: ${HIVSEQINR_DIR} has not been configured yet." >&2
    echo "Open R_HIVSeqinR_Combined_ver*.R in RStudio, set MyBlastnDir and" >&2
    echo "the 2nd-round PCR primer sequences for your protocol, then run:" >&2
    echo "  touch ${CONFIGURED_MARKER}" >&2
    echo "to confirm you've done this before rerunning." >&2
    exit 1
fi

mkdir -p "${OUTDIR}"
RAW_FASTA_DIR="${HIVSEQINR_DIR}/RAW_FASTA"
mkdir -p "${RAW_FASTA_DIR}"
rm -f "${RAW_FASTA_DIR}"/*.seq

# Split into one file per record (id sanitized, sequence unwrapped), then
# resolve IUPAC ambiguity codes to a single concrete base (arbitrary,
# alphabetically-first choice per code -- R/Y/S/W/K/M/B/D/H/V/N below).
# Non-ambiguous bases (A/C/G/T, either case) pass through untouched.
awk -v outdir="${RAW_FASTA_DIR}" '
function flush(   safe_id, outfile) {
    if (id == "") return
    safe_id = id
    gsub(/-/, "_", safe_id)
    gsub(/\*/, "_", safe_id)
    outfile = outdir "/" safe_id ".seq"
    print ">" safe_id > outfile
    print seq > outfile
    close(outfile)
}
/^>/ {
    flush()
    header = substr($0, 2)
    split(header, tok, /[ \t]/)
    id = tok[1]
    seq = ""
    next
}
{ seq = seq $0 }
END { flush() }
' "${IN}"

for f in "${RAW_FASTA_DIR}"/*.seq; do
    [ -e "${f}" ] || continue
    HEADER=$(head -1 "${f}")
    SEQ=$(tail -n +2 "${f}" | tr 'RYSWKMBDHVNryswkmbdhvn' 'ACCAGACAAAAACCAGACAAAA')
    printf '%s\n%s\n' "${HEADER}" "${SEQ}" > "${f}"
done

cd "${HIVSEQINR_DIR}" || exit 1
RSCRIPT_FILE=$(compgen -G "R_HIVSeqinR_Combined_ver*.R" | head -1)
[ -n "${RSCRIPT_FILE}" ] || { echo "ERROR: could not find R_HIVSeqinR_Combined_ver*.R" >&2; exit 1; }

Rscript "${RSCRIPT_FILE}" > "${OUTDIR}/hivseqinr.log" 2>&1

RESULT_CSV="${HIVSEQINR_DIR}/Results_Final/Output_MyBigSummary_DF_FINAL.csv"
if [ -s "${RESULT_CSV}" ]; then
    cp "${RESULT_CSV}" "${OUTDIR}/"
else
    echo "ERROR: HIVSeqinR did not produce ${RESULT_CSV}, see ${OUTDIR}/hivseqinr.log" >&2
    exit 1
fi
