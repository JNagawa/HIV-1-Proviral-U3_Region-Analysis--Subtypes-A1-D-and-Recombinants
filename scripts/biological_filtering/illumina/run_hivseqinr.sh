#!/bin/bash
# HIVSeqinR is NOT a normal CLI tool -- confirmed from its own README: it's
# an R script meant to be run inside RStudio ("highlight all, run all"),
# with configuration (blast DB directory, 2nd-round PCR primer sequences,
# amino-acid length cutoffs) hardcoded at the top of the script rather than
# passed as arguments. This wrapper does NOT attempt to auto-patch that
# config, since guessing at variable names risks silently misconfiguring a
# tool that classifies genome intactness. Instead it requires you to have
# manually edited scripts/tools/HIVSeqinR/R_HIVSeqinR_Combined_ver*.R yourself and
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
# -u errors on unset vars, pipefail fails a pipe if any stage fails
set -uo pipefail
# positional args: input FASTA and where to copy results
IN="$1" OUTDIR="$2"

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"           # absolute path of this script's own dir
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"     # repo top-level (three dirs up)
# the cloned HIVSeqinR repo (where the R script runs from)
HIVSEQINR_DIR="${REPO_ROOT}/scripts/tools/HIVSeqinR"
# marker file the user creates to attest they've edited the hardcoded config
CONFIGURED_MARKER="${HIVSEQINR_DIR}/.CONFIGURED"

# refuse to run until the user has confirmed manual configuration
if [ ! -f "${CONFIGURED_MARKER}" ]; then
    echo "ERROR: ${HIVSEQINR_DIR} has not been configured yet." >&2  # explain the problem
    echo "Open R_HIVSeqinR_Combined_ver*.R in RStudio, set MyBlastnDir and" >&2  # ...
    # ...tell them exactly what to edit...
    echo "the 2nd-round PCR primer sequences for your protocol, then run:" >&2
    echo "  touch ${CONFIGURED_MARKER}" >&2           # ...and how to create the marker
    echo "to confirm you've done this before rerunning." >&2  # final instruction
    exit 1                                           # bail; tool not configured
fi

mkdir -p "${OUTDIR}"                                 # ensure the output dir exists
# HIVSeqinR reads one .seq file per sequence from this fixed input dir
RAW_FASTA_DIR="${HIVSEQINR_DIR}/RAW_FASTA"
mkdir -p "${RAW_FASTA_DIR}"                          # create it if needed
# clear any .seq files from a previous run so results aren't mixed
rm -f "${RAW_FASTA_DIR}"/*.seq

# Split into one file per record (id sanitized, sequence unwrapped), then
# resolve IUPAC ambiguity codes to a single concrete base (arbitrary,
# alphabetically-first choice per code -- R/Y/S/W/K/M/B/D/H/V/N below).
# Non-ambiguous bases (A/C/G/T, either case) pass through untouched.
awk -v outdir="${RAW_FASTA_DIR}" '                   # split the multi-FASTA into per-record .seq files (outdir passed in)
function flush(   safe_id, outfile) {                # write the currently-buffered record to its own file
    if (id == "") return                             # nothing buffered yet, skip
    safe_id = id                                     # copy id so we can sanitize it for use as a filename
    gsub(/-/, "_", safe_id)                          # replace dashes with underscores (safe filenames)
    gsub(/\*/, "_", safe_id)                         # replace asterisks with underscores too
    outfile = outdir "/" safe_id ".seq"              # per-record output path
    print ">" safe_id > outfile                      # write the sanitized FASTA header
    print seq > outfile                              # write the unwrapped sequence on one line
    close(outfile)                                   # close so we do not hit the open-file-descriptor limit
}
/^>/ {                                               # on each FASTA header line...
    flush()                                          # ...write out the previous record first
    header = substr($0, 2)                           # strip the leading ">"
    split(header, tok, /[ \t]/)                      # split on whitespace to isolate the id from any description
    id = tok[1]                                      # use the first token as the record id
    seq = ""                                         # reset the sequence buffer for this record
    next                                             # done with the header line
}
{ seq = seq $0 }                                     # accumulate (unwrap) sequence lines into the buffer
END { flush() }                                      # flush the final buffered record at end of input
' "${IN}"

# process each per-record file to strip IUPAC ambiguity codes
for f in "${RAW_FASTA_DIR}"/*.seq; do
    # skip if the glob matched nothing (no .seq files)
    [ -e "${f}" ] || continue
    HEADER=$(head -1 "${f}")                          # keep the header line as-is
    # map each ambiguity code to its alphabetically-first base
    SEQ=$(tail -n +2 "${f}" | tr 'RYSWKMBDHVNryswkmbdhvn' 'ACCAGACAAAAACCAGACAAAA')
    printf '%s\n%s\n' "${HEADER}" "${SEQ}" > "${f}"  # rewrite the file with the resolved sequence
done

# the R script uses relative paths, so run from the repo dir
cd "${HIVSEQINR_DIR}" || exit 1
# find the versioned main R script (version number varies)
RSCRIPT_FILE=$(compgen -G "R_HIVSeqinR_Combined_ver*.R" | head -1)
# bail if the script is missing
[ -n "${RSCRIPT_FILE}" ] || { echo "ERROR: could not find R_HIVSeqinR_Combined_ver*.R" >&2; exit 1; }

# run the R pipeline headlessly, capturing all output to the log
Rscript "${RSCRIPT_FILE}" > "${OUTDIR}/hivseqinr.log" 2>&1

# where HIVSeqinR writes its final summary
RESULT_CSV="${HIVSEQINR_DIR}/Results_Final/Output_MyBigSummary_DF_FINAL.csv"
if [ -s "${RESULT_CSV}" ]; then                       # if the run produced a non-empty summary...
    cp "${RESULT_CSV}" "${OUTDIR}/"                   # ...copy it into our per-run output dir
else
    # otherwise report failure and point at the log
    echo "ERROR: HIVSeqinR did not produce ${RESULT_CSV}, see ${OUTDIR}/hivseqinr.log" >&2
    exit 1                                           # signal failure to the caller
fi
