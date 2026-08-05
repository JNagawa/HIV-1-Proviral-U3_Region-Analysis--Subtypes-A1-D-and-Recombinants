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
# HIVSeqinR also requires input FASTA with no dashes or IUPAC ambiguity codes
# ("No dashes, *, or IUPAC mixture symbols are allowed (eg. N, R, S, Y etc)",
# from the script's own header).
#
# This wrapper used to satisfy that by mapping every ambiguity code to a
# concrete base with tr, which sent N to A. That is unsafe here: since the
# PacBio assembly step began N-masking unsequenced positions, the consensuses
# carry 8-88% N, so the substitution would have fabricated up to 8552 A bases
# for a single sample. Worse, N->A synthesises precisely the APOBEC3G G->A
# hypermutation signature that intactness callers look for, so the result would
# not be noisy -- it would be a systematic false "hypermutated" call.
#
# Instead each record is SPLIT on runs of N and only the sequenced segments are
# passed through, so nothing is invented. A segment shorter than MIN_SEGMENT_LEN
# is dropped as too short to classify. Genuine heterozygous IUPAC codes (R/Y/W
# etc, which are real observations rather than absent data) are still resolved to
# a concrete base, and the count of such substitutions is reported.
#
# Consequence to carry into the write-up: a sample whose genome arrives in
# several short segments cannot receive a whole-genome intactness verdict, and
# HIVSeqinR will classify each segment on its own. That is a limit of the data,
# not of this wrapper.
#
# NOTE ON PRIMERS: Primer2ndF_HXB2/Primer2ndR_HXB2 in the R script are the
# author's own 2nd-round PCR primers. The primers for this SMRTcap dataset are
# not documented in its SRA metadata, so the defaults are left in place and
# autotrim is expected to find no flanks. Treat trimming as not performed.
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

# shortest sequenced segment worth classifying; below this there is no ORF
# structure to assess and HIVSeqinR would only report it as truncated
MIN_SEGMENT_LEN="${MIN_SEGMENT_LEN:-500}"
# per-record accounting of what was split out and what was substituted
SEGMENT_REPORT="${OUTDIR}/segments.tsv"

# Split every record on runs of N, writing one .seq file per sequenced segment.
# Nothing is invented: unsequenced stretches are dropped, not filled. Real
# heterozygous IUPAC codes are resolved to a concrete base and counted.
awk -v outdir="${RAW_FASTA_DIR}" -v minlen="${MIN_SEGMENT_LEN}" -v report="${SEGMENT_REPORT}" '
function emit(   safe_id, nseg, parts, i, seg, resolved, nsub, outid, outfile, kept, dropn, dropbp) {
    if (id == "") return
    safe_id = id
    gsub(/[^A-Za-z0-9_]/, "_", safe_id)              # HIVSeqinR only tolerates "_" as a special character
    # split on runs of one or more N; every remaining piece was actually sequenced
    nseg = split(toupper(seq), parts, /N+/)
    kept = 0; dropn = 0; dropbp = 0
    for (i = 1; i <= nseg; i++) {
        seg = parts[i]
        if (length(seg) < minlen) { dropn++; dropbp += length(seg); continue }
        kept++
        # Resolve genuine ambiguity codes. Unlike N these are real observations
        # (a mixed base that was sequenced), so collapsing them to one allele
        # loses information but invents nothing. Counted so the loss is visible.
        resolved = seg
        nsub = gsub(/[RYSWKMBDHV]/, "", resolved)     # count them
        resolved = seg                               # then do the real substitution
        gsub(/R/, "A", resolved); gsub(/Y/, "C", resolved); gsub(/S/, "C", resolved)
        gsub(/W/, "A", resolved); gsub(/K/, "G", resolved); gsub(/M/, "A", resolved)
        gsub(/B/, "C", resolved); gsub(/D/, "A", resolved); gsub(/H/, "A", resolved)
        gsub(/V/, "A", resolved)
        # a record yielding one segment keeps its plain id; several get _seg<N>
        outid = (nseg == 1) ? safe_id : safe_id "_seg" kept
        outfile = outdir "/" outid ".seq"
        print ">" outid > outfile
        print resolved > outfile
        close(outfile)
        printf "%s\t%s\t%d\t%d\n", id, outid, length(seg), nsub >> report
    }
    printf "  %-12s %d segment(s) kept, %d dropped as shorter than %dbp (%dbp total)\n", \
           id, kept, dropn, minlen, dropbp > "/dev/stderr"
}
BEGIN { printf "record\tsegment_id\tsegment_len\tiupac_substitutions\n" > report }
/^>/ { emit(); header = substr($0, 2); split(header, tok, /[ \t]/); id = tok[1]; seq = ""; next }
{ seq = seq $0 }
END { emit(); close(report) }
' "${IN}"

# HIVSeqinR requires at least one reference/positive control per run (its own
# README, April 2021 update), so append the bundled 8E5/HXB2 control if the
# caller did not already include a control sequence.
CONTROL_FASTA="${HIVSEQINR_DIR}/Examples_8E5_HXB2.fasta"
if [ -s "${CONTROL_FASTA}" ]; then
    awk -v outdir="${RAW_FASTA_DIR}" '
    function emit(   safe_id, outfile) {
        if (id == "") return
        safe_id = "CTRL_" id
        gsub(/[^A-Za-z0-9_]/, "_", safe_id)
        outfile = outdir "/" safe_id ".seq"
        print ">" safe_id > outfile
        print toupper(seq) > outfile
        close(outfile)
    }
    /^>/ { emit(); header = substr($0, 2); split(header, tok, /[ \t]/); id = tok[1]; seq = ""; next }
    { seq = seq $0 }
    END { emit() }
    ' "${CONTROL_FASTA}"
    echo "  appended positive control from $(basename "${CONTROL_FASTA}")" >&2
fi

# nothing classifiable means nothing to run -- fail loudly rather than let
# HIVSeqinR produce an empty summary that would read as a clean result
N_SEG=$(ls -1 "${RAW_FASTA_DIR}"/*.seq 2>/dev/null | wc -l)
if [ "${N_SEG}" -eq 0 ]; then
    echo "ERROR: no sequenced segment of at least ${MIN_SEGMENT_LEN}bp survived N-splitting; nothing to classify." >&2
    exit 1
fi
echo "  ${N_SEG} sequence file(s) staged for HIVSeqinR (see ${SEGMENT_REPORT})" >&2

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
