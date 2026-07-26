#!/bin/bash
# Extract the 5' LTR U3 region from a MAFFT-aligned multi-FASTA, anchored to
# HXB2's own annotated genome coordinates.
#
# Why this exists: `seqkit subseq -r 8677:9121` (an earlier approach) applied
# a raw HXB2 *linear genome* coordinate directly to *alignment columns*. That
# only works if the alignment introduces zero gaps before that point, which
# MAFFT does not guarantee once divergent subtype A1/D/recombinant sequences
# are included. This script instead:
#
#   1. Requires HXB2 (K03455.1) to be included as one of the aligned records
#      (the production pipelines prepend it to the pre-MAFFT combined FASTA).
#   2. Fetches HXB2's own GenBank annotation (cached locally) and reads its
#      real, curated 5' LTR and R-region boundaries instead of hardcoding
#      numbers from memory:
#          repeat_region   1..634      /note="5' LTR"
#          repeat_region   454..551    /note="R repeat 5' copy"
#      U3 is defined as the LTR up to (not including) the R region, i.e.
#      [LTR_start, R_start) in genome coordinates -- this matches the
#      proposal's own definition of U3 as the region immediately upstream of
#      the transcription start site.
#   3. Walks HXB2's own row in the alignment to convert those genome
#      coordinates into alignment-column coordinates (a standard ungapped ->
#      gapped coordinate liftover), then slices every record in the
#      alignment at those columns.
#   4. Strips gaps per-record to produce each sample's own ungapped U3
#      sequence, ready for motif scanning (FIMO/TFBSTools/etc).
#   5. Runs a cheap outlier safety net: flags (does not silently drop) any
#      extracted sequence whose length falls outside a tolerance band, or
#      whose identity to HXB2's own U3 falls outside a tolerance band.
#      Identity is computed directly from the shared alignment columns
#      (both sequences are already aligned to each other via the input MSA,
#      so a fresh pairwise realignment would just re-derive what MAFFT
#      already determined).
#
# Usage:
#   extract_u3_by_hxb2_anchor.sh \
#       --alignment aligned_consensus.fasta \
#       --hxb2-id K03455.1 \
#       --gb-cache reference/K03455.1.gb \
#       --out-gapped motifs/U3_aligned.fasta \
#       --out motifs/U3_extracted.fasta \
#       --warnings-log motifs/u3_extraction_warnings.log
set -uo pipefail                                     # -u errors on unset vars, pipefail fails a pipe if any stage fails

EUTILS_URL_TMPL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%s&rettype=gb&retmode=text"  # NCBI efetch URL template to pull a GenBank record by accession

# Length/identity tolerance for the outlier safety net.
MIN_U3_LEN=400                                       # shortest plausible U3; anything below is flagged for review
MAX_U3_LEN=500                                       # longest plausible U3; anything above is flagged for review
MIN_IDENTITY_PCT=70.0                                # minimum % identity to HXB2 U3 before a sequence is flagged

HXB2_ID="K03455.1"                                   # default HXB2 accession; overridable via --hxb2-id

while [ $# -gt 0 ]; do                               # parse long-form CLI options
    case "$1" in                                     # dispatch on the current flag
        --alignment) ALIGNMENT="$2"; shift 2 ;;      # MAFFT-aligned multi-FASTA input
        --hxb2-id) HXB2_ID="$2"; shift 2 ;;          # override the HXB2 record id
        --gb-cache) GB_CACHE="$2"; shift 2 ;;        # local path to cache HXB2's GenBank record
        --out-gapped) OUT_GAPPED="$2"; shift 2 ;;    # output: U3 slice keeping alignment gaps
        --out) OUT="$2"; shift 2 ;;                  # output: per-sample ungapped U3 sequences
        --warnings-log) WARNINGS_LOG="$2"; shift 2 ;;  # output: outlier warnings log
        *) echo "ERROR: unknown argument $1" >&2; exit 1 ;;  # reject anything unrecognised
    esac
done

for var in ALIGNMENT GB_CACHE OUT_GAPPED OUT WARNINGS_LOG; do  # every one of these is mandatory
    if [ -z "${!var:-}" ]; then                      # indirect-expand each name; empty = not supplied
        echo "ERROR: missing required argument for ${var}" >&2  # say which one is missing
        exit 1                                       # and abort
    fi
done

if [ ! -s "${ALIGNMENT}" ]; then                     # the alignment must exist and be non-empty
    echo "ERROR: no alignment found at ${ALIGNMENT}" >&2  # otherwise there's nothing to slice
    exit 1
fi

if [ ! -s "${GB_CACHE}" ]; then                      # fetch HXB2's GenBank record only if not already cached
    mkdir -p "$(dirname "${GB_CACHE}")"              # ensure the cache dir exists
    URL=$(printf "${EUTILS_URL_TMPL}" "${HXB2_ID}")  # fill the accession into the efetch URL
    curl -s "${URL}" -o "${GB_CACHE}"                # download the GenBank record quietly to the cache
    if [ ! -s "${GB_CACHE}" ]; then                  # guard against an empty/failed download
        echo "ERROR: empty GenBank response for ${HXB2_ID}" >&2  # report it
        exit 1                                       # and abort rather than parse garbage
    fi
fi

# --- Find HXB2's own annotated 5' LTR and R-region start coordinates ---
# Feature-key lines (e.g. "     repeat_region   1..634") have the feature
# type as awk's $1 and location as $2 regardless of leading whitespace;
# qualifier lines (e.g. "                     /note=\"5' LTR\"") always
# start with optional whitespace then '/'. "3'" copies are excluded by
# requiring the note NOT contain "3" -- the only two repeat_region notes in
# HXB2's record are the 5'/3' LTR pair and the 5'/3' R-repeat pair.
read -r LTR_START R_START <<EOF                      # capture the two 0-based coordinates awk prints below
$(awk '
/^[[:space:]]*\// {                                  # a qualifier line (starts with optional whitespace then /)
    if ($0 ~ /\/note=/ && pending_type == "repeat_region") {  # a /note on the repeat_region feature we just saw
        note = $0                                    # copy the line to extract the note text
        sub(/.*\/note="/, "", note)                  # drop everything up to the opening quote
        sub(/".*/, "", note)                         # drop the closing quote and beyond, leaving the note text
        notelow = tolower(note)                      # lowercase for case-insensitive matching
        if (index(notelow, "ltr") > 0 && index(notelow, "3") == 0 && ltr_loc == "") {  # 5-prime LTR note (excludes the 3-prime copy), first only
            ltr_loc = pending_loc                    # remember the LTR location string (e.g. 1..634)
        } else if (index(notelow, "r repeat") > 0 && index(notelow, "3") == 0 && r_loc == "") {  # 5-prime R-repeat note (excludes 3-prime copy), first only
            r_loc = pending_loc                      # remember the R-repeat location string
        }
    }
    next                                             # qualifier line consumed; get the next line
}
NF >= 2 && $1 ~ /^[A-Za-z_]+$/ {                     # a feature-key line: type in $1, location in $2
    pending_type = $1                                # stash the feature type for the following qualifier lines
    pending_loc = $2                                 # stash its location string
    next
}
END {
    if (ltr_loc == "" || r_loc == "") {              # both features must have been found
        print "ERROR: could not find both a 5-prime LTR and R repeat 5-prime repeat_region feature" > "/dev/stderr"  # else report failure
        exit 1
    }
    split(ltr_loc, a, "\\.\\.")                      # split "start..end" of the LTR
    split(r_loc, b, "\\.\\.")                        # split "start..end" of the R-repeat
    print a[1]-1, b[1]-1                              # print both starts converted 1-based -> 0-based
}
' "${GB_CACHE}")
EOF
if [ -z "${LTR_START:-}" ] || [ -z "${R_START:-}" ]; then  # awk must have produced both coordinates
    echo "ERROR: failed to parse 5' LTR / R-repeat coordinates from ${GB_CACHE}. Inspect its FEATURES table and adjust the note-matching logic if this reference's annotation conventions differ from HXB2's." >&2  # guide the user if parsing failed
    exit 1
fi
if [ "${R_START}" -le "${LTR_START}" ]; then         # U3 = [LTR_start, R_start), so R must come after the LTR start
    echo "ERROR: R-region start (${R_START}) is not after LTR start (${LTR_START}) -- unexpected annotation, refusing to guess." >&2  # sanity-fail rather than emit nonsense coords
    exit 1
fi
echo "HXB2 U3 genome coordinates (0-based, half-open): [${LTR_START}, ${R_START}) ($((R_START - LTR_START)) nt)" >&2  # log the resolved U3 window and its length

mkdir -p "$(dirname "${OUT_GAPPED}")" "$(dirname "${OUT}")" "$(dirname "${WARNINGS_LOG}")"  # ensure all three output dirs exist

# --- Liftover HXB2's genome coordinates to alignment columns, slice every
# record, strip gaps, and run the outlier safety net -- all in one pass
# over the buffered alignment. ---
awk -v hxb2_id="${HXB2_ID}" -v u3_start="${LTR_START}" -v u3_end="${R_START}" \
    -v min_len="${MIN_U3_LEN}" -v max_len="${MAX_U3_LEN}" -v min_identity="${MIN_IDENTITY_PCT}" \
    -v out_gapped="${OUT_GAPPED}" -v out="${OUT}" -v warnings_log="${WARNINGS_LOG}" '  # pass all resolved coords/thresholds/paths into awk
/^>/ {                                               # a FASTA header line
    n++                                              # count of records seen so far
    header = substr($0, 2)                           # header text without the leading >
    split(header, htok, /[ \t]/)                     # split on whitespace to isolate the id token
    ids[n] = htok[1]                                 # store this record id
    seqs[n] = ""                                     # start with an empty sequence buffer for it
    next
}
{ seqs[n] = seqs[n] $0 }                             # non-header line: append to the current record sequence
END {
    hxb2_idx = 0                                     # index of the HXB2 record in the alignment (0 = not yet found)
    for (i = 1; i <= n; i++) {                        # scan all records for HXB2
        if (ids[i] == hxb2_id || index(ids[i], hxb2_id ".") == 1) { hxb2_idx = i; break }  # exact id or versioned prefix match
    }
    if (hxb2_idx == 0) {                             # HXB2 is the coordinate anchor, so it must be present
        print "ERROR: HXB2 record " hxb2_id " not found in the alignment. The pipeline must include HXB2 in the pre-MAFFT combined FASTA so the alignment shares a common column-coordinate space." > "/dev/stderr"  # explain the requirement
        exit 1
    }

    hxb2_seq = seqs[hxb2_idx]                         # HXB2 aligned (gapped) row
    L = length(hxb2_seq)                              # alignment width in columns
    k = 0                                             # count of ungapped HXB2 bases walked so far (genome coordinate)
    col_start = 0                                     # alignment column mapping to U3 start (found below)
    col_end = 0                                       # alignment column mapping to U3 end (found below)
    for (c = 1; c <= L; c++) {                        # walk HXB2 column by column to lift genome coords to columns
        if (substr(hxb2_seq, c, 1) != "-") {          # only real bases advance the genome coordinate
            if (k == u3_start) col_start = c          # column where the U3 start base sits
            if (k == u3_end - 1) col_end = c          # column of the last base before R (inclusive U3 end)
            k++                                       # one more ungapped base consumed
        }
    }
    if (col_start == 0 || col_end == 0) {            # both boundaries must have been located
        print "ERROR: HXB2 aligned row has only " k " ungapped bases, fewer than the U3 end coordinate (" u3_end "). The alignment'\''s HXB2 record looks truncated or wrong." > "/dev/stderr"  # HXB2 row too short/wrong
        exit 1
    }
    print "Alignment-column range for U3 (0-based, inclusive): [" col_start-1 ", " col_end-1 "]" > "/dev/stderr"  # log the resolved column window

    printf "" > out_gapped                            # truncate/create the gapped output file
    printf "" > out                                   # truncate/create the ungapped output file
    n_warnings = 0                                    # count of flagged sequences

    for (i = 1; i <= n; i++) {                        # slice the U3 window out of every record
        gslice[i] = substr(seqs[i], col_start, col_end - col_start + 1)  # gapped U3 slice (same columns for all records)
        uslice = gslice[i]                            # copy to make an ungapped version
        gsub(/-/, "", uslice)                         # remove alignment gaps -> this sample own U3 sequence
        uslice_arr[i] = uslice                        # keep the ungapped slice for the safety-net checks below

        print ">" ids[i] >> out_gapped                # write gapped record header
        print toupper(gslice[i]) >> out_gapped        # write gapped U3 (uppercased)
        print ">" ids[i] >> out                        # write ungapped record header
        print toupper(uslice) >> out                  # write ungapped U3 (uppercased)
    }
    close(out_gapped)                                 # flush/close so file sizes are final
    close(out)

    hxb2_gapped = gslice[hxb2_idx]                    # HXB2 own gapped U3, the identity reference
    hxb2_ungapped_len = length(uslice_arr[hxb2_idx])  # denominator for the identity percentage

    for (i = 1; i <= n; i++) {                        # outlier safety net over every extracted U3
        flags = ""                                    # accumulate reasons this sequence looks off
        seqlen = length(uslice_arr[i])                # its ungapped U3 length
        if (seqlen < min_len || seqlen > max_len) {   # length outside the tolerance band?
            flags = flags "; length " seqlen "nt outside [" min_len "," max_len "]"  # note it
        }
        if (i != hxb2_idx && hxb2_ungapped_len > 0) { # compute identity to HXB2 (skip HXB2 vs itself)
            matches = 0                               # count of matching aligned positions
            for (p = 1; p <= length(hxb2_gapped); p++) {  # compare column by column over the shared U3 slice
                hb = substr(hxb2_gapped, p, 1)        # HXB2 base at this column
                if (hb == "-") continue               # skip columns where HXB2 has a gap
                ob = substr(gslice[i], p, 1)          # this sample base at the same column
                if (toupper(ob) == toupper(hb)) matches++  # count a match (case-insensitive)
            }
            identity_pct = 100.0 * matches / hxb2_ungapped_len  # % identity relative to HXB2 U3 length
            if (identity_pct < min_identity + 0) {    # below threshold? (+0 forces numeric compare)
                flags = flags sprintf("; only %.1f%% identity to HXB2 U3 (threshold %s%%)", identity_pct, min_identity)  # note it
            }
        }
        if (flags != "") {                            # if anything was flagged...
            sub(/^; /, "", flags)                     # tidy the leading separator
            n_warnings++                              # bump the warning count
            warnings[n_warnings] = ids[i] ": " flags  # store the per-sequence warning
        }
    }

    printf "" > warnings_log                          # truncate/create the warnings log
    for (w = 1; w <= n_warnings; w++) print warnings[w] >> warnings_log  # write each flagged sequence
    close(warnings_log)

    if (n_warnings > 0) {                             # report the safety-net outcome
        print "WARNING: " n_warnings " sequence(s) flagged for manual review, see " warnings_log > "/dev/stderr"  # some flagged (not dropped)
    } else {
        print "All extracted U3 sequences passed the outlier safety net." > "/dev/stderr"  # all clean
    }
    print "Wrote " n " U3 sequences to " out > "/dev/stderr"  # final count written
}
' "${ALIGNMENT}"                                     # run the whole liftover/slice/check pass over the alignment
