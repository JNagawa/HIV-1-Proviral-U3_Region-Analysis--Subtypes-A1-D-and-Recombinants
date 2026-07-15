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
set -uo pipefail

EUTILS_URL_TMPL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%s&rettype=gb&retmode=text"

# Length/identity tolerance for the outlier safety net.
MIN_U3_LEN=400
MAX_U3_LEN=500
MIN_IDENTITY_PCT=70.0

HXB2_ID="K03455.1"

while [ $# -gt 0 ]; do
    case "$1" in
        --alignment) ALIGNMENT="$2"; shift 2 ;;
        --hxb2-id) HXB2_ID="$2"; shift 2 ;;
        --gb-cache) GB_CACHE="$2"; shift 2 ;;
        --out-gapped) OUT_GAPPED="$2"; shift 2 ;;
        --out) OUT="$2"; shift 2 ;;
        --warnings-log) WARNINGS_LOG="$2"; shift 2 ;;
        *) echo "ERROR: unknown argument $1" >&2; exit 1 ;;
    esac
done

for var in ALIGNMENT GB_CACHE OUT_GAPPED OUT WARNINGS_LOG; do
    if [ -z "${!var:-}" ]; then
        echo "ERROR: missing required argument for ${var}" >&2
        exit 1
    fi
done

if [ ! -s "${ALIGNMENT}" ]; then
    echo "ERROR: no alignment found at ${ALIGNMENT}" >&2
    exit 1
fi

if [ ! -s "${GB_CACHE}" ]; then
    mkdir -p "$(dirname "${GB_CACHE}")"
    URL=$(printf "${EUTILS_URL_TMPL}" "${HXB2_ID}")
    curl -s "${URL}" -o "${GB_CACHE}"
    if [ ! -s "${GB_CACHE}" ]; then
        echo "ERROR: empty GenBank response for ${HXB2_ID}" >&2
        exit 1
    fi
fi

# --- Find HXB2's own annotated 5' LTR and R-region start coordinates ---
# Feature-key lines (e.g. "     repeat_region   1..634") have the feature
# type as awk's $1 and location as $2 regardless of leading whitespace;
# qualifier lines (e.g. "                     /note=\"5' LTR\"") always
# start with optional whitespace then '/'. "3'" copies are excluded by
# requiring the note NOT contain "3" -- the only two repeat_region notes in
# HXB2's record are the 5'/3' LTR pair and the 5'/3' R-repeat pair.
read -r LTR_START R_START <<EOF
$(awk '
/^[[:space:]]*\// {
    if ($0 ~ /\/note=/ && pending_type == "repeat_region") {
        note = $0
        sub(/.*\/note="/, "", note)
        sub(/".*/, "", note)
        notelow = tolower(note)
        if (index(notelow, "ltr") > 0 && index(notelow, "3") == 0 && ltr_loc == "") {
            ltr_loc = pending_loc
        } else if (index(notelow, "r repeat") > 0 && index(notelow, "3") == 0 && r_loc == "") {
            r_loc = pending_loc
        }
    }
    next
}
NF >= 2 && $1 ~ /^[A-Za-z_]+$/ {
    pending_type = $1
    pending_loc = $2
    next
}
END {
    if (ltr_loc == "" || r_loc == "") {
        print "ERROR: could not find both a 5-prime LTR and R repeat 5-prime repeat_region feature" > "/dev/stderr"
        exit 1
    }
    split(ltr_loc, a, "\\.\\.")
    split(r_loc, b, "\\.\\.")
    print a[1]-1, b[1]-1
}
' "${GB_CACHE}")
EOF
if [ -z "${LTR_START:-}" ] || [ -z "${R_START:-}" ]; then
    echo "ERROR: failed to parse 5' LTR / R-repeat coordinates from ${GB_CACHE}. Inspect its FEATURES table and adjust the note-matching logic if this reference's annotation conventions differ from HXB2's." >&2
    exit 1
fi
if [ "${R_START}" -le "${LTR_START}" ]; then
    echo "ERROR: R-region start (${R_START}) is not after LTR start (${LTR_START}) -- unexpected annotation, refusing to guess." >&2
    exit 1
fi
echo "HXB2 U3 genome coordinates (0-based, half-open): [${LTR_START}, ${R_START}) ($((R_START - LTR_START)) nt)" >&2

mkdir -p "$(dirname "${OUT_GAPPED}")" "$(dirname "${OUT}")" "$(dirname "${WARNINGS_LOG}")"

# --- Liftover HXB2's genome coordinates to alignment columns, slice every
# record, strip gaps, and run the outlier safety net -- all in one pass
# over the buffered alignment. ---
awk -v hxb2_id="${HXB2_ID}" -v u3_start="${LTR_START}" -v u3_end="${R_START}" \
    -v min_len="${MIN_U3_LEN}" -v max_len="${MAX_U3_LEN}" -v min_identity="${MIN_IDENTITY_PCT}" \
    -v out_gapped="${OUT_GAPPED}" -v out="${OUT}" -v warnings_log="${WARNINGS_LOG}" '
/^>/ {
    n++
    header = substr($0, 2)
    split(header, htok, /[ \t]/)
    ids[n] = htok[1]
    seqs[n] = ""
    next
}
{ seqs[n] = seqs[n] $0 }
END {
    hxb2_idx = 0
    for (i = 1; i <= n; i++) {
        if (ids[i] == hxb2_id || index(ids[i], hxb2_id ".") == 1) { hxb2_idx = i; break }
    }
    if (hxb2_idx == 0) {
        print "ERROR: HXB2 record " hxb2_id " not found in the alignment. The pipeline must include HXB2 in the pre-MAFFT combined FASTA so the alignment shares a common column-coordinate space." > "/dev/stderr"
        exit 1
    }

    hxb2_seq = seqs[hxb2_idx]
    L = length(hxb2_seq)
    k = 0
    col_start = 0
    col_end = 0
    for (c = 1; c <= L; c++) {
        if (substr(hxb2_seq, c, 1) != "-") {
            if (k == u3_start) col_start = c
            if (k == u3_end - 1) col_end = c
            k++
        }
    }
    if (col_start == 0 || col_end == 0) {
        print "ERROR: HXB2 aligned row has only " k " ungapped bases, fewer than the U3 end coordinate (" u3_end "). The alignment'\''s HXB2 record looks truncated or wrong." > "/dev/stderr"
        exit 1
    }
    print "Alignment-column range for U3 (0-based, inclusive): [" col_start-1 ", " col_end-1 "]" > "/dev/stderr"

    printf "" > out_gapped
    printf "" > out
    n_warnings = 0

    for (i = 1; i <= n; i++) {
        gslice[i] = substr(seqs[i], col_start, col_end - col_start + 1)
        uslice = gslice[i]
        gsub(/-/, "", uslice)
        uslice_arr[i] = uslice

        print ">" ids[i] >> out_gapped
        print gslice[i] >> out_gapped
        print ">" ids[i] >> out
        print uslice >> out
    }
    close(out_gapped)
    close(out)

    hxb2_gapped = gslice[hxb2_idx]
    hxb2_ungapped_len = length(uslice_arr[hxb2_idx])

    for (i = 1; i <= n; i++) {
        flags = ""
        seqlen = length(uslice_arr[i])
        if (seqlen < min_len || seqlen > max_len) {
            flags = flags "; length " seqlen "nt outside [" min_len "," max_len "]"
        }
        if (i != hxb2_idx && hxb2_ungapped_len > 0) {
            matches = 0
            for (p = 1; p <= length(hxb2_gapped); p++) {
                hb = substr(hxb2_gapped, p, 1)
                if (hb == "-") continue
                ob = substr(gslice[i], p, 1)
                if (toupper(ob) == toupper(hb)) matches++
            }
            identity_pct = 100.0 * matches / hxb2_ungapped_len
            if (identity_pct < min_identity + 0) {
                flags = flags sprintf("; only %.1f%% identity to HXB2 U3 (threshold %s%%)", identity_pct, min_identity)
            }
        }
        if (flags != "") {
            sub(/^; /, "", flags)
            n_warnings++
            warnings[n_warnings] = ids[i] ": " flags
        }
    }

    printf "" > warnings_log
    for (w = 1; w <= n_warnings; w++) print warnings[w] >> warnings_log
    close(warnings_log)

    if (n_warnings > 0) {
        print "WARNING: " n_warnings " sequence(s) flagged for manual review, see " warnings_log > "/dev/stderr"
    } else {
        print "All extracted U3 sequences passed the outlier safety net." > "/dev/stderr"
    }
    print "Wrote " n " U3 sequences to " out > "/dev/stderr"
}
' "${ALIGNMENT}"
