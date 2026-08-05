#!/bin/bash
# Extract the U3 region from a MAFFT-aligned multi-FASTA, anchored to HXB2's
# own annotated genome coordinates. Both LTR copies are considered and each
# record uses whichever one it actually sequenced (see step 4 below).
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
#      real, curated LTR and R-region boundaries for BOTH copies instead of
#      hardcoding numbers from memory:
#          repeat_region   1..634        /note="5' LTR"
#          repeat_region   454..551      /note="R repeat 5' copy"
#          repeat_region   9086..9719    /note="3' LTR"
#          repeat_region   9540..9636    /note="R repeat 3' copy"
#      U3 is defined as the LTR up to (not including) the R region, i.e.
#      [LTR_start, R_start) in genome coordinates -- this matches the
#      proposal's own definition of U3 as the region immediately upstream of
#      the transcription start site.
#   3. Walks HXB2's own row in the alignment to convert those genome
#      coordinates into alignment-column coordinates (a standard ungapped ->
#      gapped coordinate liftover), then slices every record in the
#      alignment at both column windows.
#   4. Picks, per record, whichever LTR copy carries more sequenced (ACGT)
#      bases, ties going to the 5' copy. The two LTRs are identical in an
#      integrated provirus, so either is a valid U3 source -- but an amplicon
#      library may only cover one of them, and taking the 5' copy blindly
#      returns reference-padded sequence for samples that only sequenced the
#      3' copy. The choice per record is written to <out>_ltr_choice.tsv.
#   5. Strips gaps per-record to produce each sample's own ungapped U3
#      sequence, ready for motif scanning (FIMO/TFBSTools/etc).
#   6. Runs a cheap outlier safety net: flags (does not silently drop) any
#      extracted sequence whose length falls outside a tolerance band, whose
#      identity to HXB2's own U3 falls outside a tolerance band, or which has
#      no sequenced base in either LTR. Identity is scored against HXB2's copy
#      from the SAME window the record used, and is computed directly from the
#      shared alignment columns (both sequences are already aligned to each
#      other via the input MSA, so a fresh pairwise realignment would just
#      re-derive what MAFFT already determined).
#
# Why anchor on HXB2 at all, given the samples are not subtype B?
# ---------------------------------------------------------------
# U3 position and content genuinely vary between subtypes, so this needed
# checking rather than assuming. Literature reviewed 2026-08-05:
#
#   * Mbondji-Wonje 2018 (PLoS One, 10.1371/journal.pone.0195661) sequenced U3R
#     across clades A1, B, C, D and F2 plus CRFs: the R region is "very well
#     conserved" across strains while U3 inter-strain dissimilarity reaches 25%.
#     Defining U3 as [LTR_start, R_start) and projecting HXB2's R boundary is
#     therefore sound; the U3 interior is what needs care.
#   * Jeeninga 2000 (J Virol, 10.1128/jvi.74.8.3740-3751.2000) found "a unique
#     LTR enhancer-promoter configuration for each subtype", NF-kB site count
#     ranging from one (subtype E) to three (subtype C); Naghavi 1999
#     (10.1089/088922299310197) showed subtype C's third site arises from an
#     INSERTION. Because the variation is partly insertional, a raw genome
#     coordinate slice cannot be correct -- which is why the earlier
#     `seqkit subseq -r 8677:9121` approach was abandoned. Lifting HXB2
#     coordinates through alignment columns follows insertions automatically.
#   * Parreira 2006 (Microbes Infect, 10.1016/j.micinf.2006.05.005) found only
#     63.3% of Mozambican subtype C viruses carried three NF-kB sites, so the
#     variation is not even fully subtype-determined and a per-subtype
#     coordinate table would be wrong for a third of that subtype's isolates.
#
# Consequently no per-subtype coordinate system exists to switch to (Los Alamos
# publishes HXB2-only coordinates, plus Mac239 for SIV), and adopting one would
# be incorrect anyway. Per-sequence, alignment-anchored derivation is the
# published method, not a workaround. Subtype variation is handled here as
# measured offsets relative to HXB2 rather than by swapping coordinate frames.
#
# Known limitation: the window ends at the column of HXB2's LAST U3 base, so a
# sample insertion sitting between HXB2's U3 end and HXB2's R start falls outside
# the window and is clipped, even though it is upstream of that sample's own R
# region and so belongs to its U3. Anchoring the end on HXB2's R-start column
# minus one would include it.
#
# Usage:
#   extract_u3_by_hxb2_anchor.sh \
#       --alignment aligned_consensus.fasta \
#       --hxb2-id K03455.1 \
#       --gb-cache reference/K03455.1.gb \
#       --out-gapped motifs/U3_aligned.fasta \
#       --out motifs/U3_extracted.fasta \
#       --warnings-log motifs/u3_extraction_warnings.log
# -u errors on unset vars, pipefail fails a pipe if any stage fails
set -uo pipefail

# NCBI efetch URL template to pull a GenBank record by accession
EUTILS_URL_TMPL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%s&rettype=gb&retmode=text"

# Length/identity tolerance for the outlier safety net.
# shortest plausible U3; anything below is flagged for review
MIN_U3_LEN=400
# longest plausible U3; anything above is flagged for review
MAX_U3_LEN=500
# minimum % identity to HXB2 U3 before a sequence is flagged
MIN_IDENTITY_PCT=70.0

# default HXB2 accession; overridable via --hxb2-id
HXB2_ID="K03455.1"

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
    # indirect-expand each name; empty = not supplied
    if [ -z "${!var:-}" ]; then
        echo "ERROR: missing required argument for ${var}" >&2  # say which one is missing
        exit 1                                       # and abort
    fi
done

if [ ! -s "${ALIGNMENT}" ]; then                     # the alignment must exist and be non-empty
    echo "ERROR: no alignment found at ${ALIGNMENT}" >&2  # otherwise there's nothing to slice
    exit 1
fi

# fetch HXB2's GenBank record only if not already cached
if [ ! -s "${GB_CACHE}" ]; then
    mkdir -p "$(dirname "${GB_CACHE}")"              # ensure the cache dir exists
    URL=$(printf "${EUTILS_URL_TMPL}" "${HXB2_ID}")  # fill the accession into the efetch URL
    # download the GenBank record quietly to the cache
    curl -s "${URL}" -o "${GB_CACHE}"
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
# Both LTR copies are parsed, not just the 5' one. In an integrated provirus the
# two LTRs are identical, so either is a valid U3 source -- but a given library
# may only have read coverage over one of them. Measured 2026-08-04 on the
# PacBio subset: 124_4 covers the 3' LTR at 100% and the 5' at 0%, while 203_3
# and 211_0 are the other way round. Taking the 5' copy unconditionally handed
# back reference-padded sequence for 124_4, whose extracted U3 came out
# byte-identical to HXB2 and produced HXB2's own motif hits. Each record now
# picks whichever copy it actually sequenced.
# capture the four 0-based coordinates awk prints below
read -r LTR5_START R5_START LTR3_START R3_START <<EOF
$(awk '
/^[[:space:]]*\// {                                  # a qualifier line (starts with optional whitespace then /)
    if ($0 ~ /\/note=/ && pending_type == "repeat_region") {  # a /note on the repeat_region feature we just saw
        note = $0                                    # copy the line to extract the note text
        sub(/.*\/note="/, "", note)                  # drop everything up to the opening quote
        sub(/".*/, "", note)                         # drop the closing quote and beyond, leaving the note text
        notelow = tolower(note)                      # lowercase for case-insensitive matching
        # HXB2 annotates exactly four relevant repeat_regions: the 5'/3' LTR
        # pair and the 5'/3' R-repeat pair, each note carrying its own "5" or "3"
        is5 = (index(notelow, "5") > 0)              # note refers to the 5-prime copy
        is3 = (index(notelow, "3") > 0)              # note refers to the 3-prime copy
        if (index(notelow, "ltr") > 0) {             # an LTR note
            if (is5 && ltr5_loc == "") ltr5_loc = pending_loc
            else if (is3 && ltr3_loc == "") ltr3_loc = pending_loc
        } else if (index(notelow, "r repeat") > 0) { # an R-repeat note
            if (is5 && r5_loc == "") r5_loc = pending_loc
            else if (is3 && r3_loc == "") r3_loc = pending_loc
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
    if (ltr5_loc == "" || r5_loc == "" || ltr3_loc == "" || r3_loc == "") {  # all four must have been found
        print "ERROR: could not find both LTR copies and both R repeat copies as repeat_region features" > "/dev/stderr"
        exit 1
    }
    split(ltr5_loc, a, "\\.\\."); split(r5_loc, b, "\\.\\.")   # split "start..end" of the 5-prime pair
    split(ltr3_loc, c, "\\.\\."); split(r3_loc, d, "\\.\\.")   # split "start..end" of the 3-prime pair
    print a[1]-1, b[1]-1, c[1]-1, d[1]-1              # all four starts converted 1-based -> 0-based
}
' "${GB_CACHE}")
EOF
for coord in LTR5_START R5_START LTR3_START R3_START; do  # awk must have produced all four
    if [ -z "${!coord:-}" ]; then
        # guide the user if parsing failed
        echo "ERROR: failed to parse LTR / R-repeat coordinates (${coord}) from ${GB_CACHE}. Inspect its FEATURES table and adjust the note-matching logic if this reference's annotation conventions differ from HXB2's." >&2
        exit 1
    fi
done
# U3 = [LTR_start, R_start) for each copy, so R must come after its own LTR start
if [ "${R5_START}" -le "${LTR5_START}" ] || [ "${R3_START}" -le "${LTR3_START}" ]; then
    # sanity-fail rather than emit nonsense coords
    echo "ERROR: an R-region start is not after its LTR start (5': ${LTR5_START}/${R5_START}, 3': ${LTR3_START}/${R3_START}) -- unexpected annotation, refusing to guess." >&2
    exit 1
fi
# log both resolved U3 windows and their lengths
echo "HXB2 U3 genome coordinates (0-based, half-open): 5' [${LTR5_START}, ${R5_START}) ($((R5_START - LTR5_START)) nt), 3' [${LTR3_START}, ${R3_START}) ($((R3_START - LTR3_START)) nt)" >&2

# ensure all three output dirs exist
mkdir -p "$(dirname "${OUT_GAPPED}")" "$(dirname "${OUT}")" "$(dirname "${WARNINGS_LOG}")"
# Per-record record of which LTR copy was used and how many bases each carried.
# Derived from --out rather than taking a new flag, so existing callers (the
# Illumina pipeline included) pick it up without changing their invocation.
LTR_REPORT="${OUT%.*}_ltr_choice.tsv"

# --- Liftover HXB2's genome coordinates to alignment columns, slice every
# record, strip gaps, and run the outlier safety net -- all in one pass
# over the buffered alignment. ---
awk -v hxb2_id="${HXB2_ID}" \
    -v u3s5="${LTR5_START}" -v u3e5="${R5_START}" -v u3s3="${LTR3_START}" -v u3e3="${R3_START}" \
    -v min_len="${MIN_U3_LEN}" -v max_len="${MAX_U3_LEN}" -v min_identity="${MIN_IDENTITY_PCT}" \
    -v out_gapped="${OUT_GAPPED}" -v out="${OUT}" -v warnings_log="${WARNINGS_LOG}" \
    -v ltr_report="${LTR_REPORT}" '  # pass all resolved coords/thresholds/paths into awk
# count bases that were actually sequenced (N and gaps are not evidence)
function n_called(s,   i, ch, c) {
    c = 0
    for (i = 1; i <= length(s); i++) {
        ch = toupper(substr(s, i, 1))
        if (ch == "A" || ch == "C" || ch == "G" || ch == "T") c++
    }
    return c
}
# strip alignment gaps from a slice to give the record own ungapped sequence
function degap(s) { gsub(/-/, "", s); return s }
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
    cs5 = 0; ce5 = 0; cs3 = 0; ce3 = 0                # alignment columns for both U3 windows (found below)
    for (c = 1; c <= L; c++) {                        # walk HXB2 column by column to lift genome coords to columns
        if (substr(hxb2_seq, c, 1) != "-") {          # only real bases advance the genome coordinate
            if (k == u3s5) cs5 = c                    # column where the 5-prime U3 start base sits
            if (k == u3e5 - 1) ce5 = c                # column of the last base before the 5-prime R
            if (k == u3s3) cs3 = c                    # column where the 3-prime U3 start base sits
            if (k == u3e3 - 1) ce3 = c                # column of the last base before the 3-prime R
            k++                                       # one more ungapped base consumed
        }
    }
    if (cs5 == 0 || ce5 == 0 || cs3 == 0 || ce3 == 0) {  # all four boundaries must have been located
        print "ERROR: HXB2 aligned row has only " k " ungapped bases, fewer than the 3-prime U3 end coordinate (" u3e3 "). The HXB2 record in this alignment looks truncated or wrong." > "/dev/stderr"
        exit 1
    }
    print "Alignment-column ranges for U3 (0-based, inclusive): 5-prime [" cs5-1 ", " ce5-1 "], 3-prime [" cs3-1 ", " ce3-1 "]" > "/dev/stderr"

    # HXB2 own U3 from each window. Needed up front as the per-window
    # denominator below, and reused later as the identity reference.
    hxb2_g5 = substr(seqs[hxb2_idx], cs5, ce5 - cs5 + 1)
    hxb2_g3 = substr(seqs[hxb2_idx], cs3, ce3 - cs3 + 1)
    den5 = n_called(degap(hxb2_g5))                   # fully-sequenced 5-prime U3 length (453 for HXB2)
    den3 = n_called(degap(hxb2_g3))                   # fully-sequenced 3-prime U3 length (454 for HXB2)
    if (den5 == 0 || den3 == 0) {                     # HXB2 must carry both copies to anchor either
        print "ERROR: HXB2 record has no sequenced bases in one of its LTR copies -- cannot anchor U3." > "/dev/stderr"
        exit 1
    }

    printf "" > out_gapped                            # truncate/create the gapped output file
    printf "" > out                                   # truncate/create the ungapped output file
    printf "record\tltr_used\tcalled_5prime\tcalled_3prime\n" > ltr_report  # header for the choice report
    n_warnings = 0                                    # count of flagged sequences
    n_three = 0                                       # how many records ended up using the 3-prime copy

    for (i = 1; i <= n; i++) {                        # slice BOTH U3 windows out of every record
        g5 = substr(seqs[i], cs5, ce5 - cs5 + 1)      # gapped 5-prime U3 slice
        g3 = substr(seqs[i], cs3, ce3 - cs3 + 1)      # gapped 3-prime U3 slice
        c5 = n_called(degap(g5))                      # sequenced bases in the 5-prime copy
        c3 = n_called(degap(g3))                      # sequenced bases in the 3-prime copy
        # Both LTRs are identical in an integrated provirus, so prefer whichever
        # copy this record actually sequenced. Compared as a FRACTION of each
        # window, not as raw counts: HXB2 annotates the two U3s at 453 and 454nt,
        # so a raw count would hand every fully-covered record to the 3-prime
        # copy on a one-base technicality. Ties keep the 5-prime copy, which
        # preserves the previous behaviour for fully-covered inputs.
        if (c3 / den3 > c5 / den5) { gslice[i] = g3; used[i] = "3-prime"; hxb2_ref_win[i] = "3"; n_three++ }
        else                       { gslice[i] = g5; used[i] = "5-prime"; hxb2_ref_win[i] = "5" }
        uslice_arr[i] = degap(gslice[i])              # this record own ungapped U3 sequence
        printf "%s\t%s\t%d\t%d\n", ids[i], used[i], c5, c3 >> ltr_report

        print ">" ids[i] >> out_gapped                # write gapped record header
        print toupper(gslice[i]) >> out_gapped        # write gapped U3 (uppercased)
        print ">" ids[i] >> out                        # write ungapped record header
        print toupper(uslice_arr[i]) >> out           # write ungapped U3 (uppercased)
    }
    close(out_gapped)                                 # flush/close so file sizes are final
    close(out)
    close(ltr_report)

    for (i = 1; i <= n; i++) {                        # outlier safety net over every extracted U3
        flags = ""                                    # accumulate reasons this sequence looks off
        seqlen = length(uslice_arr[i])                # its ungapped U3 length
        if (seqlen < min_len || seqlen > max_len) {   # length outside the tolerance band?
            flags = flags "; length " seqlen "nt outside [" min_len "," max_len "]"  # note it
        }
        # a record with no sequenced base in either LTR carries no U3 evidence at all
        if (n_called(uslice_arr[i]) == 0) {
            flags = flags "; no sequenced bases in either LTR (U3 not recovered)"
        }
        hxb2_gapped = (hxb2_ref_win[i] == "3") ? hxb2_g3 : hxb2_g5   # match the window this record used
        hxb2_ungapped_len = length(degap(hxb2_gapped))               # denominator for the identity percentage
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
    print "Wrote " n " U3 sequences to " out " (" n_three " used the 3-prime LTR copy; per-record detail in " ltr_report ")" > "/dev/stderr"  # final count written
}
' "${ALIGNMENT}"                                     # run the whole liftover/slice/check pass over the alignment
