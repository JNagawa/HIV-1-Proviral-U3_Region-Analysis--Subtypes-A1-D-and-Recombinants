#!/bin/bash
# Downloads JASPAR 2024 PWMs for the six TFs the proposal targets (NF-kB
# p65/p50, SP1, NFAT, AP-1, TBP -- standing in for TFIID, which has no
# dedicated JASPAR PWM since it's a large basal multi-subunit complex, not
# a sequence-specific DNA-binding factor). Matrix IDs below were looked up
# via JASPAR's own REST API (https://jaspar.elixir.no/api/v1/matrix/?name=...),
# not guessed:
#   NF-kB p65 (RELA):  MA0107.1
#   NF-kB p50 (NFKB1): MA0105.4
#   SP1:               MA0079.5
#   NFAT (NFATC2):     MA0152.3
#   AP-1 (FOS::JUN):   MA0099.4
#   TBP:               MA0108.3
# -e exit on error, -u error on unset vars, pipefail catch pipeline failures
set -euo pipefail
# absolute path of this script's dir, so it runs from anywhere
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
# repo top-level (two dirs up), base for the reference path
REPO_ROOT="$(cd "${STAGE_DIR}/../.." && pwd)"
mkdir -p "${REPO_ROOT}/data/reference/jaspar"         # ensure the JASPAR download dir exists
# work inside it so all outputs land there with plain names
cd "${REPO_ROOT}/data/reference/jaspar"

# the six core-TF JASPAR matrix IDs (see header)
MATRIX_IDS="MA0107.1 MA0105.4 MA0079.5 MA0152.3 MA0099.4 MA0108.3"

# truncate/create the combined MEME file before appending
> core6_pfms.meme
# flag: first matrix keeps the MEME header, rest are appended headerless
FIRST=1
for id in ${MATRIX_IDS}; do                           # download each matrix in MEME format
    # first matrix: take the whole response (includes the MEME header)
    if [ "${FIRST}" = "1" ]; then
        # fetch and write full MEME (header + motif)
        curl -s -H "Accept: text/meme" "https://jaspar.elixir.no/api/v1/matrix/${id}/" > core6_pfms.meme
        FIRST=0                                       # subsequent matrices skip the header
    else
        # append only the MOTIF block (skip the repeated MEME header) for
        # subsequent matrices, so the file is one valid multi-motif MEME file
        # print from the first MOTIF line onward, dropping the duplicate header
        curl -s -H "Accept: text/meme" "https://jaspar.elixir.no/api/v1/matrix/${id}/" | \
            awk '/^MOTIF/{p=1} p' >> core6_pfms.meme
    fi
done
# report how many motifs landed in the combined file
echo "Wrote $(grep -c '^MOTIF' core6_pfms.meme) motifs to data/reference/jaspar/core6_pfms.meme (for FIMO)"

# Raw JASPAR-format flat file (for TFBSTools::readJASPARMatrix and, split
# per-matrix below, for MOODS).
# download the full CORE vertebrates flat file; retry without the Accept header if the first attempt
# fails
curl -s -H "Accept: text/plain" \
    "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt" \
    -o all_vertebrates_pfms.jaspar 2>/dev/null || \
curl -s "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt" \
    -o all_vertebrates_pfms.jaspar
# report the total motif count (>-lines) downloaded
echo "Wrote all_vertebrates_pfms.jaspar ($(grep -c '^>' all_vertebrates_pfms.jaspar) motifs, full CORE vertebrates set)"

# Split the 6 core matrices into individual .pfm files (4 lines of counts,
# no header) for MOODS, which expects one matrix per file. JASPAR raw format
# is repeating blocks: a ">ID name" header line followed by 4 count rows,
# one per base, shaped like "A  [ 1 2 3 ... ]".
awk -v ids="MA0107.1 MA0105.4 MA0079.5 MA0152.3 MA0099.4 MA0108.3" '  # pass the 6 wanted IDs into awk to extract them from the full set
BEGIN { n_ids = split(ids, idlist, " ") }             # split the ID string into an array we can iterate at the end
/^>/ {                                                # a ">ID name" header line starts a new matrix block
    curid = $1                                        # first field is ">ID"
    sub(/^>/, "", curid)                              # strip the leading ">" to get the bare matrix ID
    headers[curid] = $0                               # remember the full header line, keyed by ID
    counts[curid] = 0                                 # reset the count-row counter for this matrix
    next                                              # header consumed, move to next line
}
{
    counts[curid]++                                   # this is a count row for the current matrix
    lines[curid, counts[curid]] = $0                  # store it, indexed by matrix ID and row number
}
END {
    for (i = 1; i <= n_ids; i++) {                    # write one .pfm file per wanted matrix
        mid = idlist[i]                               # current wanted matrix ID
        if (!(mid in headers)) {                      # matrix not present in the downloaded set...
            print "WARNING: " mid " not found in all_vertebrates_pfms.jaspar"  # ...warn so the gap is visible
            continue                                  # ...and skip it
        }
        pfmfile = mid ".pfm"                           # per-matrix output filename MOODS expects
        printf "" > pfmfile                            # truncate/create it before writing
        for (k = 1; k <= counts[mid]; k++) {           # emit each of the 4 base count rows
            row = lines[mid, k]                        # the stored raw row, shaped like "A  [ 1 2 3 ]"
            open_bracket = index(row, "[")             # find the "[" that opens the count list
            close_bracket = index(row, "]")            # find the closing "]"
            row_counts = substr(row, open_bracket + 1, close_bracket - open_bracket - 1)  # pull out just the numbers between the brackets
            gsub(/^[ \t]+|[ \t]+$/, "", row_counts)    # trim leading/trailing whitespace
            print row_counts >> pfmfile                # append the bare count row (MOODS wants no header/base labels)
        }
        close(pfmfile)                                 # close so the file is flushed before ls counts it
    }

    combined = "core6_pfms.jaspar"                     # single JASPAR flat file holding just the 6 matrices, for TFBSTools
    printf "" > combined                               # truncate/create it
    for (i = 1; i <= n_ids; i++) {                     # re-emit each wanted matrix in full JASPAR format
        mid = idlist[i]                                # current wanted matrix ID
        if (!(mid in headers)) continue                # skip any that were missing from the download
        print headers[mid] >> combined                 # write the ">ID name" header
        for (k = 1; k <= counts[mid]; k++) print lines[mid, k] >> combined  # then its 4 original count rows verbatim
    }
    close(combined)                                    # flush the combined file
}
' all_vertebrates_pfms.jaspar                          # feed the full downloaded flat file into the awk program above
# report how many per-matrix .pfm files were produced
echo "Wrote per-matrix .pfm files for MOODS: $(ls MA*.pfm 2>/dev/null | wc -l) files"
# confirm the combined TFBSTools input was written
echo "Wrote data/reference/jaspar/core6_pfms.jaspar for TFBSTools"
