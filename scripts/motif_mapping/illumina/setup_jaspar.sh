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
set -euo pipefail
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../.." && pwd)"
mkdir -p "${REPO_ROOT}/data/reference/jaspar"
cd "${REPO_ROOT}/data/reference/jaspar"

MATRIX_IDS="MA0107.1 MA0105.4 MA0079.5 MA0152.3 MA0099.4 MA0108.3"

> core6_pfms.meme
FIRST=1
for id in ${MATRIX_IDS}; do
    if [ "${FIRST}" = "1" ]; then
        curl -s -H "Accept: text/meme" "https://jaspar.elixir.no/api/v1/matrix/${id}/" > core6_pfms.meme
        FIRST=0
    else
        # append only the MOTIF block (skip the repeated MEME header) for
        # subsequent matrices, so the file is one valid multi-motif MEME file
        curl -s -H "Accept: text/meme" "https://jaspar.elixir.no/api/v1/matrix/${id}/" | \
            awk '/^MOTIF/{p=1} p' >> core6_pfms.meme
    fi
done
echo "Wrote $(grep -c '^MOTIF' core6_pfms.meme) motifs to data/reference/jaspar/core6_pfms.meme (for FIMO)"

# Raw JASPAR-format flat file (for TFBSTools::readJASPARMatrix and, split
# per-matrix below, for MOODS).
curl -s -H "Accept: text/plain" \
    "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt" \
    -o all_vertebrates_pfms.jaspar 2>/dev/null || \
curl -s "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt" \
    -o all_vertebrates_pfms.jaspar
echo "Wrote all_vertebrates_pfms.jaspar ($(grep -c '^>' all_vertebrates_pfms.jaspar) motifs, full CORE vertebrates set)"

# Split the 6 core matrices into individual .pfm files (4 lines of counts,
# no header) for MOODS, which expects one matrix per file. JASPAR raw format
# is repeating blocks: a ">ID name" header line followed by 4 count rows,
# one per base, shaped like "A  [ 1 2 3 ... ]".
awk -v ids="MA0107.1 MA0105.4 MA0079.5 MA0152.3 MA0099.4 MA0108.3" '
BEGIN { n_ids = split(ids, idlist, " ") }
/^>/ {
    curid = $1
    sub(/^>/, "", curid)
    headers[curid] = $0
    counts[curid] = 0
    next
}
{
    counts[curid]++
    lines[curid, counts[curid]] = $0
}
END {
    for (i = 1; i <= n_ids; i++) {
        mid = idlist[i]
        if (!(mid in headers)) {
            print "WARNING: " mid " not found in all_vertebrates_pfms.jaspar"
            continue
        }
        pfmfile = mid ".pfm"
        printf "" > pfmfile
        for (k = 1; k <= counts[mid]; k++) {
            row = lines[mid, k]
            open_bracket = index(row, "[")
            close_bracket = index(row, "]")
            row_counts = substr(row, open_bracket + 1, close_bracket - open_bracket - 1)
            gsub(/^[ \t]+|[ \t]+$/, "", row_counts)
            print row_counts >> pfmfile
        }
        close(pfmfile)
    }

    combined = "core6_pfms.jaspar"
    printf "" > combined
    for (i = 1; i <= n_ids; i++) {
        mid = idlist[i]
        if (!(mid in headers)) continue
        print headers[mid] >> combined
        for (k = 1; k <= counts[mid]; k++) print lines[mid, k] >> combined
    }
    close(combined)
}
' all_vertebrates_pfms.jaspar
echo "Wrote per-matrix .pfm files for MOODS: $(ls MA*.pfm 2>/dev/null | wc -l) files"
echo "Wrote data/reference/jaspar/core6_pfms.jaspar for TFBSTools"
