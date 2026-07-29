#!/bin/bash
# Extract the proviral core from HIV-SMRTcap reads/consensus sequences whose
# host-genome flanks have been masked with N (the SMRTcap convention: host is
# N-masked, integrated provirus stays in ACGT). The virus integrates at a
# different host position in every read, so the flanks are of variable length
# and sit at BOTH ends, e.g.
#     NNNNNACGGTTCTAGGGTTTCCACTANNNNN  ->  ACGGTTCTAGGGTTTCCACTA
#     NNNCCGTATGCGCCCTAANNNNNN         ->  CCGTATGCGCCCTAA
# so the extraction is simply: strip the leading N-run and the trailing N-run,
# keep the ACGT core. Only the anchored end-runs are removed -- any internal
# bases (including a stray ambiguous base inside the provirus) are preserved,
# so a genuine internal N never truncates the genome.
#
# Handles BOTH FASTA and FASTQ input, detected from the first record rather
# than the filename (the SMRTcap files are named "*.fastq.hiv.unmasked.fa" but
# hold FASTA, so the extension cannot be trusted). For FASTQ the quality string
# is sliced with exactly the same coordinates as the sequence, so bases and
# their qualities stay in register and the output is a valid FASTQ.
#
# Bash/awk/seqkit only (no Python), matching this repo's convention. seqkit is
# used purely to linearise records (fx2tab) and re-wrap them (tab2fx); all the
# stripping logic is awk.
#
# Usage: extract_provirus_strip_hostN.sh <in.fast[aq][.gz]> <out_provirus> [coords.tsv]
#   in            host-N-masked reads/consensus (FASTA or FASTQ, plain or .gz)
#   out_provirus  proviral cores, one record per input record that had a
#                 non-empty core (all-N / empty records are dropped); written
#                 in the SAME format as the input
#   coords.tsv    optional provenance log; one row per INPUT record:
#                 name  orig_len  lead_N  trail_N  provirus_len  prov_start  prov_end  status
#                 (prov_start/prov_end are 1-based, inclusive, in original
#                 read coordinates; status = kept | dropped_all_N | dropped_empty)
# -u errors on unset vars, pipefail fails a pipe if any stage fails
set -uo pipefail

# arg 1 = host-N-masked input; :? prints usage and aborts if missing
IN="${1:?usage: extract_provirus_strip_hostN.sh <in.fast[aq][.gz]> <out> [coords.tsv]}"
OUT="${2:?missing output path}"                      # arg 2 = where to write the proviral cores
# arg 3 = optional coords/provenance TSV (empty if not given)
COORDS="${3:-}"

if [ ! -s "${IN}" ]; then                            # nothing to do without a non-empty input...
    echo "ERROR: input '${IN}' not found or empty." >&2  # ...report the problem...
    exit 1                                           # ...and fail
fi
# seqkit does the record linearise/re-wrap, so it must be present
if ! command -v seqkit >/dev/null 2>&1; then
    # tell the user how to get it
    echo "ERROR: seqkit not on PATH (conda activate HIV_U3analysis)." >&2
    exit 1                                           # fail if it's missing
fi

mkdir -p "$(dirname "${OUT}")"                       # make sure the output dir exists
# and the coords dir too, only if a coords path was given
[ -n "${COORDS}" ] && mkdir -p "$(dirname "${COORDS}")"
# truncate any coords file from a previous run: the awk below APPENDS rows, so
# without this a rerun would stack new rows onto the old ones
[ -n "${COORDS}" ] && : > "${COORDS}"

# Detect the record format from the first non-blank character of the data
# itself ('>' = FASTA, '@' = FASTQ). zcat -f reads plain and gzipped files
# alike, so this works for .gz input without a separate branch.
FIRST_CHAR=$(zcat -f "${IN}" 2>/dev/null | awk 'NF{print substr($0,1,1); exit}')
case "${FIRST_CHAR}" in
    '>') FORMAT="fasta" ;;
    '@') FORMAT="fastq" ;;
    *)
        # anything else is not sequence data we can strip
        echo "ERROR: '${IN}' does not start with '>' or '@' -- not FASTA or FASTQ." >&2
        exit 1 ;;
esac

# fx2tab emits "name<TAB>sequence" for FASTA. Adding -q appends the quality
# string, and then an average-quality column that we must NOT pass on -- the
# awk below prints exactly 2 or 3 fields so tab2fx re-wraps the right format.
# -w0 (in tab2fx below) disables line wrapping so downstream length checks are
# unambiguous. Case-insensitive N-stripping ([Nn]) covers masks written in
# either case. The core is bases [lead+1 .. len-trail] of the original read.
if [ "${FORMAT}" = "fastq" ]; then
    FX2TAB_OPTS="-q"                                 # carry the quality string through
else
    FX2TAB_OPTS=""                                   # FASTA has no quality column
fi

# compress the output when the caller asks for a .gz path, so the QC harness can
# hand this straight to tools that expect gzipped reads
if [ "${OUT%.gz}" != "${OUT}" ]; then
    OUT_CMD="gzip -c"                                # OUT ends in .gz
else
    OUT_CMD="cat"                                    # plain-text output
fi

# linearise each record, strip the N-flanks in awk, then re-wrap to the input format
# shellcheck disable=SC2086
zcat -f "${IN}" 2>/dev/null | seqkit fx2tab ${FX2TAB_OPTS} 2>/dev/null | awk -F'\t' -v coords="${COORDS}" -v fmt="${FORMAT}" '
    {
        name = $1                    # record name (first tab field)
        seq  = $2                    # linearised sequence (second tab field)
        qual = (fmt == "fastq" ? $3 : "")  # quality string, FASTQ only
        orig = length(seq)           # original read length before stripping

        # leading N run
        lead = 0                     # count of leading N bases (the 5-prime host flank)
        while (lead < orig && substr(seq, lead+1, 1) ~ /[Nn]/) lead++  # advance past every leading N (case-insensitive)

        if (lead == orig) {            # sequence is entirely N -> no provirus
            if (coords != "")        # log it as dropped only if a coords file was requested
                printf "%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n", name, orig, orig, 0, 0, "NA", "NA", "dropped_all_N" >> coords  # provenance row for an all-N read
            next                     # emit no core, move to next record
        }

        # trailing N run
        trail = 0                    # count of trailing N bases (the 3-prime host flank)
        while (trail < orig && substr(seq, orig-trail, 1) ~ /[Nn]/) trail++  # advance inward past every trailing N

        core_len   = orig - lead - trail  # length of the ACGT proviral core left after removing both flanks
        prov_start = lead + 1             # 1-based start of the core in original coordinates
        prov_end   = orig - trail         # 1-based inclusive end of the core in original coordinates
        core       = substr(seq, prov_start, core_len)  # the proviral core itself

        if (core_len <= 0) {           # defensive: nothing left after stripping
            if (coords != "")        # log the empty result if coords requested
                printf "%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n", name, orig, lead, trail, 0, "NA", "NA", "dropped_empty" >> coords  # provenance row for an empty core
            next                     # emit no core
        }

        # emit the core for tab2fx to re-wrap; the quality string is cut with the
        # SAME start/length as the sequence so bases and qualities stay aligned
        if (fmt == "fastq")
            printf "%s\t%s\t%s\n", name, core, substr(qual, prov_start, core_len)  # name<TAB>seq<TAB>qual -> FASTQ
        else
            printf "%s\t%s\n", name, core  # name<TAB>seq -> FASTA
        if (coords != "")            # and record its provenance if coords requested
            printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", name, orig, lead, trail, core_len, prov_start, prov_end, "kept" >> coords  # provenance row for a kept core
    }
' | seqkit tab2fx -w0 2>/dev/null | ${OUT_CMD} > "${OUT}"  # re-wrap back into the input format (-w0 = no line wrapping), gzipping if asked

# Prepend the coords header (awk appended rows without one, so it stays
# valid even when run per-record above).
# only if a coords file was requested and actually got written
if [ -n "${COORDS}" ] && [ -f "${COORDS}" ]; then
    # the column header for the coords TSV
    HDR="name\torig_len\tlead_N\ttrail_N\tprovirus_len\tprov_start\tprov_end\tstatus"
    TMP_C="$(mktemp)"                                # scratch file to prepend the header
    # write header then existing rows, then swap it in
    { printf "%b\n" "${HDR}"; cat "${COORDS}"; } > "${TMP_C}" && mv "${TMP_C}" "${COORDS}"
fi

# number of input records (num_seqs column)
N_IN=$(seqkit stats -T "${IN}" 2>/dev/null | awk -F'\t' 'NR==2{print $4}')
# number of output records that kept a core
N_OUT=$(seqkit stats -T "${OUT}" 2>/dev/null | awk -F'\t' 'NR==2{print $4}')
# summary line to stderr
echo "extract_provirus_strip_hostN: ${N_OUT:-0}/${N_IN:-0} records had a non-empty proviral core (${FORMAT}) -> ${OUT}" >&2
