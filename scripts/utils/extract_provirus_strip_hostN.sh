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
# Bash/awk/seqkit only (no Python), matching this repo's convention. seqkit is
# used purely to linearise multi-line FASTA (fx2tab) and re-wrap it (tab2fx);
# all the stripping logic is awk.
#
# Usage: extract_provirus_strip_hostN.sh <in.fasta[.gz]> <out_provirus.fasta> [coords.tsv]
#   in.fasta      host-N-masked reads/consensus (FASTA or FASTA.gz)
#   out_provirus  proviral cores, one record per input record that had a
#                 non-empty core (all-N / empty records are dropped)
#   coords.tsv    optional provenance log; one row per INPUT record:
#                 name  orig_len  lead_N  trail_N  provirus_len  prov_start  prov_end  status
#                 (prov_start/prov_end are 1-based, inclusive, in original
#                 read coordinates; status = kept | dropped_all_N | dropped_empty)
# -u errors on unset vars, pipefail fails a pipe if any stage fails
set -uo pipefail

# arg 1 = host-N-masked input; :? prints usage and aborts if missing
IN="${1:?usage: extract_provirus_strip_hostN.sh <in.fasta[.gz]> <out.fasta> [coords.tsv]}"
OUT="${2:?missing output FASTA path}"                # arg 2 = where to write the proviral cores
# arg 3 = optional coords/provenance TSV (empty if not given)
COORDS="${3:-}"

if [ ! -s "${IN}" ]; then                            # nothing to do without a non-empty input...
    echo "ERROR: input '${IN}' not found or empty." >&2  # ...report the problem...
    exit 1                                           # ...and fail
fi
# seqkit does the FASTA linearise/re-wrap, so it must be present
if ! command -v seqkit >/dev/null 2>&1; then
    # tell the user how to get it
    echo "ERROR: seqkit not on PATH (conda activate HIV_U3analysis)." >&2
    exit 1                                           # fail if it's missing
fi

mkdir -p "$(dirname "${OUT}")"                       # make sure the output dir exists
# and the coords dir too, only if a coords path was given
[ -n "${COORDS}" ] && mkdir -p "$(dirname "${COORDS}")"

# fx2tab emits "name<TAB>sequence" (one line per record, sequence linearised).
# -w0 (in tab2fx below) disables line wrapping so downstream length checks are
# unambiguous. Case-insensitive N-stripping ([Nn]) covers masks written in
# either case. The core is bases [lead+1 .. len-trail] of the original read.
# scratch file (reserved for temp work; cleaned on exit)
TMP_TAB="$(mktemp)"
# always remove the temp file when the script exits
trap 'rm -f "${TMP_TAB}"' EXIT

seqkit fx2tab "${IN}" 2>/dev/null | awk -F'\t' -v coords="${COORDS}" '  # linearise each record to name<TAB>seq, then strip N-flanks in awk
    {
        name = $1                    # record name (first tab field)
        seq  = $2                    # linearised sequence (second tab field)
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

        # emit the core as name<TAB>seq for tab2fx to re-wrap into FASTA
        printf "%s\t%s\n", name, core  # pass the kept core downstream to tab2fx
        if (coords != "")            # and record its provenance if coords requested
            printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", name, orig, lead, trail, core_len, prov_start, prov_end, "kept" >> coords  # provenance row for a kept core
    }
' | seqkit tab2fx -w0 > "${OUT}" 2>/dev/null         # re-wrap name<TAB>seq back into FASTA (-w0 = no line wrapping) as the output

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
echo "extract_provirus_strip_hostN: ${N_OUT:-0}/${N_IN:-0} records had a non-empty proviral core -> ${OUT}" >&2
