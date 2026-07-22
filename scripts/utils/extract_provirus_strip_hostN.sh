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
set -uo pipefail

IN="${1:?usage: extract_provirus_strip_hostN.sh <in.fasta[.gz]> <out.fasta> [coords.tsv]}"
OUT="${2:?missing output FASTA path}"
COORDS="${3:-}"

if [ ! -s "${IN}" ]; then
    echo "ERROR: input '${IN}' not found or empty." >&2
    exit 1
fi
if ! command -v seqkit >/dev/null 2>&1; then
    echo "ERROR: seqkit not on PATH (conda activate HIV_U3analysis)." >&2
    exit 1
fi

mkdir -p "$(dirname "${OUT}")"
[ -n "${COORDS}" ] && mkdir -p "$(dirname "${COORDS}")"

# fx2tab emits "name<TAB>sequence" (one line per record, sequence linearised).
# -w0 (in tab2fx below) disables line wrapping so downstream length checks are
# unambiguous. Case-insensitive N-stripping ([Nn]) covers masks written in
# either case. The core is bases [lead+1 .. len-trail] of the original read.
TMP_TAB="$(mktemp)"
trap 'rm -f "${TMP_TAB}"' EXIT

seqkit fx2tab "${IN}" 2>/dev/null | awk -F'\t' -v coords="${COORDS}" '
    {
        name = $1
        seq  = $2
        orig = length(seq)

        # leading N run
        lead = 0
        while (lead < orig && substr(seq, lead+1, 1) ~ /[Nn]/) lead++

        if (lead == orig) {            # sequence is entirely N -> no provirus
            if (coords != "")
                printf "%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n", name, orig, orig, 0, 0, "NA", "NA", "dropped_all_N" >> coords
            next
        }

        # trailing N run
        trail = 0
        while (trail < orig && substr(seq, orig-trail, 1) ~ /[Nn]/) trail++

        core_len   = orig - lead - trail
        prov_start = lead + 1
        prov_end   = orig - trail
        core       = substr(seq, prov_start, core_len)

        if (core_len <= 0) {           # defensive: nothing left after stripping
            if (coords != "")
                printf "%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n", name, orig, lead, trail, 0, "NA", "NA", "dropped_empty" >> coords
            next
        }

        # emit the core as name<TAB>seq for tab2fx to re-wrap into FASTA
        printf "%s\t%s\n", name, core
        if (coords != "")
            printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", name, orig, lead, trail, core_len, prov_start, prov_end, "kept" >> coords
    }
' | seqkit tab2fx -w0 > "${OUT}" 2>/dev/null

# Prepend the coords header (awk appended rows without one, so it stays
# valid even when run per-record above).
if [ -n "${COORDS}" ] && [ -f "${COORDS}" ]; then
    HDR="name\torig_len\tlead_N\ttrail_N\tprovirus_len\tprov_start\tprov_end\tstatus"
    TMP_C="$(mktemp)"
    { printf "%b\n" "${HDR}"; cat "${COORDS}"; } > "${TMP_C}" && mv "${TMP_C}" "${COORDS}"
fi

N_IN=$(seqkit stats -T "${IN}" 2>/dev/null | awk -F'\t' 'NR==2{print $4}')
N_OUT=$(seqkit stats -T "${OUT}" 2>/dev/null | awk -F'\t' 'NR==2{print $4}')
echo "extract_provirus_strip_hostN: ${N_OUT:-0}/${N_IN:-0} records had a non-empty proviral core -> ${OUT}" >&2
