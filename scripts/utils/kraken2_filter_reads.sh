#!/bin/bash
# Classify paired-end reads with Kraken2 against a database containing
# human + bacterial + viral genomes, then filter OUT any read whose
# assigned taxID descends from Homo sapiens (9606) or Bacteria (taxid 2),
# keeping unclassified reads and reads classified as viral (retaining real
# HIV-1 signal). Descent is resolved via the database's own bundled
# nodes.dmp with a bash/awk ancestor walk -- no external taxonomy
# tool (e.g. KrakenTools, which is Python) needed.
#
# A viral-only database can't do this: it has no human/bacterial genomes to
# match against, so it can only confirm "this read matches a known virus",
# not detect host contamination -- the dominant contamination source in
# HIV proviral sequencing. This is why a size-capped Standard database
# (bacteria+archaea+viral+human, ~16GB) is used instead of the much
# smaller viral-only option.
#
# Usage: kraken2_filter_reads.sh <R1.fastq.gz> <R2.fastq.gz> <OUTDIR> <SAMPLE_NAME> <KRAKEN2_DB_DIR>
# Produces:
#   <OUTDIR>/<SAMPLE_NAME>_1.kraken_filtered.fastq.gz
#   <OUTDIR>/<SAMPLE_NAME>_2.kraken_filtered.fastq.gz
#   <OUTDIR>/<SAMPLE_NAME>.kreport   (Kraken2's own hierarchical report)
# -u errors on unset vars, pipefail fails a pipe if any stage fails
set -uo pipefail
# positional args: R1/R2 inputs, output dir, sample name, Kraken2 DB dir
R1="$1" R2="$2" OUTDIR="$3" SAMPLE="$4" KRAKEN2_DB="$5"

# taxonomy tree bundled with the DB, used for the ancestor walk
NODES_DMP="${KRAKEN2_DB}/nodes.dmp"
# bail if the DB isn't a real extracted Kraken2 database
[ -s "${NODES_DMP}" ] || { echo "ERROR: ${NODES_DMP} not found -- is KRAKEN2_DB (${KRAKEN2_DB}) an extracted Kraken2 database?" >&2; exit 1; }

mkdir -p "${OUTDIR}"                                  # ensure the output dir exists

# classify read pairs against the DB; per-pair output + hierarchical report
kraken2 --db "${KRAKEN2_DB}" --paired --gzip-compressed --threads "${THREADS:-4}" \
    --output "${OUTDIR}/${SAMPLE}.kraken" \
    --report "${OUTDIR}/${SAMPLE}.kreport" \
    "${R1}" "${R2}"

# nodes.dmp fields are separated by "\t|\t" (taxid | parent_taxid | rank | ...);
# Kraken2's own --output is plain-tab-separated (C/U, read_id, taxid,
# length, LCA mapping). These two files need DIFFERENT field separators --
# a single global FS (via -F) applied to both was a real bug caught during
# testing: it silently left every kraken-output field empty, so no read was
# ever excluded and the "keep" list was full of blank IDs. Explicit
# per-file split() calls avoid that.
#
# First file (nodes.dmp) builds the taxid->parent map; second file is
# walked per read, climbing ancestors until root (1), Bacteria (2), or
# Human (9606) is hit. Unclassified reads (taxid 0) have no entry in
# parent[], so the walk stops immediately and they're kept, same as any
# non-human/non-bacterial classified read (e.g. Viruses, Archaea).
awk '
    NR==FNR {                                         # first file = nodes.dmp: build the taxid->parent map
        split($0, f, "\t\\|\t")                       # nodes.dmp fields are separated by "\t|\t"
        parent[f[1]] = f[2]                           # parent[taxid] = parent_taxid
        next
    }
    {
        split($0, f, "\t")                            # second file = kraken output, plain-tab separated
        read_id = f[2]                                # field 2 = read id
        t = f[3]                                      # field 3 = assigned taxid
        keep = 1                                       # keep unless we hit human/bacteria
        depth = 0                                      # loop guard against cycles/broken trees
        while (t != "1" && t != "" && depth < 50) {   # climb ancestors until root (1) or missing
            if (t == "9606" || t == "2") { keep = 0; break }  # descends from Homo sapiens (9606) or Bacteria (2) -> drop
            if (!(t in parent)) break                 # no parent (e.g. unclassified taxid 0) -> stop, keep
            t = parent[t]                             # step up to the parent taxid
            depth++                                   # count the step
        }
        if (keep) print read_id                       # emit ids of pairs to retain
    }
' "${NODES_DMP}" "${OUTDIR}/${SAMPLE}.kraken" > "${OUTDIR}/${SAMPLE}.keep_read_ids.txt"  # write the keep-list of read ids

N_TOTAL=$(wc -l < "${OUTDIR}/${SAMPLE}.kraken")              # total classified read pairs
# pairs kept after host/bacterial removal
N_KEEP=$(wc -l < "${OUTDIR}/${SAMPLE}.keep_read_ids.txt")
# progress summary to stderr
echo "Kraken2 filtering for ${SAMPLE}: ${N_KEEP}/${N_TOTAL} read pairs retained (human/bacterial reads discarded)." >&2

# No -n/--by-name here: seqkit's default match target is the ID (text
# before the first space), which is exactly what keep_read_ids.txt
# contains -- -n switches to matching the FULL header line instead (an
# exact whole-string match per seqkit's own docs, not substring), which
# silently matched nothing against these bare IDs and was a real bug
# caught during testing (0-byte output despite a correct ID list).
# subset R1 to the kept ids (same id list keeps pairs in sync)
seqkit grep -f "${OUTDIR}/${SAMPLE}.keep_read_ids.txt" "${R1}" -o "${OUTDIR}/${SAMPLE}_1.kraken_filtered.fastq.gz"
# subset R2 to the same kept ids
seqkit grep -f "${OUTDIR}/${SAMPLE}.keep_read_ids.txt" "${R2}" -o "${OUTDIR}/${SAMPLE}_2.kraken_filtered.fastq.gz"
