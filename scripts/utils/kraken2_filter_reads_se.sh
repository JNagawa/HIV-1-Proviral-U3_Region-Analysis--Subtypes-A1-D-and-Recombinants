#!/bin/bash
# Single-end counterpart of kraken2_filter_reads.sh, for long-read (PacBio
# HiFi / Nanopore) data where each read is one record with no mate. Same
# logic: classify with Kraken2 against a human+bacterial+viral database, then
# drop any read whose assigned taxID descends from Homo sapiens (9606) or
# Bacteria (2), keeping unclassified and viral reads (real HIV-1 signal). The
# only differences from the paired script are `--paired` -> single input and a
# single seqkit grep at the end. See kraken2_filter_reads.sh for the detailed
# rationale on the two different field separators and the seqkit -n pitfall.
#
# Usage: kraken2_filter_reads_se.sh <reads.fastq.gz> <OUTDIR> <SAMPLE_NAME> <KRAKEN2_DB_DIR>
# Produces:
#   <OUTDIR>/<SAMPLE_NAME>.kraken_filtered.fastq.gz
#   <OUTDIR>/<SAMPLE_NAME>.kreport
set -uo pipefail                                     # -u errors on unset vars, pipefail fails a pipe if any stage fails
READS="$1" OUTDIR="$2" SAMPLE="$3" KRAKEN2_DB="$4"   # positional args: input reads, output dir, sample name, Kraken2 DB dir

NODES_DMP="${KRAKEN2_DB}/nodes.dmp"                   # taxonomy tree bundled with the DB, used for the ancestor walk
[ -s "${NODES_DMP}" ] || { echo "ERROR: ${NODES_DMP} not found -- is KRAKEN2_DB (${KRAKEN2_DB}) an extracted Kraken2 database?" >&2; exit 1; }  # bail if the DB isn't a real extracted Kraken2 database

mkdir -p "${OUTDIR}"                                  # ensure the output dir exists

kraken2 --db "${KRAKEN2_DB}" --gzip-compressed --threads "${THREADS:-4}" \
    --output "${OUTDIR}/${SAMPLE}.kraken" \
    --report "${OUTDIR}/${SAMPLE}.kreport" \
    "${READS}"                                        # classify each read against the DB; per-read output + hierarchical report

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
        if (keep) print read_id                       # emit ids of reads to retain
    }
' "${NODES_DMP}" "${OUTDIR}/${SAMPLE}.kraken" > "${OUTDIR}/${SAMPLE}.keep_read_ids.txt"  # write the keep-list of read ids

N_TOTAL=$(wc -l < "${OUTDIR}/${SAMPLE}.kraken")              # total classified reads
N_KEEP=$(wc -l < "${OUTDIR}/${SAMPLE}.keep_read_ids.txt")   # reads kept after host/bacterial removal
echo "Kraken2 filtering for ${SAMPLE}: ${N_KEEP}/${N_TOTAL} reads retained (human/bacterial reads discarded)." >&2  # progress summary to stderr

seqkit grep -f "${OUTDIR}/${SAMPLE}.keep_read_ids.txt" "${READS}" -o "${OUTDIR}/${SAMPLE}.kraken_filtered.fastq.gz"  # subset the FASTQ to just the kept read ids
