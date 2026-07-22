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
set -uo pipefail
READS="$1" OUTDIR="$2" SAMPLE="$3" KRAKEN2_DB="$4"

NODES_DMP="${KRAKEN2_DB}/nodes.dmp"
[ -s "${NODES_DMP}" ] || { echo "ERROR: ${NODES_DMP} not found -- is KRAKEN2_DB (${KRAKEN2_DB}) an extracted Kraken2 database?" >&2; exit 1; }

mkdir -p "${OUTDIR}"

kraken2 --db "${KRAKEN2_DB}" --gzip-compressed --threads "${THREADS:-4}" \
    --output "${OUTDIR}/${SAMPLE}.kraken" \
    --report "${OUTDIR}/${SAMPLE}.kreport" \
    "${READS}"

awk '
    NR==FNR {
        split($0, f, "\t\\|\t")
        parent[f[1]] = f[2]
        next
    }
    {
        split($0, f, "\t")
        read_id = f[2]
        t = f[3]
        keep = 1
        depth = 0
        while (t != "1" && t != "" && depth < 50) {
            if (t == "9606" || t == "2") { keep = 0; break }
            if (!(t in parent)) break
            t = parent[t]
            depth++
        }
        if (keep) print read_id
    }
' "${NODES_DMP}" "${OUTDIR}/${SAMPLE}.kraken" > "${OUTDIR}/${SAMPLE}.keep_read_ids.txt"

N_TOTAL=$(wc -l < "${OUTDIR}/${SAMPLE}.kraken")
N_KEEP=$(wc -l < "${OUTDIR}/${SAMPLE}.keep_read_ids.txt")
echo "Kraken2 filtering for ${SAMPLE}: ${N_KEEP}/${N_TOTAL} reads retained (human/bacterial reads discarded)." >&2

seqkit grep -f "${OUTDIR}/${SAMPLE}.keep_read_ids.txt" "${READS}" -o "${OUTDIR}/${SAMPLE}.kraken_filtered.fastq.gz"
