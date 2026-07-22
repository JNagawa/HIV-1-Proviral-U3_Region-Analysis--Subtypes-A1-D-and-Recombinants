# Kraken2 size-capped Standard database

Used by `scripts/utils/kraken2_filter_reads.sh` for host (human) and
bacterial contamination filtering, per the tools review's recommended
`fastp -> Kraken2 -> SHIVER` chain.

A viral-only database was considered first (much smaller, ~0.5GB) but
rejected: it contains no human or bacterial genomes, so it cannot detect
host contamination at all -- the dominant contamination source in HIV
proviral sequencing (human genomic DNA co-purified with the provirus). The
size-capped Standard database below still contains bacteria, archaea,
viral, and human genomes (so it can actually do the job), just at a
smaller memory footprint (~16GB in-memory index) than the full uncapped
Standard database (50GB+, the RAM footprint the tools review itself flags
as Kraken2's main limitation).

- **Source**: https://benlangmead.github.io/aws-indexes/k2 (Ben Langmead's
  pre-built Kraken2/Bracken index collection)
- **Download URL**: `https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16_GB_20260626.tar.gz`
- **Downloaded**: 2026-07-21
- **Archive size**: 11.2 GB compressed
- **Contents**: RefSeq archaea, bacteria, viral, plasmid, human, and
  UniVec_Core sequences, minimizer-subsampled to cap the in-memory index
  at ~16GB
- **Build date**: 2026-06-26

Extracted here (not gitignored the way `*.fasta` reference files are,
since this is a directory of Kraken2's own binary index files, not a
FASTA -- see `.gitignore` if this needs excluding too; at ~16GB it likely
shouldn't be committed to git regardless, re-download instead using the
URL above).
