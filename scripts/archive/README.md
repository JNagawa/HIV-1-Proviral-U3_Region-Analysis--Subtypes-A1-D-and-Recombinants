# Archived scripts

These two scripts are earlier, superseded versions of the production pipelines
(`scripts/pipelines/illumina_u3analysis.sh`, `scripts/pipelines/oxnano_u3analysis.sh`). They only
cover steps 1-6 (SRA download through pre/post-trim QC + MultiQC) and predate
steps 7-11 (reference mapping, biological filtering, MSA, subtyping, U3
extraction / motif mapping). Kept here for reference rather than deleted,
since `qc_oxnano_u3analysis.sh` documents a genuinely different, stricter
Nanopore filtering choice worth citing in the methods write-up.

## Open discrepancy — not resolved, needs a decision

`qc_oxnano_u3analysis.sh` filters at `MIN_QUALITY=10` / `MIN_LENGTH=7000`
(strict — keeps only near-full-length reads), while the current
`oxnano_u3analysis.sh` uses `MIN_QUALITY=7` / `MIN_LENGTH=200` (permissive —
keeps partial reads too). This directly affects what fraction of each
Nanopore sample's reads reach assembly. Given the thesis's explicit framing
around *near-full-length* proviral genomes, the stricter legacy threshold may
be the more defensible choice methodologically — but this changes real
results, so it's a decision for Jovita to make deliberately, not something
resolved silently during the pipeline restructure.
