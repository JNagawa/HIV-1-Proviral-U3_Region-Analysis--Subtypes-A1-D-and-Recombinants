# SHIVER reference alignment -- manual acquisition step

`shiver_init.sh` requires a curated alignment of existing HIV-1 reference
sequences (its 3rd positional argument) covering the diversity you might
find in your samples. This is **not scriptable** from this environment --
LANL's alignment tool is form-driven with no static download URL -- so this
is a one-time manual step:

1. Go to the LANL HIV Sequence Database alignment tool:
   `https://www.hiv.lanl.gov/content/sequence/HIV/mainpage.html` -> the
   "Genome Alignments" / compendium alignment page (SHIVER's own README
   points at "the 2021 compendium genome alignment for HIV-1 group M
   including recombinants").
2. Select: organism = HIV-1, region = complete genome, subtype = all Group M
   subtypes + recombinants (CRF/URF), one sequence per patient, format =
   FASTA.
3. Download the result and place it here as
   `scripts/tool_comparison/02_assembly_illumina/shiver_setup/HIV1_COM_ref_alignment.fasta`
   (this filename is what `run_shiver.sh` expects; it's gitignored like all
   other `*.fasta` files in this repo, so it stays untracked -- record here
   once downloaded: exact URL used, selection filters, and download date,
   so this is reproducible without redistributing LANL data through git).

This alignment is Group-M-wide with recombinants, so it already spans the
A1/D/CRF-URF diversity relevant to this cohort -- no further curation
needed.

## Adapters / primers

`run_shiver.sh` uses SHIVER's own bundled example files as a starting point:
- Adapters: `scripts/tools/shiver/data/example_input/adapters_Illumina.fasta` (same
  TruSeq family already used by Trimmomatic in the production pipeline).
- Primers: `scripts/tools/shiver/data/example_input/primers_GallEtAl2012.fasta` -- **flagged
  as an approximation**, not confirmed to match PRJNA207834's actual library
  prep (not documented anywhere in this repo). Primer trimming is optional
  in SHIVER; if the real protocol's primers can't be sourced, an empty/
  minimal primers file is an acceptable fallback.

## Download log

(Fill in once you've done step 1-3 above:)
- Date downloaded:
- Exact URL / tool page used:
- Selection filters applied:
