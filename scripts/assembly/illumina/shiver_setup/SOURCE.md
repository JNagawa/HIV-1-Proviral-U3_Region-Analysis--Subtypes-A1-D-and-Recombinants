# SHIVER reference alignment -- manual acquisition step

**STATUS (2026-07-22): resolved.** The alignment has been downloaded and
placed at `HIV1_COM_ref_alignment.fasta` (180 sequences, HXB2 + Group M
subtypes/CRFs). Correction to the earlier note below: this turned out to
be scriptable after all -- see "How this was actually obtained" for the
reusable curl recipe (e.g. for pulling a newer year's compendium later).

`shiver_init.sh` requires a curated alignment of existing HIV-1 reference
sequences (its 3rd positional argument) covering the diversity you might
find in your samples. The original plan below assumed LANL's alignment
tool required a human in a browser -- it doesn't; its dropdowns are
populated via a plain AJAX endpoint (`list.comp`) and the form itself
posts to a plain CGI script (`align.cgi`) that returns a `download.cgi`
link, all scriptable with curl (no login/session/CAPTCHA involved):

1. Go to the LANL HIV Sequence Database alignment tool:
   `https://www.hiv.lanl.gov/content/sequence/HIV/mainpage.html` -> the
   "Genome Alignments" / compendium alignment page (SHIVER's own README
   points at "the 2021 compendium genome alignment for HIV-1 group M
   including recombinants").
2. Select: organism = HIV-1, region = complete genome, subtype = all Group M
   subtypes + recombinants (CRF/URF), one sequence per patient, format =
   FASTA.
3. Download the result and place it here as
   `scripts/assembly/illumina/shiver_setup/HIV1_COM_ref_alignment.fasta`
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

## How this was actually obtained

The form at `https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html`
posts (as `multipart/form-data`) to `/cgi-bin/NEWALIGN/align.cgi`. Its
cascading dropdowns are populated by `list.comp?col=<field>&where=<prior
selections>` -- queried directly below to confirm valid values instead of
guessing:

```bash
BASE="https://www.hiv.lanl.gov/content/sequence/NEWALIGN"
curl -s "$BASE/list.comp?submit=Retrieve&server=HIV&ebolaIF=&col=al_align_type"
# -> COM = Compendium
curl -s "$BASE/list.comp?submit=Retrieve&server=HIV&ebolaIF=&col=al_region&where=al_align_type:COM,al_organism:HIV1"
# -> GENOME = complete genome
# GENO_SUB=ALLM ("M group with CRFs") comes from a client-side JS map in
# align.html (fill_GENO_SUB), not list.comp -- not an AJAX-backed field.

curl -s -X POST "https://www.hiv.lanl.gov/cgi-bin/NEWALIGN/align.cgi" \
  -F "ORGANISM=HIV" -F "ALIGN_TYPE=COM" -F "SUBORGANISM=HIV1" \
  -F "PRE_USER=predefined" -F "REGION=GENOME" -F "START=" -F "END=" \
  -F "GENO_SUB=ALLM" -F "BASETYPE=DNA" -F "YEAR=2021" -F "alignmentID=" \
  -F "down_acc=1" -F "FORMAT=fasta" -F "submit=Get Alignment" \
  -o align_response.html
# response HTML contains a link:
#   /cgi-bin/common_code/download.cgi?/tmp/NEWALIGN/<session-id>/HIV1_COM_2021_genome_DNA.fasta
# (session-id is per-request and short-lived -- fetch immediately)
curl -s "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/tmp/NEWALIGN/<session-id>/HIV1_COM_2021_genome_DNA.fasta" \
  -o HIV1_COM_ref_alignment.fasta
```

## Download log

- Date downloaded: 2026-07-22
- Exact URL / tool page used: `https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html` (via the curl recipe above, not the browser form)
- Selection filters applied: ALIGN_TYPE=Compendium, SUBORGANISM=HIV-1/SIVcpz, REGION=GENOME (complete genome), GENO_SUB=ALLM (M group with CRFs), BASETYPE=DNA, YEAR=2021, FORMAT=fasta
- Result: 180 sequences (HXB2 K03455 + 179 Group M subtype/CRF representatives)
