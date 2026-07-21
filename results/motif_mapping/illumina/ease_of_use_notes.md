# Ease-of-use notes: motif_mapping (Illumina)

Filled in 2026-07-21 after re-running `compare_motif_mapping_illumina.sh`
following two harness bugs (not tool bugs): the shared U3-extraction
utility wrote lowercase sequences that all three scanners silently
treated as masked, and `run_fimo.sh`/`run_moods.sh`/`run_tfbstools.sh`
were missing their executable bit, so every prior run failed instantly
with an empty log and got miscounted as "0 hits" instead of "didn't run."
With both fixed, all three tools produced real results.

## FIMO

- **Setup friction** (1-5, 1=trivial): 1 -- available via `meme` on
  bioconda, no extra config beyond the JASPAR `.meme` motif file already
  set up by `setup_jaspar.sh`.
- **Documentation quality**: Good -- MEME Suite docs are thorough and the
  TSV output format is documented.
- **Failure modes encountered**: None once the harness bugs above were
  fixed. Fastest of the three (0.14s positive control, 0.17s subset).
- **Output usability**: Clean TSV with p-value/q-value columns, easy to
  filter/parse with awk.
- **Would you use this again for the full cohort? Why/why not**: Yes --
  fastest, cleanest output, found real hits on the positive control (4
  hits) and subset (16 hits).

## MOODS

- **Setup friction** (1-5, 1=trivial): 1 -- bioconda `moods` package,
  works out of the box.
- **Documentation quality**: Sparser than MEME's, but the CLI is simple
  enough not to need much.
- **Failure modes encountered**: None once harness bugs were fixed. Runs
  under a second even on the full subset (0.99s), slower relative
  startup on the tiny positive-control input (0.73s) -- likely fixed
  per-invocation overhead, not scaling with input size.
- **Output usability**: Plain-text hit list, one row per match --
  workable but less self-describing than FIMO's TSV (no q-values).
- **Would you use this again for the full cohort? Why/why not**: Yes,
  as a fast cross-check against FIMO -- found fewer hits (2 positive
  control, 10 subset) which is worth understanding (stricter default
  threshold?) before treating the two tools' counts as directly
  comparable.

## TFBSTools

- **Setup friction** (1-5, 1=trivial): 2 -- R/Bioconductor package,
  heavier dependency footprint than the other two, but installed
  cleanly via the conda env.
- **Documentation quality**: Good Bioconductor vignette, though the R
  API (PWM objects, `searchSeq`) has a steeper learning curve than a
  CLI tool.
- **Failure modes encountered**: None once harness bugs were fixed. By
  far the slowest of the three (20s positive control, 16s subset vs.
  sub-second for FIMO/MOODS) -- R startup + PWM scanning overhead.
- **Output usability**: GFF3 output, most hits of the three tools (14
  positive control, 70 subset) -- likely a more permissive default
  threshold; needs a like-for-like threshold comparison before reading
  too much into the raw hit-count difference vs FIMO/MOODS.
- **Would you use this again for the full cohort? Why/why not**: Yes,
  but budget for its slower runtime, and align its scoring threshold
  with FIMO/MOODS before comparing hit counts directly.

## Overall recommendation for this stage

All three tools are viable and now produce real, positive-control-
validated hits (FIMO 4/16, MOODS 2/10, TFBSTools 14/70 on
positive-control/subset respectively). **FIMO is the current pick** for
the full-cohort run: fastest, cleanest standard output format
(p-value/q-value), and found a sensible middle number of hits. Before
finalizing, the three tools' hit-count spread should be reconciled by
matching significance thresholds (`--thresh` for FIMO, `-p` for MOODS,
`pvalueCutoff` for TFBSTools) rather than comparing each tool's
differing defaults.
