# Tool comparison harness (Objective 1 methodology)

Separate from the production pipelines (`scripts/pipelines/illumina_u3analysis.sh`,
`scripts/pipelines/oxnano_u3analysis.sh`, and the newer per-step
`scripts/pipelines/illumina/`). This harness runs every automatable
candidate tool from `Tools_review_HIV_U3_analysis.pdf` for a given pipeline
step, side by side, on a small sample subset (`subset_samples.tsv`), so tool
choices for the thesis methodology are evidence-backed rather than taken on
the review's recommendation alone.

Web-only tools with no CLI/API (COMET, REGA v3, QGRS Mapper) are
deliberately **excluded** here -- run those yourself and fold the results
into the relevant step's `summary.tsv` / `ease_of_use_notes.md` by hand.

## Layout

Each step gets its own top-level `scripts/<step>_illumina/` directory (a
sibling of this `scripts/common/` directory, not nested under any wrapper
folder) with the same shape:
- `run_<tool>.sh` -- one script per candidate tool, wraps a single
  tool-on-one-sample invocation and does that step's validity check
  (`download_qc_illumina` and `assembly_illumina` instead inline every
  tool as a function in their one combined script -- no separate
  `run_<tool>.sh` files).
- one combined `<step>_illumina.sh` or `compare_<step>_illumina.sh`
  script -- orchestrator: loops tool x sample over the subset, calls
  `measure_and_run` from `common/lib_compare.sh`, appends to
  `results/<step>_illumina/summary.tsv`.

The `_illumina` suffix is there because this whole harness only covers the
Illumina arm so far -- a Nanopore equivalent, if built, would live in
sibling `scripts/<step>_oxnano/` folders (matching the production
pipeline's own `data/raw/oxnano`, `oxnano_u3analysis.sh` naming).

Each step's actual output lives in the mirrored `results/<step>_illumina/`
directory at the repo root (kept separate so `scripts/` stays code-only):
per-run logs/timing files, `summary.tsv`, and `ease_of_use_notes.md`
(copied in from `common/ease_of_use_template.md` and filled in by hand
after running).

Third-party tool clones (SHIVER, Poplars, HIVSeqinR, HIVIntact, jpHMM) live
under `scripts/tools/`, not inside the step folder that uses them -- a
shared location for anything that's someone else's cloned/compiled code
rather than something written for this repo. Downloaded reference data (the
JASPAR PWM set `motif_mapping_illumina` scans against) lives under
`data/reference/jaspar/` for the same reason: it's data, not code.

## Running a step

```bash
conda activate HIV_U3analysis
cd scripts/<step>_illumina
./<step>_illumina.sh          # or ./compare_<step>_illumina.sh, see table below
# then fill in results/<step>_illumina/ease_of_use_notes.md by hand
```

## Steps

| Dir | Script | Compares | Notes |
|---|---|---|---|
| `download_qc_illumina` | `download_qc_illumina.sh` | fastp vs Trimmomatic | review recommends fastp; production currently uses Trimmomatic. Also runs FastQC pre/post as part of the same script. |
| `assembly_illumina` | `compare_assembly_illumina.sh` | BWA+bcftools-consensus vs SPAdes vs SHIVER | see `assembly_illumina/shiver_setup/SOURCE.md` for the one manual step SHIVER needs |
| `msa_illumina` | `compare_msa_illumina.sh` | MAFFT vs MUSCLE vs Clustal Omega | |
| `biological_filtering_illumina` | `compare_biological_filtering_illumina.sh` | Poplars vs HIVSeqinR vs HIVIntact | all three are GitHub-hosted, cloned locally by `setup_*.sh` into `scripts/tools/` |
| `subtyping_illumina` | `compare_subtyping_illumina.sh` | jpHMM vs IQ-TREE2 | COMET/REGA excluded (web-only) |
| `motif_mapping_illumina` | `compare_motif_mapping_illumina.sh` | FIMO vs TFBSTools vs MOODS | U3 extraction (via `scripts/utils/extract_u3_by_hxb2_anchor.sh`) then positive control against HXB2's own U3, then the subset |
| (none yet) | | gquad vs pqsfinder | G-quadruplex prediction -- not yet built |

Run in this order, since each step after `download_qc_illumina` consumes
the previous one's `results/<step>_illumina/` output:
`download_qc_illumina` -> `assembly_illumina` -> `msa_illumina` ->
`biological_filtering_illumina` / `subtyping_illumina` (either order, both
only need `assembly_illumina`'s output) -> `motif_mapping_illumina` (needs
`msa_illumina`'s alignment).

No Nanopore-specific trim/assembly comparison step exists: the tools-review
document offers no competing alternative to Porechop_ABI/NanoFilt or to
minimap2 for Nanopore data -- its assembly recommendations are Illumina
(SHIVER/SPAdes/BWA) and PacBio (SMRT Link ccs) only. The Nanopore arm only
needed the production-script correctness fixes, not a comparison.

Once you've reviewed each step's `summary.tsv` + notes and picked a winner,
that decision gets wired into the production pipelines (steps 8/10/11 of
the still-current monolithic scripts, currently stubs pending exactly this).
