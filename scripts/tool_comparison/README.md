# Tool comparison harness (Objective 1 methodology)

Separate from the production pipelines (`scripts/pipelines/illumina_u3analysis.sh`,
`scripts/pipelines/oxnano_u3analysis.sh`). This workspace runs every
automatable candidate tool from `Tools_review_HIV_U3_analysis.pdf` for a
given pipeline stage, side by side, on a small sample subset
(`subset_samples.tsv`), so tool choices for the thesis methodology are
evidence-backed rather than taken on the review's recommendation alone.

Web-only tools with no CLI/API (COMET, REGA v3, QGRS Mapper) are
deliberately **excluded** here -- run those yourself and fold the results
into the relevant stage's `summary.tsv` / `ease_of_use_notes.md` by hand.

## Layout

Scripts live here under `scripts/tool_comparison/`; each stage's actual
output lives in the mirrored `results/tool_comparison/0N_<stage>/` directory
at the repo root (kept separate so `scripts/` stays code-only). Each
`0N_<stage>/` script directory has the same shape:
- `run_<tool>.sh` -- one script per candidate tool, wraps a single
  tool-on-one-sample invocation and does that stage's validity check.
- `compare_<stage>.sh` -- orchestrator: loops tool x sample over the subset,
  calls `measure_and_run` from `common/lib_compare.sh`, appends to
  `results/tool_comparison/0N_<stage>/summary.tsv`.

Each stage's `results/tool_comparison/0N_<stage>/` directory holds per-run
logs/timing files, `summary.tsv`, and `ease_of_use_notes.md` (copied in from
`common/ease_of_use_template.md` and filled in by hand after running).

## Running a stage

```bash
conda activate HIV_U3analysis
cd scripts/tool_comparison/0N_<stage>
./compare_<stage>.sh
# then fill in results/tool_comparison/0N_<stage>/ease_of_use_notes.md by hand
```

## Stages

| Dir | Compares | Notes |
|---|---|---|
| `01_qc_trim_illumina` | fastp vs Trimmomatic | review recommends fastp; production currently uses Trimmomatic |
| `02_assembly_illumina` | BWA+bcftools-consensus vs SPAdes vs SHIVER | see `02_assembly_illumina/shiver_setup/SOURCE.md` for the one manual step SHIVER needs |
| `03_msa` | MAFFT vs MUSCLE vs Clustal Omega | |
| `04_biological_filtering` | Poplars vs HIVSeqinR vs HIVIntact | all three are GitHub-hosted, cloned locally by `setup_*.sh` |
| `05_subtyping` | jpHMM vs IQ-TREE2 | COMET/REGA excluded (web-only) |
| `06_motif_mapping` | FIMO vs TFBSTools vs MOODS | positive control against HXB2's own U3 runs first |
| `07_gquadruplex` | gquad vs pqsfinder | |

No Nanopore-specific trim/assembly comparison stage exists: the tools-review
document offers no competing alternative to Porechop_ABI/NanoFilt or to
minimap2 for Nanopore data -- its assembly recommendations are Illumina
(SHIVER/SPAdes/BWA) and PacBio (SMRT Link ccs) only. The Nanopore arm only
needed the production-script correctness fixes, not a comparison.

Once you've reviewed each stage's `summary.tsv` + notes and picked a winner,
that decision gets wired into steps 8/10/11 of the production scripts
(currently stubs pending exactly this).
