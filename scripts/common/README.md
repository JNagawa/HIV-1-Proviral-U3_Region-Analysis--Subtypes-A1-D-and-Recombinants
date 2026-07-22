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

## PacBio (HIV-SMRTcap) arm

A parallel harness under `scripts/<step>/pacbio/` + `results/<step>/pacbio/`,
same shape as the Illumina arm, for the HIV-SMRTcap long-read (PacBio HiFi)
data. Data source: SRA BioProjects PRJNA1251704 + PRJNA1252688 (Lee et al.
2024, Nat Commun 15:5513). These are NOT pre-downloaded -- unlike the Illumina/
Nanopore raw files -- so the first step fetches them.

Subtype is not in the SRA metadata, so "A1/D/recombinants only" is applied by
cohort membership: the `RAKAIUG_donor*` runs are the Rakai Uganda cohort (the
A1/D/recombinant target population), 92UG_029 (subtype A/A1) and 93UG_065
(subtype D) are known-subtype reference anchors, and the subtype-B/C reference
strains + cell-line controls are excluded. See `subset_samples.tsv` for the
exact accessions and rationale. The actual per-sample A1/D/recombinant call is
what the subtyping step produces.

The SMRTcap data is already HiFi/CCS (the SMRT Link `ccs` step is done before
SRA deposition), so the arm starts from HiFi reads, not raw subreads. The QC
recommendation (NanoPlot) and platform-agnostic downstream tools are taken
straight from the tools review; long-read swaps replace the Illumina-specific
tools. Steps and comparisons:

| Dir | Script | Compares | Notes |
|---|---|---|---|
| (download) | `download_qc/pacbio/download_pacbio.slurm.sh` | -- | `sbatch` job: prefetch + fasterq-dump the subset into `data/raw/pacbio/`. Multi-GB HiFi, so its own job. |
| `download_qc/pacbio` | `download_qc_pacbio.sh` | NanoFilt vs fastp vs chopper (filter); fastp vs seqkit (dedup) | NanoPlot QC + single-end Kraken2 host removal (`utils/kraken2_filter_reads_se.sh`) between them. |
| `proviral_extraction/pacbio` | `extract_provirus_pacbio.sh` | (single method) | **The special step.** N-masks host flanks (minimap2->HXB2, soft-clip = host) then strips them with `utils/extract_provirus_strip_hostN.sh` to recover the ACGT proviral core. `PROVIRUS_INPUT_MASKED=1` skips masking if input is already host-N-masked. |
| `assembly/pacbio` | `compare_assembly_pacbio.sh` | minimap2->HXB2 consensus vs hifiasm | reference-guided vs de novo, mirroring Illumina bwa-vs-SPAdes. hifiasm is the HiFi de novo assembler (review names none). |
| `msa/pacbio` | `compare_msa_pacbio.sh` | MAFFT vs MUSCLE vs Clustal Omega | reuses `scripts/msa/illumina/run_*.sh`. |
| `biological_filtering/pacbio` | `compare_biological_filtering_pacbio.sh` | Poplars vs HIVSeqinR vs HIVIntact | reuses `scripts/biological_filtering/illumina/run_*.sh` + the shared `scripts/tools/` clones. |
| `subtyping/pacbio` | `compare_subtyping_pacbio.sh` | jpHMM vs IQ-TREE2 | reuses `scripts/subtyping/illumina/run_*.sh`. Determines each sample's A1/D/recombinant subtype. |
| `motif_mapping/pacbio` | `compare_motif_mapping_pacbio.sh` | FIMO vs MOODS vs TFBSTools; gquad vs pqsfinder | reuses `scripts/motif_mapping/illumina/run_*.sh` + shared U3 extraction util. |

The four downstream steps (msa, biological_filtering, subtyping, motif_mapping)
are platform-agnostic -- they operate on assembled FASTA -- so the PacBio
orchestrators reuse the Illumina arm's `run_<tool>.sh` scripts unchanged and
only differ in input/output paths. Run order:
`download_pacbio.slurm.sh` -> `download_qc_pacbio.sh` ->
`extract_provirus_pacbio.sh` -> `compare_assembly_pacbio.sh` ->
`compare_msa_pacbio.sh` -> `compare_biological_filtering_pacbio.sh` /
`compare_subtyping_pacbio.sh` -> `compare_motif_mapping_pacbio.sh`. New env
tools (`hifiasm`, `chopper`) are pinned in `HIV_U3analysis_env.yml`; COMET/
REGA v3/QGRS Mapper remain web-only and excluded here, same as the Illumina arm.

Once you've reviewed each step's `summary.tsv` + notes and picked a winner,
that decision gets wired into the production pipelines (steps 8/10/11 of
the still-current monolithic scripts, currently stubs pending exactly this).
