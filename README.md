# HIV-1 Proviral U3 Region Analysis

[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20HPC%20(Slurm)-blue)]()
[![Conda](https://img.shields.io/badge/Env-Conda-green)]()

Recovers the U3 region of HIV-1 proviral genomes and maps host transcription-factor binding
sites (TFBS) and G-quadruplex structures within it. U3 is the HIV-1 promoter, so its TFBS
content bears on transcription and latency.

The repository is organised as a **tool-comparison harness**: at every stage several
candidate tools are run over the same input and each writes a row to a shared
`summary.tsv`, so tool choices are made on measured evidence rather than convention.

**Platforms:** Illumina (short read) and PacBio HIV-SMRTcap (HiFi long read) are both
implemented. Oxford Nanopore directories exist as scaffolding but are not populated.

---

## Repository structure

```text
Internship/
├── HIV_U3analysis_env.yml         # conda environment
├── README.md
├── data/
│   ├── raw/{illumina,pacbio}/     # downloaded reads
│   ├── reference/                 # HXB2 (K03455.1) FASTA + cached GenBank record + JASPAR PWMs
│   └── processed/
├── scripts/
│   ├── <stage>/<platform>/        # one folder per stage per platform (see below)
│   ├── common/                    # lib_compare.sh (measure_and_run, append_summary_row) + README
│   ├── tools/                     # third-party clones: shiver, jpHMM, Poplars, HIVSeqinR, HIVIntact
│   ├── utils/                     # U3 extraction, report generators, deck builders, chain launcher
│   ├── pipelines/                 # older monolithic pipelines (superseded by the per-stage harness)
│   └── archive/                   # superseded scripts
├── results/<stage>/<platform>/    # summary.tsv, tool outputs, timing logs, ease_of_use_notes.md
│   └── <stage>/pacbio/<arm>/      # PacBio is additionally split by assembly arm (see "Two arms")
├── writeups/                      # reports, slide decks, proposal, tools review
├── logs/                          # Slurm stdout/err
└── checkpoints/                   # resume markers for the older pipelines
```

Stages, in pipeline order: `download_qc`, `proviral_extraction` (PacBio only), `assembly`,
`msa`, `biological_filtering`, `subtyping`, `motif_mapping`.

---

## Pipeline stages, tools compared, and tools chosen

| stage | tools compared | Illumina choice | PacBio choice |
|---|---|---|---|
| QC / host depletion | NanoFilt, chopper, fastp, seqkit, awk, Kraken2, 4 dedup strategies | Trimmomatic/fastp + Kraken2 | seqkit or awk + Kraken2 (FASTQ-only tools are `not_applicable` on FASTA input) |
| Proviral extraction | minimap2 mask + N-strip | n/a | minimap2mask+strip |
| Assembly | SHIVER, SPAdes, bwa+bcftools consensus, minimap2+bcftools, hifiasm | `bwa_consensus` (SPAdes yields only 420–1786 bp) | reference-guided minimap2+bcftools, **two arms**; hifiasm fails structurally |
| MSA | MAFFT, MUSCLE, Clustal Omega | MAFFT | **MAFFT** — Clustal Omega mis-handles N-masked input |
| U3 extraction | HXB2-anchored alignment liftover | shared | shared |
| Subtyping | jpHMM, IQ-TREE 2 | both | both |
| Motif mapping | FIMO, TFBSTools, MOODS; gquad, pqsfinder | FIMO (+ others as cross-checks) | FIMO (+ others as cross-checks) |
| Biological filtering | HIV-Intact, Poplars, HIVSeqinR | HIV-Intact + Poplars | **HIV-Intact + Poplars**; HIVSeqinR unsuitable |

### Why HIVSeqinR is not used

HIVSeqinR requires *"linear HIV genomes WITH FLANKING PRIMER binding sites at 5′ and 3′
ends"* from a specific de novo pipeline, and rejects `N` and IUPAC codes outright. Its
provisioning is now fully automated (`setup_tools.sh` unpacks the release, installs the R
`muscle` and `pwalign` packages, patches nine defunct Biostrings calls, builds the HXB2
BLAST database and sets `MyBlastnDir`), but its 5′Psi+gag filter still rejects amplicon
segments. The 2nd-round PCR primers for the SMRTcap libraries are also undocumented, so
primer autotrim cannot work. Treat it as inapplicable, not merely unconfigured.

---

## Two arms (PacBio)

The PacBio assembly step runs twice over the same reads, so the cost of the reference
choice stays measurable all the way to the motif hits:

| arm | reference | role |
|---|---|---|
| `minimap2_consensus` | HXB2 (K03455.1) | fixed subtype-B baseline; the coordinate system everything downstream expects |
| `minimap2_bestref` | closest LTR-complete Group M panel entry | measures the cost of forcing subtype A1/D reads onto a subtype B backbone |

Every downstream stage takes an `ASSEMBLY_ARM` environment variable and writes to
`results/<stage>/pacbio/<arm>/`, so the two passes never collide.

---

## Setup

```bash
conda env create -f HIV_U3analysis_env.yml
conda activate HIV_U3analysis

# third-party tools (clones, venvs, R packages, BLAST databases)
bash scripts/biological_filtering/illumina/setup_tools.sh
bash scripts/subtyping/illumina/setup_jphmm.sh
```

`setup_tools.sh` is shared by both platforms and is idempotent — safe to re-run.

---

## Running

**Never run compute on the login node.** Every stage is submitted through the generic
wrapper, which activates the conda environment and exports `THREADS`:

```bash
# one stage, Illumina
sbatch scripts/utils/run_comparison_step.slurm.sh \
       scripts/assembly/illumina/compare_assembly_illumina.sh

# one stage, PacBio, for a specific arm
sbatch --export=ALL,ASSEMBLY_ARM=minimap2_consensus \
       scripts/utils/run_comparison_step.slurm.sh \
       scripts/assembly/pacbio/compare_assembly_pacbio.sh

# the whole PacBio chain, both arms, wired with afterok dependencies
bash scripts/utils/run_pacbio_chain.sh
```

Override which tools a stage runs with the stage's `*_TOOLS` variable, e.g.
`ASSEMBLY_TOOLS="minimap2_bestref"`.

### Generated reports

| output | produced by |
|---|---|
| `results/<stage>/<platform>/summary.tsv` | every stage, one row per (tool, sample) |
| `results/assembly/pacbio/assembly_report.md` | auto after PacBio assembly |
| `results/biological_filtering/pacbio/<arm>/intactness_basis_report.md` | auto after PacBio filtering |
| `writeups/pacbio_pipeline_report.md` | `scripts/utils/report_*.sh` + hand-written synthesis |
| `writeups/PacBio_Pipeline_Report.pptx` | `scripts/utils/build_pacbio_pipeline_pptx.sh` (login node — needs `zip`) |

---

## Interpreting results — read this first

Two distinctions the pipeline goes to some length to preserve:

- **unsequenced** — no read spans the position. Masked to `N`. Nothing can be concluded
  about it in either direction. A motif absent over an unsequenced stretch is **not
  evidence of absence**.
- **deleted** — reads span the position and agree the base is absent. A called variant,
  not an `N`. This is a positive finding.

Consequences worth knowing:

- **Validity is called (non-N) bases, not length.** A 9719 bp consensus that is 88% `N` is
  not a recovered genome.
- **The PacBio libraries are amplicon, not shotgun.** Reads stack into a few discrete
  intervals rather than tiling the genome, so de novo assembly is structurally impossible
  for most samples and only some genes can be assessed for intactness.
- **U3 is taken from whichever LTR a sample actually sequenced.** Both LTRs are identical in
  an integrated provirus, but an amplicon may reach only one.
- **U3 coordinates are HXB2-anchored via alignment columns, never by raw genome offsets.**
  U3 position varies between subtypes partly through insertions, so coordinate slicing is
  incorrect and no per-subtype coordinate system exists to substitute.

---

## Current state and known caveats

### PacBio — corrected (2026-08-05)

Six of seven stages complete for both arms; subtyping was still running at last update.
Defects found and fixed: indel calling was silently off (`bcftools mpileup -X pacbio-ccs`
now set); the `N` mask conflated "deleted" with "not sequenced" (coverage now derived from
read alignment spans); U3 was taken from the wrong LTR; HIV-Intact scored every sequence as
subtype A1 including HXB2; Poplars silently required equal-length input. See
`writeups/pacbio_pipeline_report.md`.

### Illumina — carries the same defects, not yet fixed

The Illumina arm predates those corrections and shows the same signatures. **Do not treat
its current outputs as final:**

| symptom | evidence | root cause (fixed on the PacBio side) |
|---|---|---|
| Consensus reported as `9719bp, 0.00% N` for every sample | `results/assembly/illumina/summary.tsv` | uncovered positions filled with reference instead of `N` |
| MSA returns exactly 9719 columns, 0.0% gap | `results/msa/illumina/summary.tsv` | zero indels being called, so all consensuses are reference-length |
| `0 intact, 5 non-intact`, with **HXB2 itself** reported defective | `results/biological_filtering/illumina/hivintact_out/errors.json` | `--subtype A1` applied to the whole batch including the subtype-B reference |
| HIVSeqinR exits 1 | `summary.tsv` | provisioning, since fixed in the shared `setup_tools.sh` |

The shared wrappers (`run_hivseqinr.sh`, `setup_tools.sh`) already carry their fixes, so
they apply to Illumina automatically. The remaining work is Illumina-specific: N-masking and
indel calling in `compare_assembly_illumina.sh`, and per-sample subtypes in
`compare_biological_filtering_illumina.sh`.

---

## Conventions

- **Bash and awk only.** No Python analysis scripts, including for report and `.pptx`
  generation (`.pptx` is emitted as OOXML directly, since pandoc, LibreOffice and
  python-pptx are all absent here).
- **Slurm for all compute.** Nothing heavier than a file listing runs on the login node.
- **`peak_rss_mb` is always blank** — GNU `time` is unavailable on this cluster, so runtime
  is measured and memory is not.
- **Long comments sit above the line they describe**, not trailing it.
- Compute nodes lack `zip`/`unzip`; `bsdtar` is used instead where archives are handled.

## License

Developed for academic research as part of an MSc Bioinformatics internship at Makerere
University.
