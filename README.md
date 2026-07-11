# HIV-1 Proviral U3 Region — Transcription Factor Binding Site (TFBS) & G-Quadruplex Analysis Pipeline

[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20HPC%20(Slurm)-blue)]()
[![Conda](https://img.shields.io/badge/Env-Conda-green)]()
[![License](https://img.shields.io/badge/License-Academic%20Use-lightgrey)]()

> **Author:** Jovita Nagawa — MSc Bioinformatics, Makerere University  
> **Supervisors:** *(As listed in the research proposal)*  
> **Date:** July 2026  
> **Status:** Active Development — QC & Preprocessing Stage

---

## Table of Contents

- [Project Overview](#project-overview)
- [Research Objectives](#research-objectives)
- [Datasets](#datasets)
- [Analysis Pipeline](#analysis-pipeline)
- [Repository Structure](#repository-structure)
- [Installation & Environment Setup](#installation--environment-setup)
- [Usage](#usage)
- [Pipeline Stages & Tools](#pipeline-stages--tools)
- [Key Design Decisions](#key-design-decisions)
- [References](#references)

---

## Project Overview

The **U3 region** of the HIV-1 Long Terminal Repeat (LTR) promoter is a highly variable regulatory element that controls viral transcription, replication, and latency. Different HIV-1 subtypes — including **A1**, **D**, and **A1-D recombinants** prevalent in East Africa — harbour distinct transcription factor binding site (TFBS) profiles in their U3 regions that may influence viral behaviour and treatment outcomes.

This project implements a reproducible bioinformatics pipeline to:

1. Process raw high-throughput sequencing data from both **short-read (Illumina MiSeq)** and **long-read (Oxford Nanopore GridION)** platforms.
2. Perform rigorous quality control, adapter trimming, and contamination removal.
3. Assemble near-full-length HIV-1 proviral genomes.
4. Filter out defective proviruses (APOBEC-hypermutated or structurally deleted genomes).
5. Align intact genomes, extract the **U3 promoter region**, and systematically map **host transcription factor binding sites** (NF-κB, SP1, NFAT, AP-1, COUP-TF, USF, ETS-1, LEF-1) and **G-quadruplex-forming sequences**.

The pipeline is based on the methodology described in the research proposal (*Jovita_HIV1_U3_TFBS_Research_proposal_Draftv2.pdf*) and the tools review document (*Tools_review_HIV_U3_analysis.pdf*).

---

## Research Objectives

1. **Extract U3 sequences** from intact near-full-length HIV-1 proviral genomes (FLIP-SEQ and HIV SMRTcap) in the Rakai Community Cohort.
2. **Stratify by subtype** — A1, D, and CRF/URF (circulating/unique recombinant forms).
3. **Map TFBS profiles** to identify subtype-specific regulatory motifs that may affect promoter activity.
4. **Predict G-quadruplex structures** within the U3 region that could influence LTR-driven transcription and latency.
5. **Compare TFBS landscapes** across subtypes to identify potential targets for subtype-specific therapeutic or diagnostic strategies.

---

## Datasets

### Short-Read Illumina Dataset

| Property | Value |
|:---|:---|
| **BioProject** | [PRJNA207834](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA207834) |
| **Study** | *Prevalence and Clinical Impacts of HIV-1 Intersubtype Recombinants in Uganda* |
| **Source** | Treatment-naive individuals, rural Mbarara, Uganda |
| **Platform** | Illumina MiSeq — 2×251 bp paired-end |
| **Samples** | 24 paired-end samples (`SRR908430`–`SRR908453`) |
| **Subtypes** | HIV-1 A1, D, and A1-D intersubtype recombinants |
| **Data Volume** | ~3.3 GB (compressed FASTQ) |
| **Local Directory** | `illumina_raw_data/` |

### Long-Read Nanopore Dataset

| Property | Value |
|:---|:---|
| **BioProject** | [PRJNA765218](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA765218) |
| **Study** | *NanoHIV — Oxford Nanopore GridION sequencing of near-full-length HIV-1 genomes* |
| **Source** | Stellenbosch University, South Africa |
| **Platform** | Oxford Nanopore Technologies (ONT) GridION — single-end long reads |
| **Samples** | 9 single-end samples (`SRR16005710`–`SRR16005718`) |
| **Data Volume** | ~13.0 GB (compressed FASTQ) |
| **Local Directory** | `raw_data/` |

> **Why two platforms?** Short reads (Illumina) offer high per-base accuracy for SNP-level TFBS mapping, while long reads (Nanopore) can span the entire ~9 kb proviral genome in a single read, resolving complex insertions, deletions, and LTR duplications that confuse short-read assemblers.

---

## Analysis Pipeline

The full pipeline proceeds through **seven stages**, as outlined in the tools review document:

```
┌─────────────────────────────────────────────────────────────────────────┐
│  STAGE 1: Quality Control & Preprocessing                              │
│  ┌──────────────────────┐    ┌──────────────────────┐                  │
│  │ Illumina (Short-Read)│    │ Nanopore (Long-Read) │                  │
│  │ FastQC → Trimmomatic │    │ NanoPlot/NanoQC/      │                  │
│  │ → FastQC → MultiQC   │    │ NanoStat → Porechop  │                  │
│  │                      │    │ → NanoFilt → MultiQC  │                  │
│  └──────────┬───────────┘    └──────────┬───────────┘                  │
│             └──────────┬───────────────┘                               │
│                        ▼                                               │
│  STAGE 2: Genome Assembly                                              │
│  SHIVER (reference-guided iterative assembly)                          │
│                        ▼                                               │
│  STAGE 3: Biological Filtering & Intactness                            │
│  Poplars/Hypermut 3 (APOBEC hypermutation) + HIVSeqinR (structural)   │
│                        ▼                                               │
│  STAGE 4: Multiple Sequence Alignment                                  │
│  MAFFT (anchored to HXB2 reference coordinates)                       │
│                        ▼                                               │
│  STAGE 5: Subtype and Recombination Detection                          │
│  jpHMM (breakpoints) + IQ-TREE 2 (phylogeny) + COMET/REGA v3          │
│                        ▼                                               │
│  STAGE 6: U3 Region Extraction                                         │
│  SeqKit (coordinate-based sub-sequence extraction)                     │
│                        ▼                                               │
│  STAGE 7: TFBS Motif Scanning                                          │
│  FIMO (from MEME Suite) & TFBSTools — scanning for NF-κB, SP1, etc.   │
│                        ▼                                               │
│  STAGE 8: G-Quadruplex Prediction                                      │
│  R gquad & pqsfinder — predicting G4-forming sequences in U3          │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Repository Structure

```
/home/jnagawa/Internship/
│
├── README.md                          # This file
├── HIV_U3analysis_env.yml             # Conda environment specification
├── .gitignore                         # Git tracking rules (ignores large data files)
├── progress_report.md                 # Detailed project progress report
├── progress_report.docx               # Progress report (Word format)
│
├── qc_illumina_u3analysis.sh         # Illumina preprocessing pipeline (Slurm/Bash)
│                                      #   → SRA download, FastQC, Trimmomatic, MultiQC
├── qc_u3analyis.sh                    # Nanopore preprocessing pipeline (Slurm/Bash)
│                                      #   → SRA download, NanoPlot, Porechop, NanoFilt, MultiQC
├── variant_call.sh                    # Variant calling pipeline (BWA, SAMtools, BCFtools, SnpEff)
├── test_tools.sh                      # Diagnostic script to verify tool availability
│
├── raw_data/                          # Raw Nanopore FASTQ files (~13.0 GB)
├── illumina_raw_data/                 # Raw Illumina FASTQ files (~3.3 GB)
├── illumina_trimmed_data/             # Trimmed Illumina reads (pipeline output)
├── filtered_data/                     # Filtered Nanopore reads (pipeline output)
├── illumina_qc_reports/               # QC metrics & MultiQC reports for Illumina
├── qc_reports/                        # QC metrics & MultiQC reports for Nanopore
├── logs/                              # Execution stdout/stderr logs for debugging
│
└── shiver/                            # Local clone of the SHIVER assembly suite
    ├── bin/                           #   SHIVER executable scripts
    ├── data/                          #   Reference alignments and adapters
    └── docs/                          #   SHIVER manual and pipeline diagrams
```

> **Note:** Raw data, processed outputs, and log files are excluded from version control via `.gitignore`. Only scripts, environment configuration, and documentation are tracked.

---

## Installation & Environment Setup

### Prerequisites

- **Linux** (tested on Ubuntu via WSL and HPC clusters)
- **Conda** (Miniconda or Anaconda)
- **Slurm** workload manager (for HPC execution; scripts can also run locally)

### 1. Clone this repository

```bash
git clone <repository-url>
cd Internship
```

### 2. Create the Conda environment

```bash
conda env create -f HIV_U3analysis_env.yml
conda activate HIV_U3analysis
```

### 3. Verify tool installation

```bash
bash test_tools.sh
```

This script checks that all required tools (`NanoFilt`, `NanoStat`, `porechop_abi`, `trimmomatic`, `kraken2`, `NanoPlot`, `shiver_init.sh`, etc.) are available in the environment.

### 4. Environment dependencies

The `HIV_U3analysis_env.yml` installs the following from `bioconda` and `conda-forge`:

| Category | Tools |
|:---|:---|
| **General QC & Reporting** | FastQC, MultiQC |
| **Illumina Trimming** | Trimmomatic, fastp, Cutadapt |
| **Nanopore QC & Filtering** | NanoPlot, NanoQC, NanoStat, NanoFilt, Porechop_ABI |
| **Contamination Removal** | Kraken2 |
| **Data Retrieval & Manipulation** | SRA Toolkit (`prefetch`, `fasterq-dump`), SeqKit |
| **PacBio CCS** | pbccs |

---

## Usage

### Running on an HPC cluster (Slurm)

```bash
# Illumina preprocessing pipeline
sbatch qc_illumina_u3analysis.sh

# Nanopore preprocessing pipeline
sbatch qc_u3analyis.sh
```

### Running locally (without Slurm)

```bash
# Ensure the conda environment is activated
conda activate HIV_U3analysis

# Run directly
bash qc_illumina_u3analysis.sh
bash qc_u3analyis.sh
```

### Monitoring execution

```bash
# Check Slurm job status
squeue -u $USER

# View real-time logs
tail -f logs/slurm-<job-id>.out
```

---

## Pipeline Stages & Tools

The tools review document (*Tools_review_HIV_U3_analysis.pdf*) provides a comprehensive evaluation and justification for each tool selected at every stage:

### Stage 1 — Quality Control & Preprocessing

| Step | Illumina Tool | Nanopore Tool | Purpose |
|:---|:---|:---|:---|
| Raw QC | FastQC | NanoPlot, NanoQC, NanoStat | Assess read quality, length distribution, GC content |
| Adapter Trimming | Trimmomatic | Porechop_ABI | Remove sequencing adapters and chimeric reads |
| Quality Filtering | Trimmomatic (sliding window Q≥20, min length 50) | NanoFilt (Phred ≥7, length ≥200 bp) | Remove low-quality reads |
| Contamination | Kraken2 | Kraken2 | Remove human host and bacterial reads |
| Report Aggregation | MultiQC | MultiQC | Unified interactive HTML reports |

### Stage 2 — Genome Assembly

- **SHIVER** (Sequence/Haplotype Iterative Virus assemblER) — reference-guided iterative assembly using a curated set of HIV-1 reference sequences. Particularly suited for diverse viral populations.

### Stage 3 — Biological Filtering

- **Poplars / Hypermut 3** — Detects and excludes APOBEC3G/3F-induced G→A hypermutated sequences.
- **HIVSeqinR** — Identifies structural defects: major splice donor (MSD) mutations, Rev response element (RRE) deletions, and large internal deletions.

### Stage 4 — Multiple Sequence Alignment

- **MAFFT (L-INS-i)** — Progressive, iterative multiple alignment algorithm. Alignments are anchored to the **HXB2** reference genome (GenBank K03455) for consistent coordinate mapping.
- **Minimap2** — Rapid initial placement of long PacBio contigs.
- **MEGA** — Interactive alignment and tree exploration.

### Stage 5 — Subtype and Recombination Detection

- **jpHMM** — Resolves inter-subtype recombination (A1/D and CRF/URF breakpoints).
- **IQ-TREE 2** — Robust maximum-likelihood phylogenetic inference (-m MFP -B 1000) to confirm subtype calls.
- **COMET & REGA v3** — Fast first-pass classification of pure subtypes.
- **RIP** — Quick visual recombination screen.
- **ggtree** — Rendering of final annotated trees mapping molecular covariates.

### Stage 6 — U3 Region Extraction

- **SeqKit** — High-performance sequence manipulation tool used to extract the U3 sub-region based on HXB2 coordinate annotations (positions 1–454 of the 5' LTR).

### Stage 7 — Transcription Factor Binding Site (TFBS) Analysis

- **FIMO** (from the MEME Suite) — Scans extracted U3 sequences against **JASPAR 2024** position weight matrices (PWMs) for key host transcription factors:
  - **NF-κB (p50/p65)** — Master regulator of HIV-1 transcription; binding site copy number varies by subtype.
  - **SP1** — Basal promoter element; typically 3 binding sites in the U3 core promoter.
  - **NFAT** — T-cell activation-dependent factor.
  - **AP-1 (FOS/JUN), TBP, TFIID** — Additional regulatory factors evaluated.
- **TFBSTools** — Used in parallel within R to feed motif predictions directly into downstream statistical models without text parsing.

### Stage 8 — G-Quadruplex Prediction

- **R `gquad` package & `pqsfinder`** — Predicts intramolecular G-quadruplex (G4)-forming sequences in the U3 region. G4 structures are non-canonical DNA secondary structures that can modulate LTR promoter activity and influence viral latency.

---

## Key Design Decisions

### HPC Compatibility Fixes
- **Conda activation in non-interactive Slurm shells:** Scripts dynamically locate and source `$HOME/miniconda3/etc/profile.d/conda.sh` before activating the environment, as default `conda activate` commands fail in non-interactive batch job shells.
- **`set -o pipefail`:** Enabled in all pipeline scripts to correctly propagate failures through pipe chains (e.g., `tool 2>&1 | tee log.txt`). Without this, `tee` always returns exit code 0, silently masking tool failures.

### Idempotent / Resumable Execution
- All pipeline steps check for the existence of output files before running, allowing scripts to be safely re-executed after partial failures without reprocessing completed samples.

### Data Hygiene
- Corrupted zero-byte or minimal-size files from failed runs are detected and cleaned to prevent false "already completed" skips.
- SRA cache directories are automatically cleaned after successful FASTQ conversion to conserve disk space.

---

## References

1. **Research Proposal:** `Jovita_HIV1_U3_TFBS_Research_proposal_Draftv2.pdf` — Full study rationale, objectives, and expected outcomes.
2. **Tools Review:** `Tools_review_HIV_U3_analysis.pdf` — Comprehensive evaluation and justification of all bioinformatics tools used at each pipeline stage.
3. Struck, D., Lawyer, G., Ternes, A.-M., et al. (2014). *COMET: adaptive context-based modeling for ultrafast HIV-1 subtype identification.* Nucleic Acids Res.
4. Wymant, C., Blanquart, F., Golubchik, T., et al. (2018). *Easy and accurate reconstruction of whole HIV genomes from short-read sequence data with shiver.* Virus Evol.
5. Lee, G.Q., Reddy, K., Engelbrecht, S., et al. (2021). *NanoHIV — Oxford Nanopore GridION sequencing of near-full-length HIV-1 genomes.* NCBI BioProject PRJNA765218.
6. Grant, T.M., Ciccone, E.J., et al. (2013). *Prevalence and Clinical Impacts of HIV-1 Intersubtype Recombinants in Uganda.* NCBI BioProject PRJNA207834.
7. Bailey, T.L., et al. (2015). *The MEME Suite.* Nucleic Acids Res.

---

## License

This project is developed for academic research purposes as part of an MSc Bioinformatics internship at Makerere University.

---

*For questions or contributions, please contact the author.*
