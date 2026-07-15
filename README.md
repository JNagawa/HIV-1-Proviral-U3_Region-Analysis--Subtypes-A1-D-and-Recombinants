# HIV-1 Proviral U3 Region Analysis Pipeline

[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20HPC%20(Slurm)-blue)]()
[![Conda](https://img.shields.io/badge/Env-Conda-green)]()

A reproducible bioinformatics pipeline for processing, assembling, and analyzing the U3 region of HIV-1 proviral genomes. This repository contains scripts and environments to process both short-read (Illumina) and long-read (Oxford Nanopore/PacBio) sequencing data, identify intact genomes, and systematically map host transcription factor binding sites (TFBS) and G-quadruplex structures.

## Repository Structure

```text
/home/jnagawa/Internship/
├── HIV_U3analysis_env.yml        # Conda environment specification
├── scripts/
│   ├── pipelines/                # illumina_u3analysis.sh, oxnano_u3analysis.sh (Slurm/Bash),
│   │                             # plus the newer per-step illumina/01_download_qc_trim.sh etc.
│   ├── utils/                    # extract_u3_by_hxb2_anchor.sh, test_tools.sh
│   ├── archive/                  # superseded scripts, variant_call.sh (unrelated coursework script)
│   ├── tools/                    # third-party clones: shiver/, jpHMM/, Poplars/, HIVSeqinR/, HIVIntact/
│   ├── common/                   # tool-comparison harness's shared lib_compare.sh + its own README
│   └── download_qc_illumina/, assembly_illumina/, msa_illumina/, biological_filtering_illumina/,
│       subtyping_illumina/, motif_mapping_illumina/
│                                 # tool-comparison harness, one step per folder (see scripts/common/README.md);
│                                 # _illumina since there's no Nanopore side yet (would be _oxnano, matching below)
├── data/
│   ├── raw/{illumina,oxnano}/    # downloaded FASTQ
│   ├── reference/                # HXB2 (K03455.1) fasta + cached GenBank record + jaspar/ PWM downloads
│   └── processed/{illumina,oxnano}/<stage>/  # trimmed/filtered, alignments, msa, subtyping, motifs
├── results/
│   ├── reports/qc/{illumina,oxnano}/  # FastQC/MultiQC/NanoPlot/etc. QC reports
│   ├── figures/                  # plots generated during analysis
│   └── download_qc_illumina/, assembly_illumina/, msa_illumina/, biological_filtering_illumina/,
│       subtyping_illumina/, motif_mapping_illumina/
│                                 # tool-comparison harness output: per-step summary.tsv, timing logs, ease_of_use_notes.md
├── writeups/                     # proposal, tools review, progress report
├── checkpoints/                  # pipeline resume markers
└── logs/                         # Slurm stdout/err
```

## Pipeline Overview

The analysis workflow consists of the following key steps:
1. **Quality Control & Preprocessing:** FastQC/MultiQC, Trimmomatic (Illumina), NanoPlot/Porechop/NanoFilt (Nanopore), and Kraken2 (Contamination removal).
2. **Genome Assembly:** Reference-guided iterative assembly using **SHIVER**.
3. **Biological Filtering:** Hypermutation detection (Poplars) and structural intactness classification (HIVSeqinR).
4. **Alignment & Subtyping:** Multiple sequence alignment via **MAFFT**, with subtype/recombination detection using **jpHMM** and **IQ-TREE 2**.
5. **U3 Extraction & Motif Mapping:** Coordinate-based extraction with **SeqKit**, TFBS scanning with **FIMO** & **TFBSTools**, and G-Quadruplex prediction via **gquad** and **pqsfinder**.

## Installation & Setup

1. **Clone the repository and navigate to the directory:**
   ```bash
   git clone <repository-url>
   cd Internship
   ```

2. **Create and activate the Conda environment:**
   ```bash
   conda env create -f HIV_U3analysis_env.yml
   conda activate HIV_U3analysis
   ```

## Usage

The primary preprocessing pipelines are optimized for High-Performance Computing (HPC) environments using the Slurm workload manager, but can also be run locally.

**Running via Slurm (from the repo root):**
```bash
sbatch scripts/pipelines/illumina_u3analysis.sh
sbatch scripts/pipelines/oxnano_u3analysis.sh
```

**Running Locally (from the repo root):**
```bash
bash scripts/pipelines/illumina_u3analysis.sh
bash scripts/pipelines/oxnano_u3analysis.sh
```

> **Note - Smart Resuming:** If the pipeline is interrupted, you can safely run the script again. It automatically checks for files that were already processed and picks up exactly where it left off, saving time.

## License
Developed for academic research purposes as part of an MSc Bioinformatics internship at Makerere University.
