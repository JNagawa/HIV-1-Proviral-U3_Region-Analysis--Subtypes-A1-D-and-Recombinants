# HIV-1 Proviral U3 Region Analysis Pipeline

[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20HPC%20(Slurm)-blue)]()
[![Conda](https://img.shields.io/badge/Env-Conda-green)]()

A reproducible bioinformatics pipeline for processing, assembling, and analyzing the U3 region of HIV-1 proviral genomes. This repository contains scripts and environments to process both short-read (Illumina) and long-read (Oxford Nanopore/PacBio) sequencing data, identify intact genomes, and systematically map host transcription factor binding sites (TFBS) and G-quadruplex structures.

## 📂 Repository Structure

```text
/home/jnagawa/Internship/
├── HIV_U3analysis_env.yml        # Conda environment specification
├── qc_illumina_u3analysis.sh     # Illumina preprocessing pipeline (Slurm/Bash)
├── qc_oxnano_u3analysis.sh       # Nanopore preprocessing pipeline (Slurm/Bash)
├── variant_call.sh               # Variant calling pipeline (BWA, SAMtools, BCFtools, SnpEff)
├── shiver/                       # Local clone of the SHIVER assembly suite
└── [Data Directories]            # Generated during pipeline execution (raw_data, qc_reports, etc.)
```

## ⚙️ Pipeline Overview

The analysis workflow consists of the following key stages:
1. **Quality Control & Preprocessing:** FastQC/MultiQC, Trimmomatic (Illumina), NanoPlot/Porechop/NanoFilt (Nanopore), and Kraken2 (Contamination removal).
2. **Genome Assembly:** Reference-guided iterative assembly using **SHIVER**.
3. **Biological Filtering:** Hypermutation detection (Poplars) and structural intactness classification (HIVSeqinR).
4. **Alignment & Subtyping:** Multiple sequence alignment via **MAFFT**, with subtype/recombination detection using **jpHMM** and **IQ-TREE 2**.
5. **U3 Extraction & Motif Mapping:** Coordinate-based extraction with **SeqKit**, TFBS scanning with **FIMO** & **TFBSTools**, and G-Quadruplex prediction via **gquad** and **pqsfinder**.

## 🚀 Installation & Setup

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

## 💻 Usage

The primary preprocessing pipelines are optimized for High-Performance Computing (HPC) environments using the Slurm workload manager, but can also be run locally.

**Running via Slurm:**
```bash
sbatch qc_illumina_u3analysis.sh
sbatch qc_oxnano_u3analysis.sh
```

**Running Locally:**
```bash
bash qc_illumina_u3analysis.sh
bash qc_oxnano_u3analysis.sh
```

> **Note - Smart Resuming:** If the pipeline is interrupted, you can safely run the script again. It automatically checks for files that were already processed and picks up exactly where it left off, saving time.

## 📄 License
Developed for academic research purposes as part of an MSc Bioinformatics internship at Makerere University.
