# PROGRESS REPORT: HIV-1 PROVIRAL U3 SEQUENCE ANALYSIS PIPELINE

**Prepared by:** Jovita Nagawa  
**Date:** July 8, 2026  
**Project Objective:** To extract and annotate U3 sequences from intact near-full-length HIV-1 proviral genomes (FLIP-SEQ and HIV SMRTcap) in the Rakai Cohort, stratified by subtype (A1, D, CRF/URF), and to systematically map binding sites for host transcription factors and G-quadruplexes.

---

## 1. Project Background and Objective

The U3 region of the HIV-1 Long Terminal Repeat (LTR) promoter is highly variable and plays a critical role in viral transcription, replication, and latency. Different HIV-1 subtypes (e.g., A1, D, and various recombinant forms common in East Africa) exhibit distinct transcription factor binding site (TFBS) profiles. 

This project aims to:
1. Process raw high-throughput sequencing data (short-read Illumina and long-read Oxford Nanopore/PacBio).
2. Perform quality control and trim adapter/contaminant sequences.
3. Assemble near-full-length viral genomes.
4. Filter out defective proviruses (APOBEC-hypermutated or structurally deleted genomes).
5. Align intact genomes, extract the U3 region, and map regulatory motifs (TFBS and G-quadruplexes).

---

## 2. Nature and Structure of the Datasets

To conduct this comparative study, we acquired two major publicly available sequence datasets representing short-read (Illumina) and long-read (Nanopore) sequencing technologies.

### 2.1 Short-Read Illumina Dataset (BioProject: PRJNA207834)
*   **Study Title:** *Prevalence and Clinical Impacts of HIV-1 Intersubtype Recombinants in Uganda Revealed by Near-Full-Genome Population and Deep Sequencing*
*   **Source Population:** Treatment-naive individuals from rural Mbarara, Uganda.
*   **Viral Characteristics:** Covers HIV-1 subtypes A1, D, and A1-D intersubtype recombinants, representing the exact genetic diversity profile of the Rakai cohort.
*   **Sequencing Platform:** Illumina MiSeq, generating **2x251 bp paired-end reads** (metagenomic whole-genome sequencing covering the 5' LTR).
*   **Data Size & File Structure:**
    *   **Number of Samples:** 24 paired-end samples (`SRR908430` to `SRR908453`).
    *   **Directory Location:** `illumina_raw_data/`
    *   **File Format:** Gzipped FASTQ format. Each sample consists of two files:
        *   `SRRXXXXXX_1.fastq.gz` (Forward reads)
        *   `SRRXXXXXX_2.fastq.gz` (Reverse reads)
    *   **Total Data Volume:** Approximately 3.3 GB of raw compressed sequence data.

### 2.2 Long-Read Nanopore Dataset (BioProject: PRJNA765218)
*   **Study Title:** *NanoHIV - Oxford Nanopore GridION sequencing of near-full-length HIV-1 genomes*
*   **Source Institution:** Stellenbosch University, South Africa.
*   **Viral Characteristics:** Near-full-length HIV-1 proviral genomes (~9 kb).
*   **Sequencing Platform:** Oxford Nanopore Technologies (ONT) GridION, producing single-end long reads. Long-read datasets are crucial for this study because they can span the entire proviral genome in a single read, resolving complex insertions, deletions, and LTR duplications that confuse short-read assemblers.
*   **Data Size & File Structure:**
    *   **Number of Samples:** 9 single-end samples (`SRR16005710` to `SRR16005718`).
    *   **Directory Location:** `raw_data/`
    *   **File Format:** Gzipped FASTQ format.
        *   `SRRXXXXXX.fastq.gz` (Single-end ONT long reads)
    *   **Total Data Volume:** Approximately 13.0 GB of raw compressed sequence data.

---

## 3. Computational Environment & Tools Setup

A dedicated Conda environment named `HIV_U3analysis` has been established to ensure complete reproducibility of the pipeline. The environment configuration is stored in `HIV_U3analysis_env.yml` and includes the following dependencies:

| Category | Tool | Function in Pipeline |
| :--- | :--- | :--- |
| **QC & Aggregation** | `FastQC` | Raw read quality metrics (Illumina) |
| | `NanoPlot` / `NanoStat` / `nanoQC` | Read length, quality distributions, and stats (Nanopore) |
| | `MultiQC` | Interactive HTML reports aggregating QC runs |
| **Trimming & Filtering** | `Trimmomatic` | Adapter and quality-based sliding window trimming (Illumina) |
| | `porechop_abi` / `porechop` | Nanopore adapter trimming and chimera splitting |
| | `NanoFilt` | Quality-score (Phred > 7) and length (> 200 bp) filtering for Nanopore |
| **Contamination Removal** | `Kraken2` | Classification and removal of human host/bacterial reads |
| **Manipulation** | `SRA Toolkit` (`prefetch`, `fasterq-dump`) | Downloading and converting SRA accessions to FASTQ |
| | `SeqKit` | Sequence manipulation and U3 sub-sequence coordinate extraction |

---

## 4. Pipeline Development and Implementation Status

Two automated shell pipelines designed for High-Performance Computing (HPC) environments running the Slurm workload manager have been implemented:

1.  **`qc_illumina_u3analysis.sh`**: Handles SRA download, validation, conversion, pre-trimming FastQC, Trimmomatic quality/adapter trimming, post-trimming FastQC, and MultiQC reporting for the 24 Illumina samples.
2.  **`qc_u3analyis.sh`**: Handles SRA download, validation, conversion, pre-filtering NanoPlot/NanoQC/NanoStat, Porechop adapter trimming, NanoFilt filtering, post-filtering NanoPlot/NanoStat, and MultiQC reporting for the 9 Nanopore samples.
3.  **`test_tools.sh`**: A diagnostic script used to verify tool availability and path configurations in the conda environment.

### 4.1 Critical Script Adjustments Made for HPC Compatibility
During initial test runs, we identified and corrected two critical bugs to guarantee the scripts run reliably on an HPC cluster:
*   **HPC Conda Activation Fix:** Sourcing `conda.sh` using default commands often fails in non-interactive Slurm shells. We updated the scripts to dynamically locate and source the environment profile (`$HOME/miniconda3/etc/profile.d/conda.sh` or standard modules) before activating the `HIV_U3analysis` environment.
*   **Pipe Failure Capture (`set -o pipefail`):** Because QC outputs are piped to `tee` to create logs (e.g., `tool 2>&1 | tee log.txt`), standard bash scripts return the exit status of `tee` (always `0`), hiding tool failures. We activated `set -o pipefail` so that failures during trimming (e.g. `trimmomatic` or `porechop_abi`) are correctly caught by `check_exit` functions and fallback commands.
*   **Data Cleanup:** Removed corrupted 20-byte gzipped files generated in `filtered_data/` during failed initial runs to prevent the script from skipping them in future runs.

---

## 5. Directory Structure Map

The work folder structure is organized as follows:

```text
/home/jnagawa/Internship/
+-- HIV_U3analysis_env.yml        # Conda environment YAML configuration
+-- test_tools.sh                 # Environment tool validation script
+-- qc_illumina_u3analysis.sh     # Illumina preprocessing pipeline (Slurm/Bash)
+-- qc_u3analyis.sh               # Nanopore preprocessing pipeline (Slurm/Bash)
+-- raw_data/                     # Raw Nanopore FASTQ files (13.0 GB)
+-- illumina_raw_data/            # Raw Illumina FASTQ files (3.3 GB)
+-- illumina_trimmed_data/        # Directory for trimmed Illumina reads (outputs)
+-- filtered_data/                # Directory for filtered Nanopore reads (outputs)
+-- illumina_qc_reports/          # QC metrics and MultiQC reports for Illumina
+-- qc_reports/                   # QC metrics and MultiQC reports for Nanopore
+-- logs/                         # Execution stdout and stderr logs for debugging
+-- shiver/                       # Local clone of the SHIVER assembly suite
```

---

## 6. Next Steps

With raw datasets successfully downloaded, validated, and the QC scripts refactored for robust HPC execution, we are prepared to execute the following stages of the workflow:

1.  **Run Preprocessing Pipelines**: Execute the pre-processing scripts on the HPC cluster using `sbatch` to obtain clean, trimmed reads.
2.  **Genome Assembly**: Assemble the trimmed Illumina reads into draft consensus sequences using `SHIVER` (Stage 2 in the review document).
3.  **Biological Filtering & Intactness**: Apply `Poplars` (Hypermut 3) and `HIVSeqinR` on the assemblies to exclude genomes with APOBEC-induced G-to-A hypermutations or significant structural deletions (MSD or RRE mutations) (Stage 3).
4.  **Align & Extract**: Align the confirmed intact genomes using `MAFFT` anchored to `HXB2` coordinates (Stage 4) and extract the U3 region coordinates using `SeqKit` (Stage 5).
5.  **Motif & G-quadruplex Analysis**: Scan the extracted promoter sequences for host transcription factor binding sites (NF-?B, SP1, NFAT) using `FIMO` and predict G-quadruplex-forming regions using the `gquad` package in R.
