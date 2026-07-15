# HIV-1 Proviral U3 Region Analysis — Status Report

**Date:** 2026-07-15
**Prepared by:** Jovita Nagawa (MSc Bioinformatics, Makerere University)

## 1. Objective recap

Build a reproducible pipeline that processes Illumina (PRJNA207834, 24 samples,
subtypes A1/D/recombinant) and Nanopore (PRJNA765218, 9 samples) HIV-1 whole-genome
sequencing data through to 5' LTR U3 extraction and transcription-factor
binding-site (TFBS) / G-quadruplex motif mapping, with tool choices at each
step justified by a side-by-side comparison against the alternatives named in
`Tools_review_HIV_U3_analysis.pdf`, rather than taken on the review's word alone.

## 2. Work completed

### 2.1 Repository structure
- Reorganised the whole repo from a flat layout into `scripts/` (code),
  `data/` (raw + reference + processed), `results/` (QC reports + per-step
  comparison output), and `writeups/` (this document, the proposal, the
  tools review, progress reports).
- `scripts/` is split into `pipelines/` (production end-to-end scripts),
  `utils/`, `archive/` (superseded scripts), `tools/` (third-party clones:
  SHIVER, jpHMM, and — once cloned — Poplars/HIVSeqinR/HIVIntact), and the
  tool-comparison harness itself, organised as `scripts/<step>/<platform>/`
  (e.g. `scripts/assembly/illumina/`, `scripts/msa/oxfordnano/`) so the
  Nanopore side has an unambiguous home once built out further.
- Converted two Python glue scripts (U3 coordinate extraction, JASPAR
  matrix setup) to pure bash/awk, removing the Biopython dependency, so the
  whole `scripts/` tree is bash-only — verified byte-identical output
  against the originals before switching over.

### 2.2 Production pipelines
- `scripts/pipelines/illumina_u3analysis.sh` and `oxnano_u3analysis.sh`:
  monolithic, checkpointed, resumable end-to-end scripts. Steps 1-7
  (download → QC → trim/filter → reference-mapped consensus) and step 9
  (MAFFT MSA) + U3 extraction are fully implemented. **Steps 8 (biological
  filtering) and 10 (subtyping) are still stubs** (checkpoint touched, no
  real tool call), and motif/G-quadruplex mapping (rest of step 11) is not
  yet wired in — this is intentional, pending the tool-comparison results
  below feeding back into a decision.
- Started breaking the Illumina monolith into one script per step
  (`scripts/pipelines/illumina/01_download_qc_trim.sh`): downloads, runs
  FastQC, and now trims with **both** Trimmomatic and fastp as real
  production output (previously fastp was comparison-only), each held to
  the same Q20 quality bar (see §2.3), producing two separate MultiQC
  reports.

### 2.3 Tool-comparison harness
Six steps, each comparing every automatable candidate tool from the tools
review, on a fixed sample subset (4 Illumina, 3 Nanopore, chosen by
file-size quartile as a complexity proxy):

| Step | Compares | Platform coverage |
|---|---|---|
| `download_qc` | fastp vs Trimmomatic (Illumina); Porechop_ABI+NanoFilt validity check (Nanopore, no alternative exists) | both |
| `assembly` | BWA+bcftools-consensus vs SPAdes vs SHIVER (Illumina); minimap2 validity check (Nanopore) | both |
| `msa` | MAFFT vs MUSCLE vs Clustal Omega | both |
| `biological_filtering` | Poplars vs HIVSeqinR vs HIVIntact | both (not yet run, see §3) |
| `subtyping` | jpHMM vs IQ-TREE2 | both (not yet run, see §3) |
| `motif_mapping` | FIMO vs MOODS vs TFBSTools, with a positive-control check against HXB2's own U3 | both |

Illumina and Nanopore trimming thresholds are held to a documented,
justified Q20 bar (Phred Q = -10·log₁₀(P_error); Q20 = 99% base-call
accuracy) deliberately stricter than Trimmomatic's own textbook example
(Q15), since a single miscalled base can flip a motif match downstream.
Nanopore uses Q7/200bp (Nanopore's inherently higher per-base error rate
makes a Q20 floor unusable — would discard nearly all reads).

All six steps' scripts are consolidated (no separate `run_<tool>.sh` files
where feasible) with inline bash functions for multi-step tools like SHIVER,
timed and validity-checked via a shared `measure_and_run`/`summary.tsv`
harness (`scripts/common/lib_compare.sh`).

## 3. Results so far

Both production Slurm jobs (Illumina: job 109080, Nanopore: job 109081) are
currently running in the background. In parallel, the tool-comparison
harness has been run against the Illumina subset (Nanopore in progress):

- **download_qc (Illumina, n=4):** both trimmers ran successfully.
  Surviving-read rate varied a lot by sample (Trimmomatic: 12.3%-63.2%),
  reflecting genuinely different per-sample input quality, not a tool
  difference — needs a same-sample fastp-vs-Trimmomatic comparison to be
  meaningful (currently just recorded as separate metrics).
- **assembly (Illumina, n=4):** BWA+bcftools-consensus succeeded for all 4
  samples (consistent 9,719bp, 0% N consensus — full-length, no gaps).
  SPAdes failed on all 4 (see Challenges). SHIVER failed on all 4 as
  expected — its one-time manual reference-alignment download hasn't been
  done yet.
- **msa (Illumina, n=1 alignment of 5 sequences: HXB2 + 4 consensus):**
  MAFFT and Clustal Omega both succeeded (9,719 columns, 0% gaps — the 4
  consensus sequences are essentially identical in length to HXB2 here).
  MUSCLE failed (see Challenges).
- **subtyping (Illumina):** jpHMM failed on all 4 samples (binary not yet
  compiled — one-time setup not run). IQ-TREE2 ran but found the alignment
  invalid for tree-building (see Challenges — likely too few/too similar
  sequences in this 5-sequence subset, not a real IQ-TREE2 problem).
- **motif_mapping (Illumina):** U3 extraction succeeded (453bp region
  extracted from all 5 aligned sequences, anchored to HXB2's own annotated
  LTR/R-region boundaries, no outliers flagged). **All three motif
  scanners (FIMO, MOODS, TFBSTools) found zero hits, including on HXB2's
  own positive control** — see Challenges, root cause identified.
- **biological_filtering (Illumina):** not yet run — one-time tool cloning
  (`setup_tools.sh`) hasn't been executed.

## 4. Challenges encountered (with root causes, where found)

1. **Motif scanning found zero hits everywhere, including the positive
   control — root cause identified.** Traced this back through the
   pipeline: MAFFT's aligned output is lowercase (`results/msa/illumina/mafft_aligned.fasta`
   has `tggaagggc...` where the input reference `data/reference/K03455.1.fasta`
   is uppercase `TGGAAGGGC...`). This lowercase case propagates through U3
   extraction into the sequences handed to FIMO/MOODS/TFBSTools, and all
   three appear to treat lowercase as masked/invalid rather than matching
   case-insensitively. **This is a quick fix** (uppercase the sequence
   before motif scanning, e.g. `seqkit seq -u`) but wasn't caught until
   results were actually inspected — a good example of why the positive
   control step earns its keep.
2. **SPAdes crashes with an internal error** (`spades-hammer` exits with
   OS return value -4, a segfault-class failure) on every Illumina sample
   tested. This looks like an environment/binary compatibility issue
   (SPAdes's bundled binaries vs. this system), not a data or script
   problem — needs investigation independent of this pipeline (try a
   different SPAdes build/version, or run on a different node).
3. **One-time manual/setup steps not yet done**, so three comparison arms
   currently show expected, not real, failures:
   - SHIVER needs a one-time manual LANL reference-alignment download
     (`scripts/assembly/illumina/shiver_setup/SOURCE.md`).
   - jpHMM needs compiling from source (`setup_jphmm.sh`) — not yet run.
   - Poplars/HIVSeqinR/HIVIntact need cloning (`setup_tools.sh`) — not yet
     run, and HIVSeqinR additionally needs manual RStudio configuration
     before its first real run.
4. **A latent `bwa index` bug surfaced and was fixed during this run**:
   the BWA+consensus assembly path never explicitly indexed the reference,
   silently relying on the production pipeline having done it first. Once
   actually run standalone (as the comparison harness does), it failed
   instantly for every sample until this was added.
5. **MUSCLE and IQ-TREE2 both failed on this specific 5-sequence
   subset** — not yet root-caused; plausible that MUSCLE's CLI invocation
   or IQ-TREE2's model-selection step needs more sequence diversity than
   4 near-identical-to-HXB2 consensus sequences provide. Needs the fuller
   sample set (or investigation of the actual error logs) to distinguish
   a real tool problem from a too-small/too-similar test subset.
6. **Nanopore NanoQC/NanoStat are very slow** on this system (observed
   ~18 minutes for NanoQC on one sample in the production run )— worth
   profiling if this becomes a bottleneck for the full 9-sample run.

## 5. Recommendations for future work

1. **Immediate, cheap fix:** uppercase sequences before motif scanning
   (`scripts/motif_mapping/*/compare_motif_mapping_*.sh`) and re-run —
   this one change should resolve the zero-hit result across all three
   tools and finally produce a real FIMO/MOODS/TFBSTools comparison.
2. **Run the outstanding one-time setup steps** (`setup_tools.sh`,
   `setup_jphmm.sh`, SHIVER's manual reference download) so
   biological_filtering and subtyping produce real comparative results
   instead of expected-failure placeholders.
3. **Investigate the SPAdes crash** on this system independently of the
   pipeline logic — likely an environment/binary issue worth resolving
   before drawing conclusions about SPAdes vs BWA vs SHIVER for this
   thesis's methodology section.
4. **Re-run MSA/subtyping comparisons against the full sample set**
   (not just the 4/3-sample subset) once assembly succeeds for more
   samples, to check whether MUSCLE/IQ-TREE2's failures were subset-size
   artifacts or real tool issues.
5. **Extend the same-sample comparison** for download_qc: currently fastp
   and Trimmomatic are compared as separate summary rows per sample, not
   a paired same-input comparison — worth adding a direct retained-read
   and downstream-mapping-rate comparison once more samples are assembled
   with each trimmer's output.
6. **Wire the winning tools into the production pipelines** (steps 8, 10,
   and the rest of 11) once the comparisons above are complete and a
   defensible choice can be made per step, replacing the current stubs.
7. **Build the remaining Nanopore comparison scripts' first real run**
   (msa/biological_filtering/subtyping/motif_mapping/oxfordnano) once its
   assembly step completes — these reuse the same tools as the Illumina
   side, just against Nanopore-derived consensus sequences, and were
   built but not yet exercised end-to-end at the time of this report.
8. **Complete the per-step production pipeline breakdown** started for
   Illumina step 1 (`scripts/pipelines/illumina/`) for the remaining
   steps, and build the equivalent for Nanopore, once the comparison
   harness has settled on winning tools for each step.
