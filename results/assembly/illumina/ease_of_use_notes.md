# Ease-of-use notes: assembly (Illumina)

Filled in 2026-07-21, after fixing and re-running all three tools on the
4-sample subset.

## BWA + bcftools consensus

- **Setup friction** (1-5, 1=trivial): 1 -- already the production
  pipeline's step 7, no issues.
- **Documentation quality**: Standard, well-known tools.
- **Failure modes encountered**: None in this session. A latent bug
  (missing explicit `bwa index` before mapping, silently relying on the
  production pipeline having done it already) was found and fixed
  earlier in the project, before this session.
- **Output usability**: Single consensus FASTA per sample, trivial to
  use downstream.
- **Would you use this again for the full cohort? Why/why not**: Yes --
  fast (3-48s per sample here), and valid full-length (9,719bp, 0% N)
  for all 4 subset samples.

## SPAdes (de novo assembly + best-contig-by-BLAST selection)

- **Setup friction** (1-5, 1=trivial): 4 -- `spades.py --careful`
  crashed on every subset sample with `Illegal instruction` (SIGILL,
  reported by Python's subprocess as "OS return value: -4"). Root
  cause: bioconda's default SPAdes 4.3.0 build requires AVX2/FMA
  instructions that this node's CPU doesn't have (same underlying
  issue as MUSCLE 5.3, see msa notes) -- specifically in
  `spades-hammer`'s k-mer counting, introduced by SPAdes 4.x's internal
  rewrite. Fixed by pinning `spades=3.15.5` (last pre-4.x release) in
  `HIV_U3analysis_env.yml`.
- **Documentation quality**: Good, well-established tool; the crash
  itself gave no useful diagnostic beyond the OS-level return code --
  had to reproduce directly and check `/proc/cpuinfo` flags to find the
  real cause.
- **Failure modes encountered**: The SIGILL crash above (now fixed).
  Separately, `--careful` mode's mismatch-correction step is slow and
  scales with read count more than linearly in practice --
  60K reads: ~2 min; 164K: ~8 min; 345K: ~10 min; 1.3M reads: ~44 min
  on this hardware (4 threads). This will matter for the full 24-sample
  cohort, some of which are larger still.
- **Output usability**: Produces many candidate contigs per sample;
  `run_spades.sh` already handles picking the best HXB2-matching one via
  BLAST, which works but needs a full assembly run (and its runtime
  cost) just to get one candidate sequence.
- **Assembly quality on this data -- the real finding**: even once
  fixed to run without crashing, SPAdes's best-BLAST-hit contig was
  **not consistently full-length** across the 4 subset samples: 1,111bp,
  992bp, and 16,697bp (too long, likely a chimeric/duplicated contig)
  were all outside the 8,000-10,000bp valid-length window; only the
  largest/highest-coverage sample (SRR908446, 1.3M reads) produced a
  valid-length assembly (9,228bp). This suggests de novo assembly of
  this metagenomic short-read HIV data is **coverage-sensitive** --
  lower-coverage samples fragment badly -- which is a real
  methodological data point, not a tool bug.
- **Would you use this again for the full cohort? Why/why not**: Only
  as a secondary/QC check, not as the primary assembly method for this
  cohort. BWA+consensus is faster, reference-anchored, and consistently
  produces full-length output regardless of coverage; SPAdes's
  reference-free approach is valuable for catching things BWA+consensus
  would miss by construction (novel insertions, structural variants) but
  its instability at lower coverage means it can't replace the
  reference-guided approach for this cohort's shallower samples.

## SHIVER

- **Setup friction**: Not yet assessable -- blocked on the one-time
  manual LANL reference-alignment download (see
  `shiver_setup/SOURCE.md`, still outstanding as of this session). All
  4 subset samples show the expected "reference alignment not found"
  failure, not a real SHIVER problem.
- **Would you use this again for the full cohort? Why/why not**: Can't
  assess until the manual download is done.

## Overall recommendation for this stage

**BWA + bcftools consensus remains the right primary choice** for this
cohort: fast, consistently full-length regardless of per-sample
coverage, and already validated across all 4 subset samples (and the
full 24-sample production run). SPAdes, once its hardware-compatibility
crash was fixed, revealed a genuine coverage-sensitivity limitation
that makes it unsuitable as a primary method here, though it remains
useful as a secondary check once assembly is reliable. SHIVER's
comparison is still pending the one-time manual reference download --
this is the last real blocker in the Illumina comparison harness.
