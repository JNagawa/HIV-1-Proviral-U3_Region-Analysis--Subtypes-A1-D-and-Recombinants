# Ease-of-use notes: subtyping (Illumina)

Filled in 2026-07-21.

## jpHMM

- **Setup friction** (1-5, 1=trivial): 4 -- no packaged distribution
  (conda/pip); ships as C++ source with a bundled Boost copy and must be
  compiled from a small academic server's tarball
  (`http://jphmm.gobics.de/jpHMM.tar.gz`, plain HTTP). Compiled cleanly
  with g++ 14.2.0 (some old-Boost/new-g++ warnings, no errors), but the
  resulting binary then failed at runtime with
  `GLIBCXX_3.4.32 not found` -- this HPC's Lmod-managed
  `LD_LIBRARY_PATH` puts `/opt/ohpc/pub/apps/anaconda3/lib` (an older
  libstdc++) ahead of the active conda env's own lib dir, which *does*
  have the needed GLIBCXX version once `libstdcxx-ng` was installed --
  it was just being shadowed. Fixed by prepending
  `${CONDA_PREFIX}/lib` to `LD_LIBRARY_PATH` inside `run_jphmm.sh`
  itself (scoped to this invocation, not a global env change).
- **Documentation quality**: CLI flags aren't documented anywhere
  external; had to read `src/main.cpp`'s `getopt()` call directly to
  confirm real usage.
- **Failure modes encountered**: The GLIBCXX shadowing above. Also
  slow -- ~2.3 minutes per sample on this subset.
- **Output usability**: Multiple output files per sample (breakpoint
  positions, posterior probabilities, GFF3 Viterbi path); the
  `recombination.txt` file is the most directly useful summary.
- **A finding worth flagging, not just an infra note**: all 4 subset
  samples were called with **identical** breakpoints (1-789
  5'-Insertion, 790-9411 subtype B, 9412-9719 3'-Insertion) despite
  being different samples. This cohort (PRJNA207834) is documented as
  A1/D/recombinant, not B -- so this is very likely a **reference-bias
  artifact from the BWA+consensus assembly method**, not a genuine
  finding: reference-guided consensus calling fills low-confidence
  positions toward HXB2 (itself subtype B), so near-identical
  consensus sequences across samples would produce near-identical
  (mis)calls regardless of the samples' true subtype. **This needs
  cross-checking against the full 24-sample cohort and ideally a
  non-reference-biased assembly (e.g. SPAdes, once it can produce a
  usable near-full-length contig) before trusting any subtype call
  from this comparison.**
- **Would you use this again for the full cohort? Why/why not**: Yes,
  once compiled and the library path is fixed -- but the identical-call
  finding above means its output on *this* subset shouldn't be taken at
  face value; it's flagging an assembly-stage problem more than
  answering the subtyping question.

## IQ-TREE2

- **Setup friction** (1-5, 1=trivial): 1 -- already resolved via the
  `iqtree=2.*` pin in `HIV_U3analysis_env.yml`.
- **Documentation quality**: Excellent, well-established tool.
- **Failure modes encountered**: Refused to run ultrafast bootstrap
  ("It makes no sense to perform bootstrap with less than 4 sequences")
  on this 5-sequence (HXB2 + 4 near-identical consensus) subset --
  consistent with the same reference-bias issue noted above: if the 4
  consensus sequences are nearly identical to each other and to HXB2,
  IQ-TREE effectively has too few distinct taxa for bootstrap to be
  meaningful. Dropped `-B 1000` for this small-subset comparison run
  (see `run_iqtree.sh`); re-add it once run against the full,
  genuinely diverse cohort.
- **Output usability**: Standard Newick treefile plus a detailed
  `.iqtree` report; easy to parse.
- **Would you use this again for the full cohort? Why/why not**: Yes --
  works correctly, the only issue was subset size/diversity, not the
  tool itself.

## Overall recommendation for this stage

Both tools now run correctly, but **this subset's results shouldn't be
over-interpreted**: jpHMM's identical subtype-B call across all 4
samples and IQ-TREE2's refusal to bootstrap both point to the same root
cause -- the BWA+consensus assembly method producing consensus
sequences too close to the HXB2 reference to distinguish real subtype
signal from reference bias. This is a genuine methodological flag for
the full-cohort run, not just a small-subset artifact to shrug off:
before trusting subtype calls at scale, either confirm the full
24-sample cohort has genuinely more divergent consensus sequences, or
cross-validate against a non-reference-biased assembly path.
