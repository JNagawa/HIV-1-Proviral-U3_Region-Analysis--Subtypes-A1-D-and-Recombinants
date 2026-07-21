# Ease-of-use notes: msa (Illumina)

Filled in 2026-07-21, after fixing MUSCLE (see below). All three tools
now align the same 5-sequence subset (HXB2 + 4 consensus) to identical
results: 5/5 sequences retained, 9,719 columns, 0% gaps -- unsurprising
given how close these early-subset consensus sequences are to HXB2
itself (see assembly stage notes); a true differentiator between the
three tools' alignment quality will need the full, more diverse 24-sample
cohort.

## MAFFT

- **Setup friction** (1-5, 1=trivial): 1 -- already a production-pipeline
  dependency, no issues.
- **Documentation quality**: Excellent, well-known tool.
- **Failure modes encountered**: None. Fastest correct result (12.7s).
- **Output usability**: Standard aligned FASTA, straightforward.
- **Would you use this again for the full cohort? Why/why not**: Yes --
  already the production choice, no reason found here to reconsider it.

## MUSCLE

- **Setup friction** (1-5, 1=trivial): 4 -- the bioconda default (5.3)
  crashes with `Illegal instruction` (SIGILL, exit 132) on this node's
  CPU: 5.3's build requires AVX2/FMA, but this hardware only has
  AVX+SSE4.2. Downgrading to classic MUSCLE 3.8 (`muscle=3.8*` in
  `HIV_U3analysis_env.yml`) fixed it, but v3's CLI is different from
  v5's (`-in`/`-out`, single-threaded, no `-align`/`-output`/`-threads`)
  -- `run_muscle.sh` needed updating, not just a version bump.
- **Documentation quality**: Classic MUSCLE 3.8's docs are older but
  clear; the CLI is simple.
- **Failure modes encountered**: The SIGILL crash above, on *every*
  invocation regardless of input -- this is a hardware/build
  compatibility issue, not something retry or different flags would
  fix. Confirmed via direct reproduction (`/proc/cpuinfo` flags vs.
  the build's AVX2 requirement).
- **Output usability**: Standard aligned FASTA once working.
- **Would you use this again for the full cohort? Why/why not**: Yes,
  now that it runs -- but this is a hardware-specific fix. If this
  pipeline is ever run on different hardware with AVX2 available, the
  version pin to 3.8 should be revisited (5.3 is faster and more
  actively maintained) rather than left in place indefinitely.

## Clustal Omega

- **Setup friction** (1-5, 1=trivial): 1 -- no issues.
- **Documentation quality**: Good.
- **Failure modes encountered**: None, but by far the slowest of the
  three (108-116s vs. 12-15s for MAFFT/MUSCLE) on this tiny 5-sequence
  subset -- worth watching for the full 24-sample cohort, where this
  gap could become a real runtime cost.
- **Output usability**: Standard aligned FASTA.
- **Would you use this again for the full cohort? Why/why not**:
  Possibly, as a cross-check, but its runtime scaling should be
  confirmed on a larger sample before committing to it for routine use.

## Overall recommendation for this stage

**MAFFT remains the right choice** -- already the production default,
fastest correct result here, and no issues found. MUSCLE now works as a
cross-check after the 3.8 downgrade, but that downgrade is a
hardware-specific workaround (this node lacks AVX2) that should be
revisited if the pipeline ever runs on different hardware. Clustal
Omega's ~9x slower runtime on even this small subset is a real concern
for the full cohort and should be re-measured before relying on it at
scale.
