# Ease-of-use notes: biological_filtering (Illumina)

Filled in 2026-07-21. All three tools are GitHub-hosted, not
conda/PyPI-packaged, so setup friction is inherently higher here than
for the other stages -- see `setup_tools.sh` for the clone/install steps.

## Poplars (Hypermut 3)

- **Setup friction** (1-5, 1=trivial): 3 -- `pip install -e` worked
  cleanly, but the package ships with no `__main__.py` and no
  `console_scripts` entry point despite being pip-installable, so
  `python -m poplars hypermut` (the invocation its own source structure
  suggests) fails with "'poplars' is a package and cannot be directly
  executed." The real entry point is running `poplars/hypermut.py`
  directly as a script.
- **Documentation quality**: Thin -- no CLI usage doc; had to read
  `hypermut.py`'s own `argparse` block to find real usage.
- **Failure modes encountered**: A real upstream bug -- `rate_ratio()`
  returns the string `"undef"` on division-by-zero (no potential
  control mutation sites in this data), but both `pretty_print()` and
  `make_data_file()` unconditionally `round()` it, crashing with
  `TypeError: type str doesn't define __round__ method`. Patched both
  call sites locally in the vendored copy
  (`scripts/tools/Poplars/poplars/hypermut.py`) -- **note this as a
  local patch to an upstream bug in the methods write-up**, not stock
  Poplars behavior.
- **Output usability**: Plain CSV-in-txt, easy to parse; includes a
  hypermutation call, rate ratio, and Fisher's exact p-value per
  sequence.
- **Would you use this again for the full cohort? Why/why not**: Yes --
  once patched, ran cleanly and fast (7.65s, 9,527 result rows across
  the subset). The `"undef"` rate-ratio case is likely to recur on the
  full cohort too (any near-reference consensus with zero APOBEC-context
  control sites), so the patch needs to travel with the tool, not be a
  one-off.

## HIVSeqinR

- **Setup friction** (1-5, 1=trivial): 5 -- not a CLI tool at all. It's
  an R script meant to be opened in RStudio and run interactively
  ("highlight all, run all"), with the BLAST DB path and 2nd-round PCR
  primer sequences hardcoded at the top of
  `R_HIVSeqinR_Combined_ver04.R` rather than passed as arguments.
- **Documentation quality**: README describes the RStudio workflow but
  doesn't document a way to drive it non-interactively.
- **Failure modes encountered**: `run_hivseqinr.sh` correctly detects
  and reports the missing manual configuration (gated on a
  `.CONFIGURED` marker file) rather than attempting to guess at
  BLAST DB paths / primers -- this is the right call for a tool that
  classifies genome intactness, but it means **this stage is still
  blocked on a manual step from the user** (edit the R script in
  RStudio, then `touch scripts/tools/HIVSeqinR/.CONFIGURED`).
  Separately, it requires input FASTA with no IUPAC ambiguity codes;
  bcftools consensus legitimately produces these at heterozygous sites,
  so `run_hivseqinr.sh` resolves each ambiguity code to its
  alphabetically-first base as a workaround -- **a real limitation to
  note in the methods write-up, not a transparent substitution**.
- **Output usability**: Not yet observed (blocked on config).
- **Would you use this again for the full cohort? Why/why not**: Can't
  assess yet -- needs the manual RStudio configuration step first.

## HIVIntact ("proviral")

- **Setup friction** (1-5, 1=trivial): 3 -- clean `pip install` into a
  dedicated venv, but Python 3.12's `venv` module no longer bundles
  `setuptools`, and the `proviral` CLI imports the legacy
  `pkg_resources` module (from `setuptools`) at startup. Compounded by
  `setuptools>=81` having fully dropped the `pkg_resources` shim, so
  even installing setuptools fresh doesn't fix it -- had to pin
  `setuptools<81` specifically. Both fixes are now folded into
  `setup_tools.sh` for future re-clones.
- **Documentation quality**: README documents the `--subtype` argument
  and its fixed list of supported subtypes, which is good, but doesn't
  flag its subtype-B-optimised heuristics as prominently as the tools
  review does. Also, `--help` is the only reliable source for its real
  flags -- the CLI takes `--working-folder` for output, not `-o` as one
  might guess from convention, and even `--working-folder` doesn't
  fully control where output lands: `intact.fasta`/`nonintact.fasta`/
  `orfs.json`/`errors.json` are written relative to the process's
  current working directory regardless, so `run_hivintact.sh` has to
  `cd` into the output dir before invoking it.
- **Failure modes encountered**: `ModuleNotFoundError: No module named
  'pkg_resources'` until the setuptools pin above.
- **Output usability**: Clean FASTA split into intact/non-intact, plus
  JSON detail on ORFs and classification errors per sequence -- easy to
  parse.
- **A finding worth flagging prominently**: classification is **not**
  robust to subtype choice on this cohort. Running the same 5-sequence
  subset with `--subtype B` (the comparison harness's original
  hardcoded default, before this session's fixes) called all 5 as
  intact; re-running with `--subtype A1` (the cohort's actual subtype
  per PRJNA207834's own metadata) called all 5 as non-intact -- a
  complete flip, not a marginal shift. This is exactly what the tools
  review's "subtype-B-optimised heuristics" caveat warned about, now
  demonstrated on real data rather than left as an abstract caveat.
- **Would you use this again for the full cohort? Why/why not**: Only
  with per-sample subtype confirmed first (e.g. from jpHMM's own call
  for each sample, rather than one fixed subtype for the whole cohort --
  this is A1/D/recombinant, and a recombinant sample run under either
  pure subtype's alignment is itself an approximation worth flagging).
  Given how much the A1-vs-B result diverged here, HIVIntact's raw
  intact/non-intact call should not be taken as final without
  cross-checking against HIVSeqinR once that's unblocked.

## Overall recommendation for this stage

**Poplars and HIVIntact both now produce real output** on this cohort:
Poplars' hypermutation screen ran cleanly (9,527 rows across the
subset), and HIVIntact, once its CLI flags and A1 subtype default were
corrected, classified all 5 subset sequences as non-intact -- a result
that should be treated as provisional given how much it changed between
subtype B and A1 (see above), and cross-checked against HIVSeqinR.
HIVSeqinR remains blocked on a one-time manual RStudio configuration
step that only the user can complete. None of the three tools' initial
failures reflected real problems running them at cohort scale; all were
fixable setup/environment/CLI-invocation issues, not data or design
problems -- but the subtype-sensitivity finding above is a real
scientific result, not an infrastructure issue, and belongs in the
methods write-up's limitations section.
