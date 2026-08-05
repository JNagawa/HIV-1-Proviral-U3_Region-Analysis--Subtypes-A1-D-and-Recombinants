# HIV-1 Proviral U3 Analysis — Indel Calling, Mask Semantics, and the U3 Coordinate Frame
**Date:** 2026-08-05  ·  **Author:** Jovita Nagawa  ·  **Objective 1** methodology work

Follow-on to `status_report_04082026.md`, which resolved the PacBio assembly step and
left six open items. This report covers two defects found while working item 1 (U3
extraction takes the wrong LTR) and item 4 (regenerate motif results), the automated
findings report added as a result, and a literature review establishing whether U3
coordinates vary by subtype — which determines whether the extraction design is sound.

All measurements were made on 2026-08-05 on compute nodes via `srun`, against the
existing BAMs in `results/assembly/pacbio/minimap2_consensus_out/`. Nothing here is
estimated.

---

## 1. Executive summary

- **Indel calling was silently off.** `bcftools call -c` with mpileup's generic defaults
  called **zero indels in all four samples**, even where every read agreed. At 203_3
  `K03455.1:168` all 19 reads carry a 1 bp deletion and the caller emitted a plain
  reference call.
- **The N-mask was labelling confident deletions as unsequenced.** `samtools depth`
  counts *bases*, so a position deleted in every read reports depth 0 and was
  indistinguishable from a position no read reached. 106 interior bases across the four
  samples were masked this way, including two inside the U3 window.
- **`bcftools consensus` resolves a mask/variant collision in favour of the mask.**
  Applying the deletion and then masking the same position puts an `N` back and inflates
  the sequence. So the two defects had to be fixed together; fixing only the caller makes
  the output worse.
- **Masking cannot damage a substitution** — mask and variant sets are disjoint for SNVs.
  The only failure mode was deletions, which is now closed.
- **A findings report is now generated automatically** after every assembly run, carrying
  the unsequenced-vs-deleted distinction that `summary.tsv` cannot express.
- **U3 coordinates do vary by subtype, but no per-subtype coordinate system exists or
  should be adopted.** The variation is not even fully subtype-determined. The
  alignment-anchored liftover already implemented is the published method.

---

## 2. Indel calling was disabled in practice

### 2.1 The diagnostic case

203_3 at `K03455.1:167-168`, read directly off the pileup:

```
pos 167  G  19  .-1T .-1T .-1T ... (all 19 reads)
pos 168  T  19  *******************
```

All 19 reads delete the T at 168. `samtools depth` reports **0** there, because deleted
positions carry no base to count. The VCF nonetheless contained a record —
`DP=19;DP4=19,0,0,0` — i.e. the position was piled up, 19 reads supported it, and no
variant was called.

### 2.2 Configuration comparison

Three configurations over the existing BAMs:

| sample | A: `call -c` (was in use) | B: `call -m` | C: `call -m` + `mpileup -X pacbio-ccs` |
|---|---|---|---|
| 124_4 | 1097 SNV / **0** indel | 1097 SNV / **0** indel | 941 SNV / **24** indel |
| 128_5 | 538 SNV / **0** indel | 538 SNV / **0** indel | 538 SNV / **17** indel |
| 203_3 | 149 SNV / **0** indel | 149 SNV / **0** indel | 149 SNV / **13** indel |
| 211_0 | 598 SNV / **0** indel | 598 SNV / **0** indel | 598 SNV / **11** indel |

**Switching the caller was not the fix.** `-c` → `-m` in isolation changed nothing. The
blocker was in **mpileup's** indel-candidate detection defaults. The `pacbio-ccs` profile
(`--indels-cns` plus HiFi-appropriate gap, homopolymer and indel-bias parameters)
recovers them, and calls the diagnostic locus correctly:

```
203_3  POS=167  REF=GTT  ALT=GT   QUAL=80.4     # the 19/19 deletion
211_0  POS=398  REF=TGGGG ALT=TGGG QUAL=49.4    # a 1bp deletion in a GGGG homopolymer
```

`-m` is kept regardless: `-c` is deprecated and `--indels-cns` is designed against the
multiallelic caller.

### 2.3 A conclusion in the previous report is affected

124_4's SNV count drops **1097 → 941**. Those 156 "substitutions" were misaligned-indel
artifacts resolving into their correct events. Report `04082026` §8.2 used raw variant
counts to judge whether subtype-matching helps (124_4: 1097 vs 1109, called "a wash").
Those counts are now known to be inflated, so **§8.2 must be recomputed** once both arms
re-run.

---

## 3. The mask was conflating "deleted" with "unsequenced"

### 3.1 Scale of the problem

Interior masked intervals — those with sequenced bases on both sides — under the old
depth-based mask:

| sample | mask intervals | interior intervals | interior bp |
|---|---|---|---|
| 124_4 | 14 | 13 | 79 |
| 128_5 | 6 | 4 | 7 |
| 203_3 | 8 | 6 | 13 |
| 211_0 | 7 | 5 | 7 |

All 106 interior bases had a VCF record and **none carried a real ALT** — every one was a
covered position, not missing data. Two sit inside the 5' U3 window: 203_3 at 168 and
211_0 at 399. Both are homopolymer deletions (HXB2 has `TT` at 168-169 and `GGGG` at
399-402), which is why the called indel anchors one base away from where the mask
flagged it.

### 3.2 Which wins — mask or variant

Tested directly on 203_3, positions 160-180:

```
HXB2 reference               CAGAGAAGTTAGAAGAAGCCA
variants applied, no mask    CAGATGAAG T AGAGAAGGCTA     <- deletion applied
variants applied + mask      CAGATGAAG N TAGAGAAGGCT     <- deletion undone, N in its place
```

| consensus | length | N |
|---|---|---|
| HXB2 reference | 9719 | 0 |
| variants applied, no mask | 9730 | 0 |
| variants applied **+ mask** | **9743** | 8565 |

9743 − 9730 = **+13**, exactly 203_3's interior masked base count. Every deleted base is
reinstated as `N`: `bcftools consensus` applies the deletion, then the mask writes a
character back at that reference position. **The mask wins.**

Consequence: the caller fix alone would have left the `N` in place *and* added spurious
length. The two fixes are only coherent together.

### 3.3 Substitutions were never at risk

For SNVs the mask and variant sets are disjoint by construction — the mask holds only
zero-coverage positions and calling a variant requires coverage. Verified: **0 called
SNVs** fall inside any masked interval, and 203_3's `POS=173 A>G` / `POS=174 G>A` both
appear in the consensus. A sample base differing from HXB2 by substitution was never
converted to `N`. The concern applied to deletions only.

### 3.4 The fix — mask on read spans, not pileup depth

Coverage is now derived from read alignment spans: `POS .. POS + (reference-consuming
CIGAR ops) - 1`, where `M/D/N/=/X` advance the reference and `I/S/H/P` do not, so
soft-clips and insertions correctly fail to extend coverage. A position inside a read's
span is sequenced by definition.

| sample | old: intervals (interior) | new: intervals (interior) | called bases: old → new |
|---|---|---|---|
| 124_4 | 14 (13) | **1 (0)** | 8826 → **8915** |
| 128_5 | 6 (4) | **2 (0)** | 5430 → **5469** |
| 203_3 | 8 (6) | **2 (0)** | 1153 → **1178** |
| 211_0 | 7 (5) | **2 (0)** | 5969 → **5977** |

Interior masking is eliminated entirely. Every remaining mask interval is an amplicon-edge
gap, so an `N` now means "no read reached here" and nothing else.

Rebuilt consensus with both fixes:

| sample | length | N | %N | called |
|---|---|---|---|---|
| 124_4 | 9728 | 813 | 8.36% | 8915 |
| 128_5 | 9751 | 4282 | 43.91% | 5469 |
| 203_3 | 9730 | 8552 | 87.89% | 1178 |
| 211_0 | 9720 | 3743 | 38.51% | 5977 |

Consensus length is no longer fixed at 9719, because indels are now applied. This is
correct, and the downstream design tolerates it — the U3 extraction lifts HXB2 coordinates
through the MAFFT alignment rather than assuming fixed offsets.

---

## 4. Automated findings report

`summary.tsv` records whether a run succeeded. It cannot express the distinction the
motif work depends on: for a given region, was it **sequenced**, or is it **absent**?
Conflating those is how HXB2's own motif hits were once reported as sample results.

`scripts/utils/report_pacbio_assembly.sh` now runs at the end of every assembly job
(non-fatally — a failed report must not discard a completed assembly) and writes
`results/assembly/pacbio/assembly_report.md`. Per arm it reports coverage, mask
composition, SNV/insertion/deletion counts, consensus length and called bases, merged
read spans (the amplicon structure), a cross-arm divergence comparison, and:

- a **masking audit** that fails loudly if interior masking ever reappears, so the §3
  defect cannot regress silently;
- a **U3 coverage** section per LTR copy, with an explicit "usable copy" verdict.

Figures are taken from the mask BED and VCF in **reference coordinates**, not from the
consensus FASTA, whose coordinates now drift as indels are applied.

Validated by reproducing figures derived independently in report `04082026`: 124_4 at
90.8% of the genome spanned, U3 coverage 0% (5') / 100% (3'), 203_3 and 211_0 at 94% of
the 5' copy, 128_5 at 0% / 0%. Run against the stale outputs it correctly flags them as
defective.

Two bugs fixed during its own verification: reference length was read off the mask BED,
which is only correct when a sample's 3' end happens to be uncovered (124_4 reported
90.1% instead of 90.8%); and a missing `read_spans.tsv` produced a silently empty table
reading as "no amplicon structure" rather than "not measured".

---

## 5. Does U3 vary in position between subtypes?

This determines whether anchoring U3 on HXB2's annotation is defensible.

### 5.1 The question was invisible in our data until now

Every record in the current MAFFT alignment is **exactly 9719 columns with zero gaps** —
a perfect ungapped block. That is only possible because the old consensuses were all
exactly reference-length, which follows directly from §2's zero indels. **The coordinate
liftover has never been exercised on length-variant input.** After the fixes, consensus
lengths become 9728 / 9751 / 9730 / 9720, MAFFT will insert gap columns, and the liftover
starts to carry real weight.

### 5.2 Literature (retrieved via PubMed)

U3 architecture varies by subtype through **insertions and duplications**, not only
substitutions:

| study | finding |
|---|---|
| Jeeninga et al. 2000, *J Virol* | LTRs of subtypes A–G cloned; "a unique LTR enhancer-promoter configuration for each subtype"; NF-κB site count varies from **one (subtype E) to three (subtype C)** |
| De Baar et al. 2000, *AIDS Res Hum Retroviruses* | 60 viruses incl. **20 subtype A and 10 subtype D**; describes "subtype-specific differences in sequences encompassing the regulatory elements of the LTR" |
| Naghavi et al. 1999, *AIDS Res Hum Retroviruses* | an **insertion** creating a potential third NF-κB site in subtype C; USF site in the NRE also subtype-specific |
| Mbondji-Wonje et al. 2018, *PLoS One* | U3R across A1, B, C, D, F2 + CRFs. **R region very well conserved**; U3 inter-strain dissimilarity **up to 25%**. TATAA in most strains but TAAAA in CRF01_AE/CRF22_01A1; NF-κB counts vary by strain |
| van Opijnen et al. 2004, *J Virol* | database-wide subtype-specific conservation of certain TFBS |
| Parreira et al. 2006, *Microbes Infect* | Mozambican subtype C: only **63.3%** of viruses carried three NF-κB sites |
| Engin et al. 2026, *bioRxiv* (preprint) | MPRA functional atlas across clades; LTR activity differs **among proviruses from the same individual** |

DOIs: [Jeeninga](https://doi.org/10.1128/jvi.74.8.3740-3751.2000) ·
[De Baar](https://doi.org/10.1089/088922200309160) ·
[Naghavi 1999](https://doi.org/10.1089/088922299310197) ·
[Mbondji-Wonje](https://doi.org/10.1371/journal.pone.0195661) ·
[van Opijnen](https://doi.org/10.1128/jvi.78.7.3675-3683.2004) ·
[Parreira](https://doi.org/10.1016/j.micinf.2006.05.005) ·
[Engin](https://doi.org/10.64898/2026.04.03.716403)

### 5.3 Why the current design is the right one

Three points follow, and together they justify the implementation rather than requiring a
change to it:

1. **R is conserved while U3 varies.** Mbondji-Wonje found the R region "very well
   conserved" across A1/B/C/D/F2 with U3 dissimilarity up to 25%. Defining U3 as
   `[LTR_start, R_start)` and projecting HXB2's R boundary is therefore sound, while the
   U3 *interior* demands alignment-based liftover.
2. **Coordinate slicing would be wrong.** Subtype variation is partly insertional, so a
   raw genome-coordinate slice (`seqkit subseq -r`, the approach abandoned earlier) cannot
   be correct. Lifting HXB2 coordinates through alignment columns follows insertions.
3. **No per-subtype coordinate system exists, and adopting one would be wrong anyway.**
   Parreira's 63.3% figure shows element count varies *within* subtype — a "subtype C
   coordinate table" would be wrong for over a third of subtype C isolates. This is why
   the field derives elements per-sequence from an alignment. **The per-sample
   alignment-anchored liftover is the published method, not a workaround.**

Our own numbers are consistent with the literature: 203_3 and 211_0 show 85.0% and 86.5%
identity to HXB2 U3, i.e. 15.0% and 13.5% dissimilarity, comfortably inside
Mbondji-Wonje's reported inter-strain range of up to 25%. The genuine U3 sequences look
like real HIV-1 U3 diversity rather than alignment artifacts.

### 5.4 Los Alamos resources

LANL's coordinate data is **HXB2-only** (plus Mac239 for SIV/HIV-2); there is no
per-subtype coordinate table. What is available:

| resource | provides |
|---|---|
| Reference Sequence Coordinate Search | HXB2 coordinates for 5'/3' LTR split into **U3/R/U5**, plus named **TATA box, NF-κB-I/II, Sp1-I/II/III**, TAR |
| `HXB2.xlsx` (In-depth Genome Annotation) | downloadable base-by-base landmark table incl. regulatory regions |
| Sequence Locator | maps a query sequence onto HXB2 coordinates — an independent cross-check on our own liftover |
| Subtype Reference Alignments; Consensus/Ancestral Sequences | ~4 representatives per subtype/CRF; a consensus per subtype — a route to subtype-appropriate U3 baselines |

Not yet confirmed: whether **LTR** is among the selectable pre-defined regions for
alignment download. If it is, that is a better source of LTR-complete diversity than
filtering the gag→nef compendium, which in report `04082026` §7.1 cost us down to one D
and one A1.

### 5.5 Consequent change to motif reporting

Report `04082026` §5.3 compares hits across samples by U3-relative position ("RELA
350–359"). Once U3 lengths differ between samples, those positions are **not
comparable** and the table silently assumes they are.

Recommended fix, now that LANL supplies named site coordinates: **report each hit as a
named site** (NF-κB-I/II, Sp1-I/II/III) by mapping onto LANL's curated HXB2 coordinates,
keeping numeric positions as supporting detail only. A named site survives length
variation in a way a raw offset does not, and it is directly interpretable. FIMO must
still run on ungapped sequence, so hit positions map back through the gapped U3 the
extraction already writes (`--out-gapped`).

---

## 6. MSA: Clustal Omega mis-aligns N-masked input

### 6.1 State of the two arms

`minimap2_consensus` completed for all three aligners; `minimap2_bestref` was interrupted
mid-MAFFT (0-byte output, log truncated, no `.time`, no `summary.tsv`) and MUSCLE and
Clustal Omega never started, since the aligners run in sequence.

| tool | exit | valid | records | columns | gaps/record | runtime |
|---|---|---|---|---|---|---|
| mafft | 0 | 1 | 5/5 | 9719 | **0** | 112.2 s |
| muscle | 0 | 1 | 5/5 | 9719 | **0** | 159.9 s |
| clustalo | 0 | 1 | 5/5 | 9768 | **49** | 204.8 s |

### 6.2 Clustal Omega's output is wrong, and there is ground truth to prove it

Every input sequence is exactly 9719 bp, produced by `bcftools consensus` with zero indels
called (§2) — substitutions and N-masking only. The inputs are therefore **strictly
collinear with HXB2**, and the single correct alignment is the identity: 9719 columns,
zero gaps. MAFFT and MUSCLE both return exactly that. Clustal Omega returns 9768 columns
with 49 gaps in every record and **no all-gap column**, so these are not trailing padding
but 49 spurious indels per sequence.

Residues at the offending columns show the mechanism:

```
column:      8940    9037    9142    9303    1531    2389    4133
K03455.1        -       -       -       -       A       A       A
124_4           -       -       -       -       A       A       A
128_5           -       -       -       -       N       A       A
203_3           N       N       N       N       -       -       -
211_0           -       -       -       -       A       A       A
```

**Clustal Omega treats `N` as an alignable residue.** At 8940-9303 it gave 203_3 a
single-`N` "insertion" and gapped every other record; at 1531-4133 it gapped 203_3 against
real bases. It slid 203_3's N-run (that sample is 88% N) against real sequence. The same
pattern holds for 128_5's extra gaps at 7244-7756 and 211_0's at 6209-7300 — each begins
exactly where that sample's read coverage ends and its N-block starts.

### 6.3 Why it matters downstream

The U3 extraction anchors on HXB2's row and slices **all records at the same columns**. In
Clustal Omega's alignment 203_3 is out of register with HXB2 by up to 49 positions, so its
U3 slice would be silently wrong. Nothing is currently broken, because only MAFFT feeds
the subtyping and motif stages (`ALIGNMENT=.../mafft_aligned.fasta`). But since N-masking
is the correct and now permanent representation of unsequenced regions, and N fractions
here reach 88%, **Clustal Omega is systematically unsuited to this data** — this is a
result of the tool comparison, not a run failure.

### 6.4 A positive control worth preserving

Because the pre-fix consensuses are provably collinear, this dataset is an accidental
**ground-truth benchmark for aligner accuracy**: MAFFT and MUSCLE score perfectly, Clustal
Omega introduces 49 false indels per sequence. That control disappears once the §2/§3
fixes land and consensus lengths genuinely differ, so the current
`results/msa/pacbio/minimap2_consensus/` alignments should be retained as a benchmark
snapshot rather than overwritten in place.

---

## 7. Open items

Carried forward from `04082026` §9, with status:

1. ~~U3 extraction takes the wrong LTR~~ — **implemented, not yet executed.** The
   extraction now picks whichever LTR copy carries more sequenced base, compared as a
   fraction of each window (HXB2's copies are 453 and 454 nt, so raw counts would favour
   the 3' copy on a one-base technicality). Never run — no U3 output exists newer than
   2026-07-29.
2. **128_5 should be excluded from motif analysis** — confirmed by the new report: 0% U3
   coverage at both LTRs.
3. **Subtyping is circular** — still open; jpHMM must re-run against the N-masked
   consensus.
4. **Motif results need regenerating** — per-arm fan-out implemented in
   `run_pacbio_chain.sh`; MSA done for `minimap2_consensus` (all three aligners),
   **`minimap2_bestref` interrupted mid-MAFFT** and must be redone. No downstream stage
   has run per-arm yet.
5. **Optional:** iterative consensus refinement — not implemented.
6. **Optional:** reference-free amplicon consensus (`abpoa`/`spoa`) — not implemented.

New items from this report:

7. **Recompute `04082026` §8.2** — its variant-count comparison used indel-inflated SNV
   counts.
8. **Junction insertions are clipped by the U3 liftover.** The window ends at the column
   of HXB2's last U3 base, so a sample insertion between HXB2's U3 end and R start falls
   outside. Anchor the end on HXB2's R-start column minus one instead.
9. **Report per-sample U3 length and indel structure** so subtype coordinate variation is
   measured rather than assumed.
10. **Adopt named-site motif reporting** (§5.5) before regenerating the motif comparison.
11. **Vendor `HXB2.xlsx`** into `data/reference/` for offline, version-pinned site
    coordinates, consistent with the vendored compendium alignment.

---

## 8. Files changed

| item | path |
|---|---|
| Caller + mask fixes | `scripts/assembly/pacbio/compare_assembly_pacbio.sh` |
| Findings report generator (new) | `scripts/utils/report_pacbio_assembly.sh` |
| Generated report | `results/assembly/pacbio/assembly_report.md` |
| Per-read reference spans (new artifact) | `<arm>_out/<sample>.read_spans.tsv` |
| U3 extraction (LTR choice, uncommitted) | `scripts/utils/extract_u3_by_hxb2_anchor.sh` |
| Per-arm chain fan-out (uncommitted) | `scripts/utils/run_pacbio_chain.sh` |

Re-run (never on the login node):

```bash
THREADS=4 srun -p shared -c 4 --mem=8G -t 30 scripts/assembly/pacbio/compare_assembly_pacbio.sh
```

Note the stale `<sample>.covered.pos` files from the old depth-based mask are no longer
written and can be deleted.
