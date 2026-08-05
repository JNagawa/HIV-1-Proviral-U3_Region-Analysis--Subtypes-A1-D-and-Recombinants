# HIV-1 Proviral U3 Analysis — PacBio HIV-SMRTcap Pipeline Report
**Date:** 2026-08-05  ·  **Author:** Jovita Nagawa  ·  **Objective 1**

Comprehensive summary of the PacBio arm, one section per pipeline stage, each covering
purpose, methodology, results, challenges and conclusion. All figures are measured, not
estimated, and every stage was submitted as a Slurm job.

**Samples:** 124_4, 128_5, 203_3, 211_0 (Rakai cohort, SRA PRJNA1251704/1252688), plus
HXB2 / K03455.1 as reference and control.

**Status:** stages 1–6 complete for both assembly arms. Stage 7 (subtyping) was still
running at the time of writing and is marked pending; this report will be updated when it
lands.

---

## 0. Executive summary

| finding | consequence |
|---|---|
| These are **amplicon** libraries, not shotgun | de novo assembly structurally impossible for 3 of 4 samples |
| **Indel calling was silently off** (0 indels in all samples) | real deletions were being reported as unsequenced |
| **N-masking conflated "deleted" with "not sequenced"** | 106 interior bases mislabelled, 2 inside U3 |
| **U3 was taken from the wrong LTR** | 124_4 and 128_5 returned HXB2's own sequence and motifs |
| **HIV-Intact scored every sequence as subtype A1** | HXB2 itself was reported defective |
| **Poplars silently depended on equal-length input** | broke as soon as indels were applied |
| **124_4 carries a 21 bp insertion in its 3′ U3** | candidate regulatory finding, newly visible |

Three of four samples now yield genuine U3 sequence. 128_5 is excluded on evidence rather
than assumption. No stage now reports reference sequence as sample data.

---

## 1. QC and host depletion

**Purpose.** Remove host contamination and non-proviral sequence before assembly, and
establish which QC tools are usable on this input at all.

**Methodology.** Format-aware comparison of NanoFilt, chopper, fastp, seqkit and awk for
length filtering; Kraken2 for host removal; four deduplication strategies; plus a bespoke
N-strip step that extracts the proviral core and orients all reads forward.

**Results.**

| sample | reads in | proviral core kept | proviral fraction of read | reverse-complemented | after dedup |
|---|---|---|---|---|---|
| 124_4 | 41 | 41/41 | 52.7% | 25 | 27 |
| 128_5 | 75 | 75/75 | 0.2% | 32 | 60 |
| 203_3 | 19 | 19/19 | 78.3% | 10 | 12 |
| 211_0 | 30 | 30/30 | 43.9% | 21 | 25 |

Kraken2 retained 100% of reads in every sample — no host contamination detectable.

**Challenges.** The input is FASTA, not FASTQ, so NanoFilt, chopper, fastp and fastp
`--dedup` are all inapplicable; they are recorded as `not_applicable` with the reason
rather than as failures. 128_5 is striking: only **0.2%** of each read is proviral, the
rest host-derived flanking sequence — which is why the N-strip step is load-bearing rather
than cosmetic.

**Conclusion.** seqkit and awk are equivalent for length filtering; all four dedup
strategies agree exactly. Deduplication was later shown not to affect assembly outcome, so
it is not applied downstream.

---

## 2. Proviral extraction

**Purpose.** Recover the proviral portion of each read and discard host flanks.

**Methodology.** `minimap2` mask + N-strip, retaining the longest proviral core per read.

**Results.** 40/40, 75/75, 19/19 and 30/30 reads yielded a proviral core for
124_4, 128_5, 203_3 and 211_0 respectively. Exit 0, valid, ~1 s per sample.

**Challenges.** None outstanding at this stage.

**Conclusion.** Extraction is reliable and is not a limiting factor anywhere downstream.

---

## 3. Assembly

**Purpose.** Reconstruct each proviral genome, and quantify what the reference choice costs.

**Methodology.** Three arms over the same reads:

| arm | reference | role |
|---|---|---|
| `minimap2_consensus` | HXB2 (K03455.1) | fixed subtype-B baseline; the coordinate system |
| `minimap2_bestref` | closest LTR-complete Group M panel entry | cost of forcing A1/D reads onto a B backbone |
| `hifiasm` | none | reference-free comparator |

`minimap2 -x map-hifi` → `bcftools mpileup -X pacbio-ccs` → `bcftools call -m --ploidy 1` →
`bcftools consensus -m <mask>`, with unsequenced positions N-masked from read alignment
spans. Validity is **called (non-N) bases ≥ 8000**, not total length.

**Results.**

| sample | arm | length | called | %N | valid | SNV | ins | del |
|---|---|---|---|---|---|---|---|---|
| 124_4 | HXB2 | 9728 | 8914 | 8.36 | **1** | 941 | 10 | 14 |
| 124_4 | bestref (D:KU168271) | 9694 | 8829 | 8.92 | **1** | — | — | — |
| 128_5 | HXB2 | 9751 | 5469 | 43.91 | 0 | 538 | 11 | 6 |
| 128_5 | bestref (B:K03455) | 9751 | 5469 | 43.91 | 0 | — | — | — |
| 203_3 | HXB2 | 9730 | 1178 | 87.89 | 0 | 149 | 7 | 6 |
| 203_3 | bestref (A1:KU168256) | 9807 | 1134 | 88.43 | 0 | — | — | — |
| 211_0 | HXB2 | 9720 | 5976 | 38.51 | 0 | 598 | 5 | 6 |
| 211_0 | bestref (D:KU168271) | 9747 | 5958 | 38.87 | 0 | — | — | — |

hifiasm: 0 contigs for 124_4, 128_5, 203_3; 5977 bp for 211_0 only.

**The defining property of the data — amplicon structure:**

| sample | distinct read groups (HXB2 coordinates) |
|---|---|
| 124_4 | 814–8002 (×9); 1963–9719 (×7); 8986–9719 (×24) |
| 128_5 | 1594–7030 (×73) + 2 outliers |
| 203_3 | 27–1106 (×10); 27–1193 (×9) |
| 211_0 | 26–5120 (×21); 4410–6001 (×9) |

Reads stack into a few discrete intervals rather than tiling the genome, so an
overlap-layout-consensus assembler has no dovetail overlaps to work with. This single fact
explains the hifiasm result and every coverage-shaped anomaly downstream.

**Does subtype-matching help?**

| sample | vs HXB2 | vs best reference | verdict |
|---|---|---|---|
| 203_3 | 162 | **95** | helps (−41%) |
| 124_4 | 965 | 972 | a wash |
| 211_0 | 609 | 622 | a wash |
| 128_5 | 555 | 555 | selects HXB2 itself |

**Challenges.** Two defects were found here, both material:

1. **Indel calling was silently off.** `bcftools call -c` with generic mpileup defaults
   called **zero indels in all four samples**. At 203_3 position 168 all 19 reads carry a
   1 bp deletion and the caller emitted a plain reference call. Switching the caller
   `-c`→`-m` changed nothing; the fix was `mpileup -X pacbio-ccs` (`--indels-cns` plus
   HiFi-appropriate gap/homopolymer parameters), which recovered 24/17/13/11 indels.
   124_4's SNV count fell 1097 → 941 as indel-induced false substitutions resolved.

2. **The mask conflated "deleted" with "not sequenced."** `samtools depth` counts bases,
   so a position deleted in every read reports depth 0 — indistinguishable from a position
   no read reached. 106 interior bases were mislabelled, two of them inside U3 windows.
   Worse, `bcftools consensus` resolves a mask/variant collision **in favour of the mask**,
   so a called deletion was applied and then padded back to `N` (203_3: 9730 → 9743 bp).
   The fix derives coverage from read alignment spans; interior masking is now zero.

Substitutions were never at risk — mask and SNV sets are disjoint by construction, verified.

**Conclusion.** Both fixes were only coherent together: the caller fix alone leaves the `N`
in place *and* inflates length. Interior masking is now 0 in all samples and both arms, and
an `N` means "no read reached here" and nothing else. Every sample gained called bases.

---

## 4. Multiple sequence alignment

**Purpose.** Place all consensuses and HXB2 in one coordinate frame so HXB2's annotation
can be projected onto the samples.

**Methodology.** MAFFT, MUSCLE and Clustal Omega compared on the same input; MAFFT feeds
downstream.

**Results.**

| tool | arm | columns | gap % | integrity |
|---|---|---|---|---|
| mafft | HXB2 | 9855 | 1.3 | exact |
| muscle | HXB2 | 9830 | 1.0 | exact |
| clustalo | HXB2 | 9823 | 1.0 | exact |
| mafft | bestref | 9841 | 1.0 | exact |
| muscle | bestref | 9850 | 1.1 | exact |
| clustalo | bestref | 9854 | 1.1 | exact |

All three preserved every record's exact ungapped length — no sequence altered, only gaps
inserted. HXB2's row now carries 136 gaps (MAFFT), meaning samples genuinely have
insertions relative to it, so the coordinate liftover is being exercised for the first time.

**Challenges.** Clustal Omega **treats `N` as an alignable residue**. On the pre-fix data,
where every consensus was exactly reference-length and therefore provably collinear, MAFFT
and MUSCLE both returned the correct identity alignment (9719 columns, 0 gaps) while
Clustal Omega inserted **49 spurious indels per record**, sliding 203_3's 88% N-block
against real sequence. Since N-masking is now permanent and N fractions reach 88%, Clustal
Omega is unsuited to this data. That collinear dataset is retained as a ground-truth
aligner benchmark in `results/msa/pacbio/benchmark_collinear_20260804/`.

**Conclusion.** MAFFT is the correct choice and is what feeds downstream. Aligner
disagreement on gap placement (9823–9855 columns) is now real ambiguity rather than error,
since ground truth no longer exists.

---

## 5. U3 extraction

**Purpose.** Recover each sample's own U3 sequence — the region carrying the TFBS of
interest — without substituting reference sequence.

**Methodology.** HXB2's GenBank record supplies both LTR and both R-region boundaries;
U3 = `[LTR_start, R_start)`. Those genome coordinates are lifted to alignment columns by
walking HXB2's own row, then every record is sliced at the same columns and degapped.
Each record uses whichever LTR copy carries more sequenced base, compared as a **fraction**
of each window (HXB2's copies are 453 and 454 nt, so raw counts would hand every
fully-covered record to the 3′ copy on a one-base technicality).

**Results.**

| record | LTR used | called 5′ | called 3′ | extracted U3 |
|---|---|---|---|---|
| K03455.1 | 5′ | 453 | 454 | 453 bp |
| 124_4 | **3′** | 0 | **474** | 474 bp, 147 diffs from HXB2 |
| 128_5 | 5′ (tie) | 0 | 0 | **not recovered — excluded** |
| 203_3 | 5′ | 428 | 0 | 454 bp |
| 211_0 | 5′ | 428 | 0 | 453 bp |

**Two artifacts eliminated.** Previously 124_4 and 128_5 both returned U3 **byte-identical
to HXB2**, so their motif hits were HXB2's. 124_4 is now a genuine third data point (474 bp,
147 differences), and 128_5 is auto-flagged as unrecoverable rather than silently reported.

**A candidate finding.** 124_4 carries a **21 bp insertion in its 3′ U3** — real sequence,
zero N, at roughly HXB2 position 9400. It is reproduced against **both** references and by
**all three aligners**. Placement within the run differs between aligners
(`tttataagaactgaactgctg` / `ACTGAACTGCTGACACCAGAG` / `ACTGCTGACACCAGAGACTGC`), which is the
signature of a **tandem duplication** in a repetitive stretch; the inserted sequence
resembles adjacent HXB2 U3 containing `ACTGCTGACA`.

**Challenges.** U3 position genuinely varies between subtypes, so this needed verifying
rather than assuming. Literature (PubMed): Mbondji-Wonje 2018 found the R region "very well
conserved" across clades A1/B/C/D/F2 while U3 inter-strain dissimilarity reaches **25%**;
Jeeninga 2000 found "a unique LTR enhancer-promoter configuration for each subtype" with
NF-κB counts from one (E) to three (C); Naghavi 1999 showed subtype C's third site arises
from an **insertion**; Parreira 2006 found only **63.3%** of subtype C isolates carry three
NF-κB sites. Because the variation is partly insertional and not even fully
subtype-determined, no per-subtype coordinate system exists or would be correct — LANL
publishes HXB2-only coordinates. Alignment-anchored, per-sequence derivation is therefore
the published method, not a workaround.

**Known limitation.** The window ends at the column of HXB2's last U3 base, so an insertion
sitting exactly at the U3/R junction would be clipped. Checked: no junction insertions
occur in this dataset, so the limitation is real but currently inert.

**Conclusion.** Our measured U3 identities (203_3 85.0%, 211_0 86.5% — i.e. 15.0% and 13.5%
dissimilarity) sit comfortably inside the published inter-strain range, so these look like
genuine HIV-1 U3 diversity rather than alignment artifacts.

---

## 6. Motif mapping

**Purpose.** Identify TFBS in each sample's U3 and compare across samples.

**Methodology.** FIMO, TFBSTools and MOODS for TFBS; gquad and pqsfinder for
G-quadruplexes; each run against a positive control as well as the sample set.

**Results (FIMO).**

| sequence | hits |
|---|---|
| K03455.1 (HXB2) | RELA 350–359, RELA 364–373, SP1 388–396 |
| **124_4** | RELA 371–380, RELA 384–393 |
| **128_5** | *absent — U3 not recovered* |
| 203_3 | RELA 364–373, SP1 375–383, SP1 389–397 |
| 211_0 | NFKB1 65–77 (×2), RELA 350–359, RELA 364–373, SP1 389–397 |

Tool hit counts: FIMO 14, MOODS 10, TFBSTools 89, gquad 15, pqsfinder 4.

Compare against the pre-fix table, where 124_4 and 128_5 both showed HXB2's exact three
hits. HXB2's own hits are unchanged across the re-run, which serves as the control.

**Observations.**

1. **124_4's two RELA sites are displaced by ~+21** (350→371, 364→384), consistent with the
   21 bp insertion lying upstream of them. These are very likely the same two canonical
   NF-κB sites, shifted — not novel sites.
2. **124_4 shows no SP1 hit.** The expected position after the shift (~409–417) exists in
   its 474 bp U3, so the site appears genuinely absent or degraded. Worth following up.
3. **203_3 now shows two SP1 sites** where it previously showed one.
4. **203_3 still lacks RELA 350–359**, the finding flagged earlier as possibly real. It
   survives, and is no longer obscured by two artifact samples that "had" the site.

**Challenges.** Positions are currently reported in per-sample U3-relative coordinates. Once
U3 lengths differ — which they now do — those positions are **not comparable across
samples**, exactly as 124_4's +21 shift demonstrates. LANL's Reference Sequence Coordinate
Search publishes curated HXB2 coordinates for the named NF-κB-I/II and Sp1-I/II/III sites,
so hits should be reported as **named sites** with numeric positions as supporting detail.
This is not yet implemented.

**Conclusion.** Three of four samples yield genuine motif profiles. The headline biological
question — whether 124_4's insertion creates or destroys a binding site — is now answerable
and points at the missing SP1 site.

---

## 7. Biological filtering

**Purpose.** Classify each provirus as intact or defective, and record what the verdict
actually rests on.

**Methodology.** HIV-Intact (per-sample, subtype-matched), Poplars/Hypermut 3 (APOBEC
G→A screen on the aligned FASTA), and HIVSeqinR (attempted). Per-gene read coverage is
computed from the mask BED and HXB2's CDS annotation so each gene is marked assessable or
not.

**Results.**

| tool | outcome | runtime |
|---|---|---|
| HIV-Intact | **1 intact, 4 non-intact** | ~390–435 s |
| Poplars | 5135 rows (HXB2 arm), 5102 (bestref) | ~21–30 s |
| HIVSeqinR | failed — structurally incompatible input | — |

The intact sequence is **HXB2 itself**, i.e. the positive control passes.

**Subtype used per sequence (bestref arm):**

| sequence | subtype | source |
|---|---|---|
| 124_4 | D | assembly reference D:KU168271 |
| 128_5 | B | assembly reference B:K03455 |
| 203_3 | A1 | assembly reference A1:KU168256 |
| 211_0 | D | assembly reference D:KU168271 |
| K03455.1 | HXB2 | reference control |

**Per-gene read coverage — the basis for any verdict:**

| gene | CDS bp | 124_4 | 128_5 | 203_3 | 211_0 |
|---|---|---|---|---|---|
| gag | 1503 | 98% | 47% | 27% | **100%** |
| pol | 2739 | **100%** | **100%** | 0% | **100%** |
| vif | 579 | **100%** | **100%** | 0% | **100%** |
| vpr | 237 | **100%** | **100%** | 0% | **100%** |
| tat | 261 | **100%** | 82% | 0% | 66% |
| rev | 351 | **100%** | 22% | 0% | 9% |
| env | 2571 | **100%** | 31% | 0% | 0% |
| nef | 372 | **100%** | 0% | 0% | 0% |

| sample | genes assessable | verdict support |
|---|---|---|
| 124_4 | 7–8 of 8 | whole-genome verdict supported |
| 211_0 | 4 of 8 | partial |
| 128_5 | 3 of 8 | partial |
| 203_3 | **0 of 8** | **verdict unsupported — no complete gene sequenced** |

**Challenges.** Three defects found, two fixed:

1. **HIV-Intact scored every sequence as subtype A1**, including HXB2 — which is subtype B.
   HXB2 was consequently reported defective with `FrameshiftInOrf`,
   `MajorSpliceDonorSiteMutated` and two `MisplacedORF`s: the signature of a misconfigured
   run, not a finding. Fixed by scoring each sequence separately against its own subtype.
   HXB2 now returns intact.

2. **Poplars silently depended on equal-length input.** Its `hypermut` asserts every record
   is the same length; the raw consensuses satisfied that only by accident while zero indels
   were being called. Once indels were applied, lengths diverged (9719–9751 bp) and the
   assertion fired. Fixed by feeding the MAFFT alignment, which is the correct input for a
   positional comparison anyway.

3. **HIVSeqinR could not be made to work on this data.** All mechanical blockers were
   genuinely fixed — the release ZIP was never unpacked (so the R script did not exist), the
   R `muscle` and `pwalign` dependencies were absent, nine `pairwiseAlignment()` and
   `pattern()` calls are defunct in current Biostrings, `MyBlastnDir` pointed at the
   author's own machine, and the wrapper's `N→A` substitution was fabricating up to 8552
   bases per sample while mimicking the very APOBEC G→A signature the tool detects (that
   was replaced with splitting on N runs, so nothing is invented). What remains is
   structural: HIVSeqinR requires *"linear HIV genomes WITH FLANKING PRIMER binding sites at
   5′ and 3′ ends"* from a specific de novo pipeline, and its 5′Psi+gag filter rejects
   amplicon segments outright. The SMRTcap 2nd-round PCR primers are also undocumented, so
   autotrim cannot work. The summary now records the actual R error rather than a blank.

**Conclusion.** HIV-Intact plus Poplars is the supportable combination for this data.
Crucially, 203_3 cannot receive an intactness verdict from **any** tool — it has no complete
gene — and that is a data limitation, not a software one. Reporting per-gene assessability
alongside each verdict is what keeps that honest.

---

## 8. Subtyping — PENDING

**Purpose.** Assign a formal subtype per sample, independent of the closest-reference match
used in assembly.

**Methodology.** jpHMM per sample, plus IQ-TREE2 on the whole-subset MAFFT alignment.

**Status.** Running at time of writing (~55 min elapsed, on sample 2 of 4 in both arms).
jpHMM is by far the slowest tool in the pipeline.

**Expectation to test.** The assembly stage's closest-reference matches suggest three
subtype-D and one subtype-A1 sample, consistent with a Rakai cohort. These are
closest-reference matches, **not** subtype calls; jpHMM plus IQ-TREE is the formal call.
An earlier concern — that subtyping was circular because jpHMM was reading an HXB2-padded
consensus — is now addressed, since the consensus is N-masked rather than reference-filled.

**This section will be completed when the jobs finish.**

---

## 9. Cross-cutting conclusions

1. **The pipeline no longer reports reference sequence as sample data.** This was the single
   largest class of error, and it appeared in three independent places: HXB2 padding in the
   consensus, the wrong LTR in U3 extraction, and N→A fabrication in HIVSeqinR staging.
2. **"Not sequenced" is now distinguishable from "absent" at every stage** — via the mask,
   the LTR-choice report, the per-gene coverage table, and auto-flagging of 128_5.
3. **Amplicon structure is the root constraint.** It rules out de novo assembly, limits
   which genes can be assessed, and determines which LTR each sample can contribute.
4. **Positive controls earn their place.** HXB2 returning "defective" was what exposed the
   HIV-Intact subtype bug; HXB2 returning "intact" is now the check that the stage works.
5. **A tool completing successfully is not evidence it was configured correctly.** Poplars,
   HIV-Intact and the variant caller all ran to exit 0 while producing wrong or empty
   results.

## 10. Open items

1. Report motif hits as **named sites** using LANL's curated HXB2 coordinates, so
   cross-sample positions are comparable (highest value remaining).
2. Investigate **124_4's missing SP1 site** and characterise the 21 bp insertion as a
   probable tandem duplication.
3. Complete **subtyping** and reconcile jpHMM/IQ-TREE calls with the closest-reference matches.
4. Anchor the U3 window end on HXB2's R-start column so junction insertions are not clipped.
5. Report per-sample U3 length and indel structure so subtype coordinate variation is
   measured rather than assumed.
6. Vendor LANL's `HXB2.xlsx` into `data/reference/` for offline, version-pinned site coordinates.
7. Optional: iterative consensus refinement; reference-free amplicon consensus (`abpoa`/`spoa`).

---

## 11. Reproduction

| item | path |
|---|---|
| Stage summaries | `results/<stage>/pacbio/<arm>/summary.tsv` |
| Assembly findings report (auto) | `results/assembly/pacbio/assembly_report.md` |
| Intactness basis report (auto) | `results/biological_filtering/pacbio/<arm>/intactness_basis_report.md` |
| Aligner ground-truth benchmark | `results/msa/pacbio/benchmark_collinear_20260804/` |
| Slide deck | `writeups/PacBio_Pipeline_Report.pptx` |
| Chain launcher | `scripts/utils/run_pacbio_chain.sh` |

Every stage takes `ASSEMBLY_ARM` and writes to `results/<stage>/pacbio/<arm>/`, so the two
reference arms never collide. Run on compute nodes only:

```bash
sbatch --export=ALL,ASSEMBLY_ARM=minimap2_consensus \
       scripts/utils/run_comparison_step.slurm.sh scripts/<stage>/pacbio/compare_<stage>_pacbio.sh
```

**Note:** `peak_rss_mb` is blank throughout — GNU `time` is unavailable on this cluster, so
runtime is measured and memory is not.
