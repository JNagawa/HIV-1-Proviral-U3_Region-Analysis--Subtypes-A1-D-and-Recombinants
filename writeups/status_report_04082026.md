# HIV-1 Proviral U3 Analysis — PacBio Assembly Method Resolution
**Date:** 2026-08-04  ·  **Author:** Jovita Nagawa  ·  **Objective 1** methodology work

This report resolves the PacBio assembly step. It documents why hifiasm produced no
output, a previously undetected reference-contamination defect in the reference-guided
consensus, the resulting three-arm assembly design, and two implementation bugs found
and fixed during verification. All measurements were made on 2026-08-04 on compute
nodes via `srun`; nothing in this report is estimated or assumed.

---

## 1. Executive summary

- **Root cause of the hifiasm failure is the data, not the parameters.** The PacBio
  SMRTcap libraries are **amplicon**, not shotgun: reads stack into 1–3 discrete
  intervals rather than tiling the genome. Overlap-layout-consensus assembly requires
  dovetail overlaps, so for 2 of 4 samples an overlap graph cannot exist at any
  parameter setting. Confirmed by a 7-config × 4-sample sweep plus `--hom-cov`,
  dedup and contained-read-removal tests.
- **The reference-guided consensus was silently substituting HXB2 for unsequenced
  sample sequence.** `bcftools consensus` leaves uncovered positions as the reference
  base and emits no `N`, so every sample reported a confident `9719bp, 0.00% N`. In
  reality 203_3's consensus was 88% HXB2, 128_5's 44%, 211_0's 39%.
- **Two of four samples' U3 motif results are reference artifacts.** 124_4 and 128_5
  produce U3 sequences **byte-identical to HXB2**, so their FIMO hits are HXB2's hits.
- **Assembly step now has three arms** — HXB2 baseline, subtype-matched reference, and
  de novo — all running to completion with honest metrics.
- **Two implementation bugs** were caught during verification: an LTR-less reference
  panel that discarded U3 reads, and a stale `.fai` that produced a variant-free
  consensus while appearing to succeed.

---

## 2. The core finding — these are amplicon libraries

Mapping the extracted proviral reads to HXB2 and tabulating distinct alignment
intervals:

| sample | reads | depth | distinct spans on HXB2 (read count) | genome covered |
|---|---|---|---|---|
| 124_4 | 40 | 14.1× | 814–8002 (×9), 1963–9719 (×6), 8986–9719 (×24) | 90.8% |
| 128_5 | 75 | 41.8× | **1594–7030 (×74)** — single amplicon | 55.9% |
| 203_3 | 19 | 2.2× | **27–1106 / 27–1193** — single amplicon | 11.9% |
| 211_0 | 30 | 12.5× | 26–5120 (×21), 4410–6001 (×9) | 61.4% |

This single fact explains every downstream anomaly in this report. The whole PacBio arm
had been designed as if these were shotgun/WGS reads, so coverage-shaped assumptions
(de novo assembly, whole-genome consensus, "0% N" as a quality signal) were all invalid.

---

## 3. De novo arm — hifiasm

### 3.1 Why it produced nothing

hifiasm is an OLC assembler and needs *dovetail* overlaps. In 128_5 and 203_3 every read
is contained within every other read, so there is no overlap graph to build. 211_0 is the
only sample with two amplicons that genuinely dovetail (26–5120 and 4410–6001, overlapping
at 4410–5120), and it is the only sample that assembles.

For 124_4 the reason is narrower: its two long groups *do* dovetail, but hifiasm reports
**62 overlaps, 0 strong, 62 weak** and builds unitigs only from strong overlaps. By
contrast 211_0 yields 548 overlaps, all strong.

### 3.2 Parameter sweep (7 configs × 4 samples)

Contigs produced (`p_ctg` segment count):

| config | flags | 124_4 | 128_5 | 203_3 | 211_0 |
|---|---|---|---|---|---|
| A_default | (none) | 0 | 0 | 0 | 0 |
| B_current | `-f0 -l0 --n-hap 1 --hg-size 10k` | 0 | 0 | 0 | 0 |
| **C_notips** | **B + `-n 0 --ctg-n 0`** | 0 | 0 | 0 | **1 (5977bp)** |
| D_notips_noclean | C + `-a 0` | 0 | 0 | 0 | 0 |
| E_smallk | C + `-k 31 -w 15` | 0 | 0 | 0 | 1 (5977bp) |
| F_primary | `--primary` variant + `-a 0` | 0 | 0 | 0 | 0 |
| G_maxovlp | C + `-a 0 -D 10 -N 500 --max-kocc 20000` | 0 | 0 | 0 | 0 |

Additional negative results:

- **`--hom-cov` sweep (1, 2, 3, 5, 8):** no effect on any sample. The
  `adjust_utg_by_primary` coverage floor was an initial hypothesis; forcing it to 0 still
  leaves the raw unitig graph empty, so the reads are dropped *upstream* of that filter.
  The per-sample `--hom-cov` code written for this hypothesis was removed as dead.
- **Read deduplication (`seqkit rmdup -s`):** no effect (124_4 40→26 reads, 128_5 75→57,
  203_3 19→13, 211_0 30→24; outcomes unchanged).
- **Contained-read removal:** no effect. Reducing 124_4 to its 12 non-contained reads
  (7 × 814–8002, 5 × 1963–9719) still yields 0 unitigs, with all 62 overlaps classed weak.
- **`-a 0` is actively harmful** — it breaks the one sample that works.

### 3.3 Settled configuration

```
HIFIASM_OPTS = -f0 -l0 --n-hap 1 --hg-size 10k -n 0 --ctg-n 0
```

`-n 0`/`--ctg-n 0` are the operative additions: the default `-n 3` removes tip unitigs of
≤3 reads, which on a 30-read library deletes essentially the whole graph.

### 3.4 The one successful assembly

211_0 → **5977 bp, 89.8% identity to HXB2 26–6001** (BLAST) — exactly the union of its two
amplicons. It is a correct but partial assembly and scores `output_valid=0` against the
8000 bp threshold.

> **Methodological note for the write-up:** minimap2's `asm5`/`asm20` presets *misreport*
> this contig badly — `asm5` returns no hit at all and `asm20` gives 32% identity. These
> presets are tuned for near-identical assemblies and break down at the ~90% identity of
> subtype A1/D reads against a subtype B reference. **Use BLAST for all identity claims in
> this project.**

### 3.5 Reporting change

Zero contigs is now recorded as a **result** (`exit 0`, `key_metric = "0 contigs (empty
overlap graph; amplicon input)"`) rather than a tool failure (`exit 1`), so the summary
distinguishes "structurally cannot assemble this input" from "hifiasm crashed".

---

## 4. Reference-guided arm — the HXB2 padding defect

`bcftools consensus` starts from the reference and only edits positions where a variant
was called. Positions with no read coverage silently retain the reference base, and no `N`
is emitted to flag them.

| sample | covered by reads | remainder of the "consensus" | old reported metric |
|---|---|---|---|
| 124_4 | 90.8% | 9% is HXB2 | `9719bp, 0.00% N` |
| 211_0 | 61.4% | **39% is HXB2** | `9719bp, 0.00% N` |
| 128_5 | 55.9% | **44% is HXB2** | `9719bp, 0.00% N` |
| 203_3 | 11.9% | **88% is HXB2** | `9719bp, 0.00% N` |

Worked example (203_3, reads cover only 27–1193):

```
Window 300-360 (inside covered region)
HXB2      : GAGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCG
203_3 cons: GAGCTGCATCCGGAGTTTTACAAGAACTGCTGACACAGAAGTTGCTGACGGGGACTTTCAG   <- real sample data

Window 5000-5060 (no reads at all)
HXB2      : AAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATTAGGGATTATGGAAAACAGATGGCAGGT
203_3 cons: AAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATTAGGGATTATGGAAAACAGATGGCAGGT   <- verbatim HXB2
```

**Fix:** zero-coverage intervals are derived from `samtools depth`, written as a BED, and
passed to `bcftools consensus -m`. Verified exact — the N count in the output matches the
BED interval total in all 8 runs (both arms × 4 samples).

**Validity criterion changed** from total length to *called* (non-N) bases ≥ 8000, because
a 9719 bp consensus that is 88% N should not score as a recovered genome.

---

## 5. Downstream impact on U3 / motif results

### 5.1 LTR coverage is one-LTR-only per sample

| sample | 5'U3 (1–455) | 3'U3 (9086–9550) |
|---|---|---|
| 124_4 | 0% | **100% (24 reads)** |
| 203_3 | 94% | 0% |
| 211_0 | 94% | 0% |
| **128_5** | **0%** | **0%** |

Both LTRs are identical in an integrated provirus, so having one is biologically
sufficient. But the consensus FASTA carries *both*, and the uncovered one is reference.

### 5.2 The U3 extraction takes the 5' LTR unconditionally

| sample | extracted U3 vs HXB2 |
|---|---|
| 124_4 | **IDENTICAL to HXB2** |
| 128_5 | **IDENTICAL to HXB2** |
| 203_3 | 68 differences (85.0% identity) — real |
| 211_0 | 61 differences (86.5% identity) — real |

### 5.3 Consequence for the motif comparison

| sample | FIMO hits in U3 | status |
|---|---|---|
| K03455.1 (HXB2) | RELA 350–359, RELA 364–373, SP1 388–396 | reference |
| **124_4** | RELA 350–359, RELA 364–373, SP1 388–396 | **HXB2's, verbatim** |
| **128_5** | RELA 350–359, RELA 364–373, SP1 388–396 | **HXB2's, verbatim** |
| 203_3 | RELA 364–373, SP1 388–396 | real |
| 211_0 | NFKB1 65–77 ×2, RELA 350–359, RELA 364–373, SP1 388–396 | real |

Two observations for the write-up:

1. **203_3 lacks RELA 350–359 that everything else has.** Because 203_3's U3 is genuine,
   that absence is a real biological finding — possibly a subtype A1/D variant NF-κB
   site. It is currently obscured by two artifact samples that "have" the site.
2. **124_4 is recoverable.** It has 100% coverage of the 3'LTR U3 with 24 reads; the
   extraction simply takes the uncovered 5' copy. Preferring whichever LTR has coverage
   would turn 124_4 into a genuine third data point. **128_5 is not recoverable** (no LTR
   coverage at either end) and should be excluded from motif analysis rather than
   reported as carrying HXB2's motifs.

On the concern that N-masking could hide real motifs: masking cannot hide a motif that was
*observed* — it only blanks sequence that was never sequenced. The one genuine risk is a
motif straddling a coverage boundary. Checked: 203_3 and 211_0 are missing only ~positions
1–26 of U3 and every hit sits at position ≥65, so masking clips no motif in this dataset.
**Reporting requirement:** "not sequenced" must stay distinguishable from "motif absent".

---

## 6. Resolved assembly design — three arms

| arm | reference | rationale |
|---|---|---|
| `minimap2_consensus` | HXB2 (K03455.1) | fixed subtype-B baseline; the coordinate system everything downstream expects |
| `minimap2_bestref` | closest LTR-complete Group M panel entry | measures the cost of forcing subtype A1/D reads onto a subtype B backbone |
| `hifiasm` | none (de novo) | reference-free comparator |

The two reference-guided arms are the **same pipeline over different references** — one
variable changes, so their rows are directly comparable.

Panel source: `scripts/assembly/illumina/shiver_setup/HIV1_COM_ref_alignment.fasta`
(LANL 2021 Group M compendium, 180 sequences, already vendored for the Illumina SHIVER arm).

---

## 7. Bugs found and fixed during verification

### 7.1 LTR-less reference panel discarded U3 reads

LANL compendium genomes conventionally span gag→nef and **omit the LTRs — only 22 of the
180 carry both U3 regions.** The first implementation selected
`D.KE.11.DEMD11KE003.KF716476` for 124_4, which covers only HXB2 **649–9594**. That
silently discarded 124_4's 24 3'LTR reads: **40 mapped → 16**.

For a U3 project this is disqualifying. The panel build now filters to LTR-complete
entries by reading gap content directly off the alignment: HXB2's ungapped positions give
the column range for each U3, and a candidate must be non-gap across ≥90% of both.
After filtering, **all four samples map 100% of reads in both arms.**

Cost of the filter — the LTR-complete subset is thin:

```
22 of 180 entries; subtypes A1(1) A4(1) A6(2) B(3) C(2) D(1) F1(1) G(2) H(1) J(1)
  02_AG(1) 06_cpx(1) 12_BF(1) 26_A5U(1) 27_cpx(1) 45_cpx(1) 60_BC(1)
```

Only **one D and one A1** reference survive. This is a property of the LANL compendium,
not something further tuning fixes.

### 7.2 Stale `.fai` produced a variant-free consensus that looked successful

On re-run, each `<sample>_reference.fasta` was overwritten but its `.fai` was not.
samtools/bcftools trust an existing index over the file, so `mpileup` failed with
`The sequence "..." was not found`, emitted an **empty VCF**, and the "consensus" came out
as bare masked reference with zero variants applied.

It passed silently because **`bash -c` does not inherit `set -o pipefail`** from the parent
script, so the failing `mpileup` in `mpileup | call` was masked by `call`'s exit 0 and the
`|| exit 1` never fired.

Three fixes: drop and regenerate the `.fai` whenever the reference is written; run the
tool functions under `bash -o pipefail -c`; and add an explicit guard that errors when
reads mapped but the VCF has no records.

**Self-consistency check now passing:** 128_5 selects HXB2 in *both* arms, and its two
consensus outputs are byte-identical — which is only true once variants are actually applied.

---

## 8. Final results

### 8.1 Assembly summary (`results/assembly/pacbio/summary.tsv`)

| tool | sample | exit | valid | key metric |
|---|---|---|---|---|
| minimap2_consensus | 124_4 | 0 | **1** | 9719bp, 8826bp called, 9.18% N, ref=HXB2 |
| minimap2_bestref | 124_4 | 0 | **1** | 9770bp, 8806bp called, 9.86% N, ref=D:KU168271 |
| hifiasm | 124_4 | 0 | 0 | 0 contigs (empty overlap graph; amplicon input) |
| minimap2_consensus | 128_5 | 0 | 0 | 9719bp, 5430bp called, 44.13% N, ref=HXB2 |
| minimap2_bestref | 128_5 | 0 | 0 | 9719bp, 5430bp called, 44.13% N, ref=B:K03455 |
| hifiasm | 128_5 | 0 | 0 | 0 contigs (empty overlap graph; amplicon input) |
| minimap2_consensus | 203_3 | 0 | 0 | 9719bp, 1153bp called, 88.13% N, ref=HXB2 |
| minimap2_bestref | 203_3 | 0 | 0 | 9807bp, 1131bp called, 88.46% N, ref=A1:KU168256 |
| hifiasm | 203_3 | 0 | 0 | 0 contigs (empty overlap graph; amplicon input) |
| minimap2_consensus | 211_0 | 0 | 0 | 9719bp, 5969bp called, 38.58% N, ref=HXB2 |
| minimap2_bestref | 211_0 | 0 | 0 | 9770bp, 5953bp called, 39.06% N, ref=D:KU168271 |
| hifiasm | 211_0 | 0 | 0 | 5977bp, 5977bp called, 0.00% N |

Runtimes 0.5–4 s per run; `peak_rss_mb` blank (no GNU `time` on this cluster).

### 8.2 Does subtype-matching help?

| sample | variants vs HXB2 | vs best reference | chosen reference |
|---|---|---|---|
| 203_3 | 149 | **90 (−40%)** | A1:KU168256 |
| 124_4 | 1097 | 1109 | D:KU168271 |
| 211_0 | 598 | 614 | D:KU168271 |
| 128_5 | 538 | 538 | B:K03455 (= HXB2) |

**Interpretation:** subtype-matching clearly helps 203_3 and is a wash for the others. The
limiting factor is the thin LTR-complete panel (one D, one A1) — requiring LTRs costs
candidate diversity.

An important supporting observation: against the **unfiltered** 180-entry panel with
`--secondary=no`, **HXB2 received zero primary alignments for all four samples**. HXB2 is
not the closest available reference for any of these samples; it is the closest one that
contains the LTRs.

Reference-selection scores against the unfiltered panel (summed matching bases):

| sample | top candidate | score | runner-up | score |
|---|---|---|---|---|
| 124_4 | D.KE.11.DEMD11KE003.KF716476 | 15936 | D.SE.12.077UG.MF373180 | 7342 |
| 128_5 | D.KE.11.DEMD11KE003.KF716476 | 120880 | (sole hit) | — |
| 203_3 | A1.CD.02.LA01AlPr.KU168256 | 7086 | 37_cpx.CM.00.00CMNYU926 | 477 |
| 211_0 | D.KE.11.DEMD11KE003.KF716476 | 28590 | A1.ES.15.100_117.KY496622 | 5726 |

Three samples match subtype D and one A1 — consistent with a Rakai cohort, where subtype
is not recorded in SRA metadata. **These are closest-reference matches, not subtype calls**;
the formal call remains jpHMM + IQ-TREE downstream.

---

## 9. Open items

1. **U3 extraction takes the wrong LTR** (§5.2). Preferring whichever LTR has read coverage
   recovers 124_4 as a real data point. Highest-value remaining fix.
2. **128_5 should be excluded from motif analysis** — no LTR coverage at either end.
3. **Subtyping is circular** — jpHMM currently subtypes an HXB2-padded consensus. For 203_3
   that was 88% HXB2. Re-run against the N-masked consensus.
4. **Motif results need regenerating** downstream of the corrected assemblies.
5. **Optional:** iterative consensus refinement (re-map reads to the called consensus, 2–3
   rounds) to remove residual reference bias. Not implemented.
6. **Optional:** reference-free amplicon consensus (multiple-align each read stack,
   majority-vote) as a de novo arm that *can* work on amplicon data, unlike hifiasm.
   `mafft` is installed; `abpoa`/`spoa` would be better suited.

---

## 10. Files and reproduction

| item | path |
|---|---|
| Assembly harness | `scripts/assembly/pacbio/compare_assembly_pacbio.sh` |
| Shared helpers | `scripts/common/lib_compare.sh` |
| Reference panel (built) | `results/assembly/pacbio/refpanel/panel.fasta` |
| LTR-complete IDs | `results/assembly/pacbio/refpanel/ltr_complete.ids` |
| HXB2-arm output | `results/assembly/pacbio/minimap2_consensus_out/` |
| Subtype-arm output | `results/assembly/pacbio/minimap2_bestref_out/` |
| De novo output | `results/assembly/pacbio/hifiasm_out/` |
| Per-sample mask BEDs | `<arm>_out/<sample>.uncovered.bed` |
| Chosen reference per sample | `<arm>_out/<sample>_reference.txt` |
| Summary table | `results/assembly/pacbio/summary.tsv` |

Run (never on the login node):

```bash
THREADS=4 srun -p shared -c 4 --mem=8G -t 30 scripts/assembly/pacbio/compare_assembly_pacbio.sh
```

Restrict arms via `ASSEMBLY_TOOLS`, e.g.
`ASSEMBLY_TOOLS="minimap2_bestref" srun ... compare_assembly_pacbio.sh`.
