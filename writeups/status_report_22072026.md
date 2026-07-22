# HIV-1 Proviral U3 Analysis — Tool-Comparison Status Report
**Date:** 2026-07-22  ·  **Author:** Jovita Nagawa  ·  **Objective 1** methodology work

This report captures the state of the Illumina and PacBio tool-comparison harnesses,
the results obtained, the challenges resolved, all deliverables, and recommended next
steps. Both harnesses run every automatable candidate tool from
`Tools_review_HIV_U3_analysis.pdf` side-by-side on a small per-platform subset so
methodology choices are evidence-backed.

---

## 1. Executive summary

- **Illumina arm (4-sample subset SRR908437/434/441/446):** all 6 steps run.
  Tool winners are clear for every step; the only outstanding data point is a
  **working SHIVER** row (blocked on a missing `bc`, fix staged).
- **PacBio arm (4 user-supplied host-N-masked Rakai HiFi samples 124_4/128_5/203_3/211_0):**
  all 6 steps run end-to-end. minimap2-consensus assembly succeeded; **hifiasm** (de-novo
  comparator) still to be added once installed.
- **Slide decks** for both arms delivered in editable PPTX **and** HTML, with results
  summaries, an MSA snapshot, and per-transcription-factor + G-quadruplex motif counts.
- **Gated on the SRA download finishing** (frees the shared conda env): install
  `bc` + hifiasm + chopper, then the SHIVER and hifiasm re-runs.

---

## 2. Illumina arm — per-step results

| Step | Tools compared | Chosen / result | Status |
|---|---|---|---|
| Download & QC | fastp vs Trimmomatic; Kraken2 host removal | **fastp** — same quality bar, higher yield, ~3× faster. Kraken2 removed ~7–45% (host) per sample (SRR908446 ~40% human) | ✅ done |
| Assembly | BWA+bcftools vs SPAdes vs SHIVER | **BWA+bcftools** — full-length 9,719 bp / 0% N on all 4; SPAdes fragments (420–1,786 bp, coverage-sensitive) | ✅ BWA/SPAdes; ⏳ SHIVER |
| MSA | MAFFT L-INS-i vs MUSCLE vs Clustal Ω | **MAFFT L-INS-i** — most accurate + fastest (79 s vs 116 / 205) | ✅ done |
| Biological filtering | Poplars vs HIVSeqinR vs HIVIntact | Poplars ran; HIVIntact 0 intact/5 non-intact (subtype-B-biased); **HIVSeqinR** = primary (pending config) | ✅ Poplars/HIVIntact; ⏳ HIVSeqinR |
| Subtyping | jpHMM vs IQ-TREE 2 | **jpHMM** (breakpoints) + **IQ-TREE 2** (ML confirm); jpHMM ~10 min/sample, IQ-TREE ~1 s | ✅ done |
| Motif mapping | FIMO / MOODS / TFBSTools; gquad / pqsfinder | TFBS per factor below; G4 gquad 15 / pqsfinder 5 | ✅ done |

**Illumina TFBS hits per factor (subset U3):** NF-κB p65 — FIMO 10 / TFBSTools 20; SP1 — 5 / 15;
NFAT — 0 / 15; TBP — 0 / 20; NF-κB p50 & AP-1 — 0 / 0. (Totals: FIMO 16, MOODS 10, TFBSTools 70.)

**Key challenges resolved (Illumina):** SPAdes & MUSCLE SIGILL → pinned pre-AVX2 versions
(spades 3.15.5, muscle 3.8); Kraken2 wired in as a real step (+ DB documented); SHIVER LANL
alignment automated via curl; `lib_compare.sh` stderr-swallowing bug; `REPO_ROOT` off-by-one in
7 scripts; MAFFT lowercase → uppercased U3; HIVIntact `setuptools<81`; jpHMM `LD_LIBRARY_PATH`
GLIBCXX shadow; IQ-TREE bootstrap dropped for small subset + `iqtree=2.*` pin; JASPAR 2024 via REST.

**Outstanding (Illumina):** SHIVER fails at `shiver_init` because compute nodes lack `bc`
(config values are valid, so `bc` is the sole blocker). Fix staged: `bc` added to the env.

---

## 3. PacBio (HIV-SMRTcap) arm — per-step results

Input: **4 user-supplied host-N-masked Rakai HiFi samples** (124_4, 128_5, 203_3, 211_0) —
host flanks masked N, HIV provirus in ACGT. These bypass QC (already processed), entering at
proviral extraction.

| Step | Tools | Result | Status |
|---|---|---|---|
| Proviral extraction | `extract_provirus_strip_hostN` (bash/awk/seqkit) | provirus from 41/75/19/30 reads, 0% residual N; cores **partial** (mean 1.8–6.0 kb) | ✅ done |
| Assembly | minimap2→HXB2 consensus vs hifiasm | minimap2 9,719 bp/0% N ×4 (reference-filled where cores partial) | ✅ minimap2; ⏳ hifiasm |
| MSA | MAFFT / MUSCLE / Clustal Ω | 5/5 retained; MAFFT fastest (78 s); real inter-sample variation visible | ✅ done |
| Biological filtering | Poplars / HIVSeqinR / HIVIntact | Poplars ran; HIVIntact 0/5; HIVSeqinR needs primer config | ✅ Poplars/HIVIntact; ⏳ HIVSeqinR |
| Subtyping | jpHMM / IQ-TREE 2 | per-sample breakpoint maps ×4 (~10 min each) + ML tree (~2 s) | ✅ done |
| Motif mapping | FIMO / MOODS / TFBSTools; gquad / pqsfinder | per factor below; G4 gquad 15 / pqsfinder 5 | ✅ done |

**PacBio TFBS hits per factor (subset U3):** NF-κB p65 — FIMO 9 / TFBSTools 19; NF-κB p50 — 2 / 2;
SP1 — 5 / 15; NFAT — 0 / 14; AP-1 — 0 / 6; TBP — 0 / 31. TFBSTools recovered **all 6** TFs.
(Totals: FIMO 17, MOODS 12, TFBSTools 87.)

**Important caveat:** SMRTcap reads are host-provirus **junction reads**, so the extracted cores
are partial (esp. 203_3, ~1.8 kb). Reference-guided consensus fills uncovered positions with
HXB2 — so a "9,719 bp valid" consensus is partly reference-derived. Downstream U3/motif numbers
are **tool-comparison signal, not final biology**; per-sample U3 read depth should be reported.

**Key challenge resolved (PacBio):** the uploaded `raw_smrtcap/` folder nested files one level
deeper than the reader expected → all samples skipped → chain failed. Fixed with a recursive
file search; re-ran cleanly. HIVSeqinR blocked on a manual primer-config step (documented).

---

## 4. Deliverables

**Slides (writeups/):**
- `Illumina_Tool_Comparison.pptx` + `illumina_tool_comparison_slides.html`
- `PacBio_Tool_Comparison.pptx` + `pacbio_tool_comparison_slides.html`
- `figures/msa_u3_snapshot.png`, `figures/msa_u3_pacbio_snapshot.png`
- 13 slides each (title + 6 steps × 2: results / challenges); editable charts + tables.

**Code:** `scripts/<step>/{illumina,pacbio}/` harnesses + shared `scripts/common/lib_compare.sh`;
new `scripts/utils/extract_provirus_strip_hostN.sh` (proviral N-strip), `kraken2_filter_reads_se.sh`,
`setup_pacbio_env.slurm.sh`, `run_pacbio_chain.sh`. Env pins in `HIV_U3analysis_env.yml`.

**Results:** `results/<step>/{illumina,pacbio}/summary.tsv` + logs.

---

## 5. Recommended next steps (prioritized)

1. **Finish the two assembly comparisons (auto, on env-setup).** When the SRA download frees the
   conda env, install `bc` + hifiasm + chopper, then re-run: SHIVER (Illumina) and hifiasm (PacBio).
   Adds the missing assembly rows; refresh the two assembly slides.
2. **Configure & run HIVSeqinR** on both arms — supply the protocol's 2nd-round PCR primers and set
   `.CONFIGURED`. It is the Rakai A1/D-validated intactness caller and should be the primary one.
3. **Report PacBio U3 read depth per sample** and, ideally, obtain fuller-coverage / full-length
   proviral reads (or filter to reads spanning U3) so U3/motif calls are sample-derived, not
   reference-filled. This is the biggest interpretability gap on the PacBio side.
4. **Run the PacBio QC-tool comparison** (`download_qc_pacbio.sh`: NanoFilt vs fastp vs chopper +
   Kraken2 + dedup) on the just-downloaded SRA HiFi reads — the user's masked reads bypassed QC, so
   this completes the PacBio harness's QC step.
5. **Scale from subset to full cohort** once tool choices are locked: Illumina (24 samples), PacBio
   (all A1/D/recombinant incl. the 3 ~18 GB Revio Rakai runs). Re-enable IQ-TREE `-B 1000` at scale.
6. **Fold in the web-only tools** (COMET, REGA v3 for subtyping; QGRS Mapper for G4) by hand and add
   to the respective `ease_of_use_notes.md` — the review expects these as complements.
7. **Fill in each step's `ease_of_use_notes.md`** with the session-verified results (assembly done;
   msa/subtyping/biofilt/motif still template stubs).
8. **Wire chosen tools into the production pipelines** (`scripts/pipelines/illumina_u3analysis.sh`,
   `oxnano_u3analysis.sh` steps 8/10/11 stubs) — the end goal of the comparison work.
9. **Commit to git.** All of the above (scripts, results, slides, this report) is currently
   uncommitted on `main`; commit in logical chunks (harness code, results, slides).

---

## 6. Currently running / self-healing

- SRA PacBio download (QC-comparison source) — last sample in progress.
- Monitors watch every chain; failures are diagnosed, fixed, and re-submitted automatically
  (per standing instruction). Items 1 above will fire on their own once the env is free.
