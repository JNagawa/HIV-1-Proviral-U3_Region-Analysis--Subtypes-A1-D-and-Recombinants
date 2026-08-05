# Aligner accuracy benchmark — collinear input (snapshot 2026-08-04)

Snapshot of the `minimap2_consensus` MSA outputs from **before** the 2026-08-05 indel
calling and mask fixes. Retained because it is an accidental ground-truth benchmark that
cannot be reproduced once those fixes are in place.

**Why it is ground truth.** At the time these ran, `bcftools call -c` was calling zero
indels (see `writeups/status_report_05082026.md` §2), so every consensus was exactly
9719 bp — HXB2 with substitutions and N-masking only, no length change. The five input
sequences are therefore *strictly collinear* with HXB2, and the single correct alignment
is the identity: 9719 columns, zero gaps.

**Result.**

| tool | columns | gaps/record | verdict |
|---|---|---|---|
| mafft | 9719 | 0 | correct |
| muscle | 9719 | 0 | correct |
| clustalo | 9768 | 49 | 49 spurious indels per sequence |

Clustal Omega treats `N` as an alignable residue and slides long N-runs against real
sequence — 203_3 (88% N) is worst affected. Full analysis in
`writeups/status_report_05082026.md` §6.

**Do not regenerate.** After the fixes, consensus lengths differ (9728 / 9751 / 9730 /
9720) because indels are now applied, so the identity alignment is no longer the correct
answer and this control no longer exists.
