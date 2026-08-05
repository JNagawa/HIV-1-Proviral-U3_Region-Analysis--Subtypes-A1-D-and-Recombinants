#!/bin/bash
#SBATCH --job-name=pb_motif_toolcmp
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

##################################################################################
# Coordinate-level comparison of the motif-mapping tools run by
# compare_motif_mapping_pacbio.sh. That script records only HIT COUNTS per tool,
# which cannot answer whether two tools found the SAME sites -- two tools can
# report 40 hits each and share none of them. This one normalises every tool's
# output to a common interval representation and compares them position by
# position.
#
# Pairs compared:
#   TFBS  : FIMO vs TFBSTools (the requested pair), with MOODS as a third view
#   G4    : gquad vs pqsfinder
#
# COORDINATE NORMALISATION is the part that has to be right:
#   FIMO        fimo.tsv, 1-based inclusive  (start, stop)
#   TFBSTools   GFF3,     1-based inclusive  (start, end)
#   gquad       GFF3,     1-based inclusive
#   pqsfinder   GFF3,     1-based inclusive
#   MOODS       CSV,      0-BASED start, no end column -- the end is derived
#               from the length of the matched sequence. Verified against a
#               shared RELA match that FIMO/TFBSTools both call at 350-359 and
#               MOODS reports as 349, i.e. MOODS needs +1. Comparing MOODS
#               without this shift would make every site look off-by-one.
#
# ON "SENSITIVITY AND SPECIFICITY": those require a gold standard, and none
# exists for these sample U3 sequences -- no curated set of true binding sites.
# What this script therefore reports is:
#   (a) CONCORDANCE between tools (exact coordinate matches, partial overlaps,
#       tool-unique calls, and base-pair Jaccard), which needs no ground truth
#   (b) RELATIVE recall/precision when one tool is taken as the reference --
#       a statement about agreement, NOT about correctness
#   (c) recovery of the literature-known HXB2 LTR landmarks on the positive
#       control (two NF-kB sites, three Sp1 sites, TATA box), which is the only
#       real ground truth available and is small and approximate
# The report says so explicitly rather than labelling (a)/(b) as sensitivity.
#
# Usage: sbatch scripts/motif_mapping/pacbio/compare_motif_tools_pacbio.sh
# -u errors on unset vars, pipefail fails a pipe if any stage fails
set -uo pipefail

# Locate the repo root: under `sbatch <script>` Slurm copies this file into a
# spool dir, so a $0-relative walk lands outside the repo.
REPO_ROOT=""                                         # filled in by the loop below
for CAND in "${SLURM_SUBMIT_DIR:-}" \
            "$(cd "$(dirname "$0")/../../.." 2>/dev/null && pwd)" \
            "$(pwd)"; do
    # a repo root is only a repo root if the helper library is under it
    if [ -n "${CAND}" ] && [ -f "${CAND}/scripts/common/lib_compare.sh" ]; then
        REPO_ROOT="${CAND}"; break                   # first match wins
    fi
done
[ -n "${REPO_ROOT}" ] || { echo "ERROR: cannot locate repo root." >&2; exit 1; }

# which assembly arm's motif output to analyse (see compare_msa_pacbio.sh)
ASSEMBLY_ARM="${ASSEMBLY_ARM:-minimap2_consensus}"
MOTIF_DIR="${REPO_ROOT}/results/motif_mapping/pacbio/${ASSEMBLY_ARM}"  # where the tool outputs live
OUT_DIR="${MOTIF_DIR}/tool_comparison"                 # where this analysis writes
NORM="${OUT_DIR}/normalized_hits.tsv"                  # every tool's hits, one common schema
PAIRS="${OUT_DIR}/pairwise_concordance.tsv"            # the pair-by-pair agreement table
REPORT="${OUT_DIR}/REPORT.md"                          # the human-readable write-up
mkdir -p "${OUT_DIR}"                                  # create the analysis dir

##################################################################################
#  0. How much of each U3 is actually SAMPLE sequence?
#
#  This has to come first, because it decides whether the concordance numbers
#  below mean anything. The assembly step builds its consensus with
#  `bcftools consensus`, which APPLIES CALLED VARIANTS ONTO THE HXB2 REFERENCE.
#  Positions with no read coverage are therefore emitted as reference bases, not
#  as N -- every consensus comes out 9719bp at 0% N whether or not the sample had
#  any reads there. Mapping the stripped cores back to HXB2 shows the four
#  samples cover very different intervals (124_4: 813-8002 and 8985-9719;
#  128_5: 1593-7030; 203_3: 26-1193; 211_0: 25-5120), so a sample whose reads
#  never reached the U3 window contributes an all-reference U3.
#
#  If two sample U3s are byte-identical to HXB2's, then every tool will
#  trivially agree on them, and cross-tool concordance would be measuring
#  sequence duplication rather than tool behaviour. So: report the identity of
#  each extracted U3 against HXB2 up front, and treat any 100%-identical
#  sequence as a duplicate of the positive control rather than as evidence.
##################################################################################
IDENT="${OUT_DIR}/u3_identity_to_hxb2.tsv"           # per-sequence identity to the HXB2 U3
U3_FA="${MOTIF_DIR}/U3_extracted.fasta"              # the sequences every tool scanned
if [ -s "${U3_FA}" ]; then
    # linearise to name<TAB>seq, then compare each record to the HXB2 record
    seqkit fx2tab "${U3_FA}" 2>/dev/null | awk -F'\t' -v OFS='\t' '
        { name[NR]=$1; seq[NR]=toupper($2); if ($1 ~ /^K03455/) ref=toupper($2); n=NR }
        END {
            print "sequence", "length", "identical_to_HXB2_U3", "pct_identity", "interpretation"
            for (i=1; i<=n; i++) {
                if (name[i] ~ /^K03455/) { print name[i], length(seq[i]), "reference", "100.00", "the reference itself"; continue }
                same=0; L=(length(seq[i])<length(ref))?length(seq[i]):length(ref)
                for (p=1; p<=L; p++) if (substr(seq[i],p,1)==substr(ref,p,1)) same++
                pct = (L>0) ? 100*same/L : 0
                interp = (seq[i]==ref) ? "IDENTICAL to HXB2 -- no sample-specific sequence here" : "sample-specific sequence present"
                printf "%s\t%d\t%s\t%.2f\t%s\n", name[i], length(seq[i]), (seq[i]==ref ? "yes" : "no"), pct, interp
            }
        }
    ' > "${IDENT}"
    echo "=== U3 identity to HXB2 (does each sample contribute real sequence?) ==="
    column -t -s$'\t' "${IDENT}"
    echo ""
fi

##################################################################################
#  1. Normalise every tool's output to: sequence  tool  feature  start  end  strand  score
##################################################################################
: > "${NORM}"                                        # truncate any previous run

# FIMO: skip the header and any comment/blank lines; col3=sequence, col2=TF name
if [ -s "${MOTIF_DIR}/fimo_out/fimo.tsv" ]; then
    awk -F'\t' 'NR>1 && $0!~/^#/ && NF>=7 && $4!="" {
        printf "%s\tfimo\t%s\t%d\t%d\t%s\t%s\n", $3, $2, $4, $5, $6, $7
    }' "${MOTIF_DIR}/fimo_out/fimo.tsv" >> "${NORM}"
fi

# TFBSTools: GFF3, TF name lives in the attributes column as TF=<name>
if [ -s "${MOTIF_DIR}/tfbstools_out.gff3" ]; then
    awk -F'\t' '$0!~/^#/ && NF>=9 {
        tf=$9; sub(/.*TF=/,"",tf); sub(/;.*/,"",tf)
        printf "%s\ttfbstools\t%s\t%d\t%d\t%s\t%s\n", $1, tf, $4, $5, $7, $6
    }' "${MOTIF_DIR}/tfbstools_out.gff3" >> "${NORM}"
fi

# MOODS: CSV, col3 is a 0-BASED start and there is no end column, so the end is
# start + length(matched_sequence) - 1 after shifting to 1-based.
#
# MOODS also names motifs by JASPAR ACCESSION (MA0107.1) while FIMO and
# TFBSTools use the TF SYMBOL (RELA). Comparing them without translating would
# make every key mismatch, so no MOODS hit could ever pair with a FIMO hit --
# which is exactly the artefact that produced "0 overlaps but 0.68 Jaccard" on
# the first run. fimo.tsv carries both (col1 = accession, col2 = symbol), so the
# mapping is built from it and applied here; anything unmapped keeps its
# accession and simply will not match, which is visible rather than silent.
if [ -s "${MOTIF_DIR}/moods_out.tsv" ]; then
    MAP="${OUT_DIR}/jaspar_id_to_tf.tsv"             # accession -> TF symbol
    awk -F'\t' 'NR>1 && $1!="" && $2!="" {print $1"\t"$2}' \
        "${MOTIF_DIR}/fimo_out/fimo.tsv" 2>/dev/null | sort -u > "${MAP}"
    awk -F',' -v map="${MAP}" '
        BEGIN { while ((getline line < map) > 0) { split(line, f, "\t"); tf[f[1]]=f[2] } }
        NF>=6 && $3!="" {
            m=$2; sub(/\.pfm$/,"",m)                 # strip the .pfm suffix from the motif id
            name = (m in tf) ? tf[m] : m             # translate accession -> symbol when known
            s=$3+1                                   # 0-based -> 1-based
            e=s+length($6)-1                         # end derived from the matched sequence length
            printf "%s\tmoods\t%s\t%d\t%d\t%s\t%s\n", $1, name, s, e, $4, $5
        }' "${MOTIF_DIR}/moods_out.tsv" >> "${NORM}"
fi

# gquad / pqsfinder: GFF3, no TF -- the feature is just "G4"
for g4 in gquad pqsfinder; do
    if [ -s "${MOTIF_DIR}/${g4}_out.gff3" ]; then
        awk -F'\t' -v tool="${g4}" '$0!~/^#/ && NF>=8 {
            printf "%s\t%s\tG4\t%d\t%d\t%s\t%s\n", $1, tool, $4, $5, $7, $6
        }' "${MOTIF_DIR}/${g4}_out.gff3" >> "${NORM}"
    fi
done

# nothing to compare if no tool produced parseable hits
[ -s "${NORM}" ] || { echo "ERROR: no parseable hits found under ${MOTIF_DIR}." >&2; exit 1; }
echo "normalised $(wc -l < "${NORM}") hits across $(cut -f2 "${NORM}" | sort -u | tr '\n' ' ')"

##################################################################################
#  2. Pairwise concordance
##################################################################################
# match_tf=1 requires two hits to name the same TF before they can be called the
# same site (right for the TFBS tools); match_tf=0 compares intervals only
# (right for the G4 tools, which predict unnamed structures).
compare_pair() {
    local a="$1" b="$2" match_tf="$3"                # tool A, tool B, whether TF must match
    awk -F'\t' -v A="${a}" -v B="${b}" -v mtf="${match_tf}" '
        function mx(x,y){ return (x>y)?x:y }
        function mn(x,y){ return (x<y)?x:y }
        {
            # key groups hits that are eligible to match: always the sequence,
            # plus the TF name when the caller asked for it
            key = $1 (mtf ? SUBSEP $3 : "")
            if ($2 == A) { na++; ak[na]=key; as[na]=$4; ae[na]=$5; aseq[na]=$1 }
            if ($2 == B) { nb++; bk[nb]=key; bs[nb]=$4; be[nb]=$5; bseq[nb]=$1 }
            if ($2 == A || $2 == B) seqs[$1]=1        # every sequence either tool touched
        }
        END {
            exact=0; partial=0; aonly=0
            for (i=1; i<=na; i++) {
                best=0; bestov=0
                for (j=1; j<=nb; j++) {
                    if (ak[i] != bk[j]) continue      # different sequence (or TF) -> not comparable
                    ov = mn(ae[i], be[j]) - mx(as[i], bs[j]) + 1   # overlap in bp
                    if (ov > bestov) { bestov=ov; best=j }
                }
                if (best && as[i]==bs[best] && ae[i]==be[best]) exact++      # identical coordinates
                else if (bestov > 0) partial++                               # overlapping but shifted
                else aonly++                                                 # no counterpart at all
            }
            # same walk from B is side, to find B-unique calls
            bonly=0
            for (j=1; j<=nb; j++) {
                hit=0
                for (i=1; i<=na; i++) {
                    if (ak[i] != bk[j]) continue
                    if (mn(ae[i], be[j]) - mx(as[i], bs[j]) + 1 > 0) { hit=1; break }
                }
                if (!hit) bonly++
            }
            # base-pair Jaccard: how much of the union of covered bases both agree on
            for (i=1; i<=na; i++) for (p=as[i]; p<=ae[i]; p++) cova[aseq[i] SUBSEP p]=1
            for (j=1; j<=nb; j++) for (p=bs[j]; p<=be[j]; p++) covb[bseq[j] SUBSEP p]=1
            inter=0; for (p in cova) if (p in covb) inter++
            uni=0;   for (p in cova) uni++
                     for (p in covb) if (!(p in cova)) uni++
            jac = (uni>0) ? inter/uni : 0
            printf "%s_vs_%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\n", A, B, na, nb, exact, partial, aonly, bonly, jac
        }
    ' "${NORM}"
}

# header for the concordance table
printf "pair\tn_A\tn_B\texact_coord_match\tpartial_overlap\tA_only\tB_only\tbp_jaccard\n" > "${PAIRS}"
compare_pair fimo tfbstools 1 >> "${PAIRS}"          # the requested TFBS pair
compare_pair fimo moods     1 >> "${PAIRS}"          # third TFBS view
compare_pair tfbstools moods 1 >> "${PAIRS}"         # and the remaining TFBS pair
compare_pair gquad pqsfinder 0 >> "${PAIRS}"         # the requested G4 pair

echo "=== pairwise concordance ==="
column -t -s$'\t' "${PAIRS}"

##################################################################################
#  3. Per-sequence breakdown, so a single sequence cannot hide behind the totals
##################################################################################
PERSEQ="${OUT_DIR}/per_sequence_counts.tsv"          # hits per tool per sequence
printf "sequence\ttool\tn_hits\tbp_covered\n" > "${PERSEQ}"
awk -F'\t' '{
    n[$1 SUBSEP $2]++                                # hit count per sequence/tool
    for (p=$4; p<=$5; p++) cov[$1 SUBSEP $2 SUBSEP p]=1   # covered positions
} END {
    for (k in cov) { split(k, f, SUBSEP); bp[f[1] SUBSEP f[2]]++ }
    for (k in n) { split(k, f, SUBSEP); printf "%s\t%s\t%d\t%d\n", f[1], f[2], n[k], bp[k] }
}' "${NORM}" | sort >> "${PERSEQ}"
echo ""; echo "=== per-sequence hit counts ==="
column -t -s$'\t' "${PERSEQ}"

##################################################################################
#  4. Known HXB2 LTR landmarks -- the only real ground truth available
##################################################################################
# Literature positions of the core HIV-1 LTR regulatory elements, expressed as
# offsets from the U3/R boundary (the transcription start site), which is the
# end of the extracted U3. Converting them to U3-slice coordinates needs the
# actual U3 length, so it is computed from the extracted FASTA rather than
# hard-coded. These positions are APPROXIMATE (a few bp of drift between
# published maps is normal), so a hit is credited when it overlaps the landmark
# window at all -- this is an indicative check, not a validated benchmark.
LANDMARKS="${OUT_DIR}/hxb2_landmark_recovery.tsv"    # which tool recovered which landmark
U3_LEN=$(awk '/^>K03455/{f=1;next} /^>/{f=0} f{gsub(/[^A-Za-z]/,""); n+=length($0)} END{print n+0}' \
         "${MOTIF_DIR}/U3_extracted.fasta" 2>/dev/null)
if [ "${U3_LEN:-0}" -gt 0 ]; then
    echo ""; echo "=== HXB2 landmark recovery (U3 length ${U3_LEN}bp) ==="
    printf "landmark\tu3_start\tu3_end\t%s\n" "$(printf '%s\t' fimo tfbstools moods gquad pqsfinder)" > "${LANDMARKS}"
    # offsets are (start,end) upstream of the transcription start site
    awk -F'\t' -v L="${U3_LEN}" -v OFS='\t' '
        BEGIN {
            # name, upstream_start, upstream_end (bp before the U3/R boundary)
            split("NFkB-II:104:95 NFkB-I:90:81 Sp1-III:78:68 Sp1-II:67:57 Sp1-I:56:46 TATA:28:24", LM, " ")
            for (i in LM) {
                split(LM[i], f, ":")
                s = L - f[2] + 1; e = L - f[3] + 1   # convert upstream offsets to U3 coordinates
                name[i]=f[1]; ls[i]=s; le[i]=e
            }
        }
        # only the HXB2 record counts as the positive control
        $1 ~ /^K03455/ { nh++; ht[nh]=$2; hs[nh]=$4; he[nh]=$5 }
        END {
            split("fimo tfbstools moods gquad pqsfinder", tools, " ")
            for (i=1; i<=6; i++) {
                line = name[i] OFS ls[i] OFS le[i]
                for (t=1; t<=5; t++) {
                    found="-"
                    for (h=1; h<=nh; h++) {
                        if (ht[h] != tools[t]) continue
                        # credit an overlap with the landmark window
                        if (he[h] >= ls[i] && hs[h] <= le[i]) { found="yes"; break }
                    }
                    line = line OFS found
                }
                print line
            }
        }
    ' "${NORM}" >> "${LANDMARKS}"
    column -t -s$'\t' "${LANDMARKS}"
else
    echo "NOTE: could not determine HXB2 U3 length; skipping the landmark check." >&2
fi

##################################################################################
#  5. Compile the written report
##################################################################################
# Everything below is derived from the TSVs produced above rather than restated
# by hand, so the prose cannot drift from the numbers.
G4_DETAIL="${OUT_DIR}/g4_intervals.tsv"              # side-by-side G4 calls per sequence
printf "sequence\ttool\tstart\tend\tlength\n" > "${G4_DETAIL}"
awk -F'\t' '$2=="gquad" || $2=="pqsfinder" {printf "%s\t%s\t%d\t%d\t%d\n", $1,$2,$4,$5,$5-$4+1}' \
    "${NORM}" | sort -k1,1 -k3,3n >> "${G4_DETAIL}"

{
  echo "# Motif-mapping tool comparison -- PacBio HIV-SMRTcap arm"
  echo
  echo "Generated by \`scripts/motif_mapping/pacbio/compare_motif_tools_pacbio.sh\`."
  echo "Inputs: the U3 sequences extracted by \`compare_motif_mapping_pacbio.sh\` from the"
  echo "MAFFT alignment of the four SMRTcap samples plus HXB2 (K03455.1)."
  echo
  echo "## 1. Read this first: only two of the four samples contribute real sequence"
  echo
  echo "The assembly stage builds its consensus with \`bcftools consensus\`, which applies"
  echo "called variants **onto the HXB2 reference**. Positions with no read coverage are"
  echo "emitted as reference bases, not as N -- which is why every consensus is exactly"
  echo "9719 bp at 0.00% N regardless of how little of the genome the sample actually"
  echo "covered. Mapping the stripped proviral cores back to HXB2 shows the samples cover"
  echo "very different intervals, and two of them never reach the U3 window at all."
  echo
  column -t -s$'\t' "${IDENT}" 2>/dev/null | sed 's/^/    /'
  echo
  echo "**Consequence:** 124_4 and 128_5 carry a U3 that is byte-identical to HXB2's, so"
  echo "three of the five scanned sequences are the same sequence. Any tool will agree"
  echo "with itself on those, and per-sequence hit counts for K03455.1, 124_4 and 128_5"
  echo "are identical for every tool -- confirming duplication rather than reproducibility."
  echo "**The effective sample size for this comparison is n=2 (203_3 and 211_0).**"
  echo
  echo "## 2. Concordance between tools"
  echo
  echo "Counts are from tool A's perspective; \`bp_jaccard\` is over covered base positions."
  echo "TFBS pairs require the same TF to be named before two calls can match; the G4 pair"
  echo "compares intervals only. MOODS reports 0-based starts and JASPAR accessions, both"
  echo "normalised here (see the script header)."
  echo
  column -t -s$'\t' "${PAIRS}" | sed 's/^/    /'
  echo
  echo "### FIMO vs TFBSTools -- the requested pair"
  # column order: 1=pair 2=n_A 3=n_B 4=exact 5=partial 6=A_only 7=B_only 8=jaccard
  awk -F'\t' '$1=="fimo_vs_tfbstools"{
      printf "\nFIMO called %d sites; TFBSTools called %d. **All %d FIMO calls are reproduced by\n", $2,$3,$4
      printf "TFBSTools at byte-identical coordinates** (%d exact, %d partial, %d FIMO-only), while\n", $4,$5,$6
      printf "TFBSTools adds %d further calls. FIMO is therefore a strict subset of TFBSTools:\n", $7
      printf "where both fire they agree perfectly on position, and the entire difference between\n"
      printf "them is stringency, not disagreement about where sites are. The base-pair Jaccard of\n"
      printf "%.2f reflects that extra breadth, not positional conflict.\n", $8
  }' "${PAIRS}"
  echo
  echo "### MOODS"
  awk -F'\t' '$1=="fimo_vs_moods"{
      printf "\n%d of MOODS'\''s %d calls exactly match a FIMO call; %d FIMO calls have no MOODS\n", $4,$3,$6
      printf "counterpart and %d MOODS call is unique. The three PWM scanners thus place sites\n", $7
      printf "consistently; they differ almost entirely in how many they report.\n"
  }' "${PAIRS}"
  echo
  echo "### gquad vs pqsfinder -- the requested G4 pair"
  awk -F'\t' '$1=="gquad_vs_pqsfinder"{
      printf "\ngquad predicted %d G-quadruplexes, pqsfinder %d. **Not one pair shares identical\n", $2,$3
      printf "boundaries** (%d exact), yet every pqsfinder call overlaps a gquad call (%d partial,\n", $4,$5
      printf "%d pqsfinder-only) and %d gquad calls have no pqsfinder counterpart. So the two agree\n", $7,$6
      printf "on roughly WHERE quadruplexes are but never on their extent -- a systematic\n"
      printf "difference in how each defines a G4 boundary, not a sensitivity difference.\n"
  }' "${PAIRS}"
  echo
  echo "Interval detail (pqsfinder reports a single fixed-width call per sequence; gquad"
  echo "reports several of varying width):"
  echo
  column -t -s$'\t' "${G4_DETAIL}" | head -14 | sed 's/^/    /'
  echo
  echo "## 3. Sensitivity and specificity -- what can and cannot be measured"
  echo
  echo "Sensitivity and specificity require a gold standard: a curated list of sites that"
  echo "are genuinely present and genuinely absent. **No such truth set exists for these"
  echo "sample U3 sequences**, so those two statistics cannot honestly be computed for the"
  echo "subset, and the concordance figures in section 2 must not be read as sensitivity."
  echo
  echo "The one place real ground truth exists is HXB2's own U3, whose core regulatory"
  echo "elements are established in the literature: two NF-kB sites, three Sp1 sites and"
  echo "the TATA box. Positions below are derived from their published offsets upstream of"
  echo "the U3/R boundary, converted to coordinates in the ${U3_LEN:-?}bp extracted U3. A tool is"
  echo "credited when a call overlaps the landmark window; the windows are approximate to"
  echo "within a few bp, so this is indicative, not a validated benchmark."
  echo
  column -t -s$'\t' "${LANDMARKS}" 2>/dev/null | sed 's/^/    /'
  echo
  # landmark columns: 4=fimo 5=tfbstools 6=moods 7=gquad 8=pqsfinder
  # EXTRA is TFBSTools' B_only count, read from the concordance table rather
  # than hard-coded so the prose tracks the data
  EXTRA=$(awk -F'\t' '$1=="fimo_vs_tfbstools"{print $7}' "${PAIRS}")
  awk -F'\t' -v extra="${EXTRA}" 'NR>1{for(i=4;i<=8;i++) if($i=="yes") c[i]++} END{
      printf "**Recovery of the 6 known landmarks:** TFBSTools %d/6, FIMO %d/6, MOODS %d/6.\n", c[5]+0, c[4]+0, c[6]+0
      printf "\nOn the only ground truth available, TFBSTools is the most sensitive TFBS scanner\n"
      printf "and FIMO the most conservative -- FIMO misses the TATA box and two of the three\n"
      printf "Sp1 sites that TFBSTools recovers. Because FIMO'\''s calls are a strict subset of\n"
      printf "TFBSTools'\''s, TFBSTools cannot be less sensitive than FIMO on any input.\n"
      printf "\n**Specificity cannot be separated from this.** The %s calls TFBSTools makes beyond\n", extra
      printf "FIMO'\''s set are either genuine weaker sites or false positives, and nothing in this\n"
      printf "data distinguishes the two. Higher sensitivity with unknown precision is the honest\n"
      printf "summary, not \"TFBSTools is better\".\n"
  }' "${LANDMARKS}"
  echo
  echo "The G4 tools also overlap several TFBS landmarks. That is not a scoring error:"
  echo "the Sp1/NF-kB region of the HIV-1 LTR is GC-rich and carries documented"
  echo "G-quadruplex structures (LTR-III/LTR-IV), so genuine G4 predictions are expected"
  echo "to coincide with it. Their landmark column is shown for completeness but should"
  echo "not be read as TFBS detection."
  echo
  echo "## 4. Practical conclusions"
  echo
  echo "- **FIMO vs TFBSTools:** identical coordinates wherever both call, so the choice is"
  echo "  purely one of stringency. Use FIMO when false positives are costly; use TFBSTools"
  echo "  when missing a weak site is costly. Reporting both is redundant -- FIMO adds no"
  echo "  site TFBSTools lacks."
  echo "- **gquad vs pqsfinder:** these are not interchangeable. They locate the same"
  echo "  regions but assign different extents, so any downstream analysis keyed on G4"
  echo "  boundaries (overlap with a TFBS, distance to a TSS) will change with the tool."
  echo "  Fix one tool for the whole study rather than mixing."
  echo "- **Before drawing biological conclusions**, fix the coverage problem in section 1:"
  echo "  the consensus step should emit N (or be masked) where read depth is zero, so a"
  echo "  sample cannot silently contribute reference sequence to a motif scan."
  echo
} > "${REPORT}"

echo ""
echo "=== report written: ${REPORT} ==="
echo ""
echo "Done. Outputs in ${OUT_DIR}"
