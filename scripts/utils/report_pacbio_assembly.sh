#!/bin/bash
# Generate a findings report for the PacBio assembly stage.
#
# Why this exists: summary.tsv carries one row per (tool, sample) with a single
# free-text key_metric, which is enough to see whether a run succeeded but not
# enough to answer the question this project actually turns on -- for any given
# region, was it sequenced, or is it simply absent? Those two are what a
# consensus FASTA cannot distinguish on its own, and conflating them is exactly
# how HXB2's own motif hits were once reported as sample results.
#
# So this report is written in REFERENCE coordinates, read off the mask BED and
# the VCF rather than off the consensus FASTA. The consensus has indels applied,
# so its coordinates drift from the reference; the mask and VCF do not.
#
# Everything here is derived from files the assembly stage already writes:
#   <arm>_out/<sample>.uncovered.bed      unsequenced intervals (the mask)
#   <arm>_out/<sample>.read_spans.tsv     per-read reference spans
#   <arm>_out/<sample>.vcf.gz             called variants
#   <arm>_out/<sample>_reference.txt      which reference this arm chose
#   <arm>_out/<sample>_consensus.fasta    the consensus itself
#   summary.tsv                           the per-run timing/validity table
#
# Usage: report_pacbio_assembly.sh [output.md]
set -euo pipefail

# locate the repo root the same way the comparison harnesses do
REPO_ROOT=""
for CAND in "${SLURM_SUBMIT_DIR:-}" \
            "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)" \
            "$(pwd)"; do
    if [ -n "${CAND}" ] && [ -d "${CAND}/scripts" ] && [ -d "${CAND}/results" ]; then
        REPO_ROOT="${CAND}"; break
    fi
done
[ -n "${REPO_ROOT}" ] || { echo "ERROR: cannot locate repo root." >&2; exit 1; }

ASM_DIR="${REPO_ROOT}/results/assembly/pacbio"
SUMMARY="${ASM_DIR}/summary.tsv"
# HXB2's GenBank record, cached by the U3 extraction step. Parsed here rather
# than hardcoding U3 boundaries so this report and the extraction can never
# disagree about where U3 is.
GB="${REPO_ROOT}/data/reference/K03455.1.gb"
OUT="${1:-${ASM_DIR}/assembly_report.md}"
# reference-guided arms only; the de novo arm has no BAM/VCF/mask to report on
ARMS="${ASSEMBLY_ARMS:-minimap2_consensus minimap2_bestref}"

[ -s "${SUMMARY}" ] || { echo "ERROR: ${SUMMARY} not found -- run the assembly stage first." >&2; exit 1; }

# sample ids, taken from the summary table so the report covers exactly what ran
SAMPLES=$(awk -F'\t' 'NR>1 && $3!="" {print $3}' "${SUMMARY}" | sort -u)
[ -n "${SAMPLES}" ] || { echo "ERROR: no samples found in ${SUMMARY}." >&2; exit 1; }

# --- U3 window boundaries, read off HXB2's own annotation ---------------------
# Same note-matching logic as extract_u3_by_hxb2_anchor.sh: HXB2 annotates four
# relevant repeat_regions (the 5'/3' LTR pair and the 5'/3' R-repeat pair), and
# U3 is the LTR up to but not including its own R region.
U3_COORDS=""
if [ -s "${GB}" ]; then
    U3_COORDS=$(awk '
    /^[[:space:]]*\// {
        if ($0 ~ /\/note=/ && pending_type == "repeat_region") {
            note = $0
            sub(/.*\/note="/, "", note); sub(/".*/, "", note)
            notelow = tolower(note)
            is5 = (index(notelow, "5") > 0); is3 = (index(notelow, "3") > 0)
            if (index(notelow, "ltr") > 0) {
                if (is5 && l5 == "") l5 = pending_loc
                else if (is3 && l3 == "") l3 = pending_loc
            } else if (index(notelow, "r repeat") > 0) {
                if (is5 && r5 == "") r5 = pending_loc
                else if (is3 && r3 == "") r3 = pending_loc
            }
        }
        next
    }
    NF >= 2 && $1 ~ /^[A-Za-z_]+$/ { pending_type = $1; pending_loc = $2; next }
    END {
        if (l5 == "" || r5 == "" || l3 == "" || r3 == "") exit 1
        split(l5, a, "\\.\\."); split(r5, b, "\\.\\.")
        split(l3, c, "\\.\\."); split(r3, d, "\\.\\.")
        # emit 1-based inclusive windows: U3 = [LTR_start, R_start-1]
        print a[1], b[1]-1, c[1], d[1]-1
    }' "${GB}" 2>/dev/null || true)
fi
if [ -n "${U3_COORDS}" ]; then
    read -r U5S U5E U3S U3E <<EOF
${U3_COORDS}
EOF
else
    echo "WARNING: could not parse U3 windows from ${GB}; U3 section will be omitted." >&2
fi

# --- helpers -----------------------------------------------------------------
# Reference length, taken from the BAM's own @SQ header. Deliberately not read
# off the mask BED: the BED's last interval only reaches the reference end when
# the sample happens to be uncovered there, so for a sample whose 3' end IS
# covered (124_4) that would report a reference a fraction of its real length,
# silently corrupting every percentage and the interior/edge test below.
ref_length() {
    local bam="$1"
    [ -s "${bam}" ] || { echo 0; return; }
    samtools view -H "${bam}" 2>/dev/null \
        | awk '/^@SQ/ { for (i = 1; i <= NF; i++) if ($i ~ /^LN:/) { sub(/^LN:/, "", $i); print $i; exit } }'
}

# count reference positions in [w_start, w_end] that fall inside a mask BED
# (BED is 0-based half-open, so interval [s,e) covers 1-based s+1 .. e)
masked_in_window() {
    local bed="$1" ws="$2" we="$3"
    [ -s "${bed}" ] || { echo 0; return; }
    awk -v WS="${ws}" -v WE="${we}" '
        { s = $2 + 1; e = $3
          lo = (s > WS ? s : WS); hi = (e < WE ? e : WE)
          if (hi >= lo) n += hi - lo + 1 }
        END { print n+0 }' "${bed}"
}

# count called variants of a given class inside a window, in reference coords
variants_in_window() {
    local vcf="$1" ws="$2" we="$3" kind="$4"
    [ -s "${vcf}" ] || { echo 0; return; }
    bcftools view -H -v "${kind}" "${vcf}" 2>/dev/null \
        | awk -F'\t' -v WS="${ws}" -v WE="${we}" '$2 >= WS && $2 <= WE {n++} END{print n+0}'
}

# split called indels into insertion / deletion counts and total deleted bases
indel_breakdown() {
    local vcf="$1"
    [ -s "${vcf}" ] || { echo "0 0 0"; return; }
    bcftools view -H -v indels "${vcf}" 2>/dev/null \
        | awk -F'\t' '{
              rl = length($4); al = length($5)
              if (rl > al) { del++; delbp += rl - al } else if (al > rl) ins++
          } END { print ins+0, del+0, delbp+0 }'
}

# --- report ------------------------------------------------------------------
mkdir -p "$(dirname "${OUT}")"
{
printf '# PacBio assembly -- findings report\n\n'
printf '**Generated:** %s  ·  auto-generated by `scripts/utils/report_pacbio_assembly.sh`\n\n' "$(date +%Y-%m-%d)"
printf 'All figures are in **reference coordinates**, read from the mask BED and the VCF\n'
printf 'rather than from the consensus FASTA (the consensus has indels applied, so its\n'
printf 'own coordinates drift from the reference).\n\n'
printf 'The distinction this report exists to preserve:\n\n'
printf '%s **unsequenced** -- no read spans the position. Masked to `N`. Nothing can be\n' '-'
printf '  concluded about it, in either direction.\n'
printf '%s **deleted** -- reads span the position and agree the base is absent. A called\n' '-'
printf '  variant, not an `N`. This is a positive finding.\n\n'
printf 'A motif that is absent over an unsequenced stretch is **not evidence of absence**.\n\n'
printf '%s\n\n' '---'

# ============================ per-arm tables ==============================
for ARM in ${ARMS}; do
    ARMD="${ASM_DIR}/${ARM}_out"
    [ -d "${ARMD}" ] || continue
    printf '## Arm: `%s`\n\n' "${ARM}"
    printf '| sample | reference | reads | genome spanned | unsequenced | SNV | ins | del | del bp | consensus | called | %%N |\n'
    printf '|---|---|---|---|---|---|---|---|---|---|---|---|\n'
    for S in ${SAMPLES}; do
        BAM="${ARMD}/${S}.sorted.bam"; BED="${ARMD}/${S}.uncovered.bed"
        VCF="${ARMD}/${S}.vcf.gz";     CONS="${ARMD}/${S}_consensus.fasta"
        [ -s "${CONS}" ] || continue

        REFID=$(head -1 "${ARMD}/${S}_reference.txt" 2>/dev/null || echo "?")
        NREADS=$(samtools view -c -F 0x904 "${BAM}" 2>/dev/null || echo "?")
        # reference length, from the BAM header (see ref_length above)
        REFLEN=$(ref_length "${BAM}")
        MASKBP=$(awk '{n += $3 - $2} END{print n+0}' "${BED}" 2>/dev/null || echo 0)
        SNV=$(bcftools view -H -v snps "${VCF}" 2>/dev/null | wc -l || echo 0)
        read -r INS DEL DELBP <<EOF
$(indel_breakdown "${VCF}")
EOF
        # consensus length and N count, straight off the FASTA
        read -r CLEN CN <<EOF
$(awk '!/^>/{s = s $0} END{n = gsub(/[Nn]/, "N", s); print length(s), n+0}' "${CONS}")
EOF
        CALLED=$((CLEN - CN))
        SPANPCT=$(awk -v m="${MASKBP}" -v l="${REFLEN:-0}" 'BEGIN{ if (l>0) printf "%.1f%%", 100*(l-m)/l; else printf "?" }')
        NPCT=$(awk -v n="${CN}" -v l="${CLEN}" 'BEGIN{ if (l>0) printf "%.2f", 100*n/l; else printf "?" }')
        printf '| %s | %s | %s | %s | %s bp | %s | %s | %s | %s | %s bp | %s | %s |\n' \
            "${S}" "${REFID}" "${NREADS}" "${SPANPCT}" "${MASKBP}" \
            "${SNV}" "${INS}" "${DEL}" "${DELBP}" "${CLEN}" "${CALLED}" "${NPCT}"
    done
    printf '\n'

    # ---------------------- amplicon structure ----------------------
    printf '### Sequenced regions (amplicon structure)\n\n'
    printf 'Reads clustered by span endpoints (both ends agreeing within 50bp), NOT merged\n'
    printf 'by overlap. Reads that stack into a few discrete intervals rather than\n'
    printf 'tiling the reference indicate an amplicon library, which is what makes\n'
    printf 'overlap-based de novo assembly structurally impossible for some samples.\n\n'
    printf '| sample | distinct regions | regions (reference coords, read count) |\n'
    printf '|---|---|---|\n'
    for S in ${SAMPLES}; do
        SPANS="${ARMD}/${S}.read_spans.tsv"
        # An absent spans file means these outputs predate the span-based mask.
        # Say so rather than leaving a silently empty table, which would read as
        # "no amplicon structure" instead of "not measured".
        if [ ! -s "${SPANS}" ]; then
            printf '| %s | n/a | no `read_spans.tsv` -- outputs predate the span-based mask; re-run assembly |\n' "${S}"
            continue
        fi
        # Cluster reads by their span ENDPOINTS, not by merging overlaps. Merging
        # overlapping intervals collapses distinct amplicons whose ranges happen to
        # touch -- 124_4's three groups (814-8002, 1963-9719, 8986-9719) all overlap,
        # so merging reported a single 814-9719 region and made an amplicon library
        # look like continuous tiling, inverting the finding this table exists to
        # show. Two reads belong to the same amplicon when BOTH endpoints agree to
        # within TOL, which tolerates ragged ends without merging separate products.
        CLUSTERS=$(sort -k1,1n -k2,2n "${SPANS}" | awk -v TOL=50 '
            function flush() { if (c > 0) printf "%d-%d (x%d); ", smin, emax, c }
            NR == 1 { sref = $1; eref = $2; smin = $1; emax = $2; c = 1; next }
            ($1 - sref <= TOL) && (($2 - eref <= TOL) && (eref - $2 <= TOL)) {
                c++
                if ($1 < smin) smin = $1
                if ($2 > emax) emax = $2
                next
            }
            { flush(); sref = $1; eref = $2; smin = $1; emax = $2; c = 1 }
            END { if (c > 0) printf "%d-%d (x%d)", smin, emax, c }')
        NCLUST=$(printf '%s' "${CLUSTERS}" | awk -F';' '{print NF}')
        printf '| %s | %s | %s |\n' "${S}" "${NCLUST}" "${CLUSTERS}"
    done
    printf '\n'

    # ---------------------- masking audit ----------------------
    printf '### Masking audit\n\n'
    printf 'The mask must contain only stretches no read reaches. An **interior** masked\n'
    printf 'interval -- one with sequenced bases on both sides -- means a position that\n'
    printf 'reads do span was labelled unsequenced. That regression previously turned\n'
    printf 'unanimous deletions into `N`, so it is checked on every run.\n\n'
    printf '| sample | mask intervals | interior intervals | interior bp | status |\n'
    printf '|---|---|---|---|---|\n'
    for S in ${SAMPLES}; do
        BED="${ARMD}/${S}.uncovered.bed"
        [ -s "${BED}" ] || continue
        RL=$(ref_length "${ARMD}/${S}.sorted.bam")
        awk -v s="${S}" -v L="${RL}" -v OFS=' ' '
            { n++
              # an interval touching either genome end is an amplicon edge gap,
              # not an interior hole
              if ($2 != 0 && $3 != L) { ic++; ibp += $3 - $2 } }
            END {
                status = (ic+0 == 0) ? "OK" : "**CHECK -- interior masking present**"
                printf "| %s | %d | %d | %d | %s |\n", s, n+0, ic+0, ibp+0, status
            }' "${BED}"
    done
    printf '\n'

    # ---------------------- U3 assessment ----------------------
    if [ -n "${U3_COORDS:-}" ]; then
        printf '### U3 coverage (the regions of interest)\n\n'
        printf 'U3 windows taken from HXB2 annotation in `%s`:\n' "data/reference/K03455.1.gb"
        printf '5%s U3 = %s-%s, 3%s U3 = %s-%s (1-based inclusive).\n\n' "'" "${U5S}" "${U5E}" "'" "${U3S}" "${U3E}"
        printf 'Both LTRs are identical in an integrated provirus, so **either** copy is a\n'
        printf 'valid U3 source -- but an amplicon library may only reach one of them.\n\n'
        printf 'Windows are HXB2 coordinates, and are only reported for samples mapped against\n'
        printf 'HXB2. U3 position varies between subtypes (Mbondji-Wonje 2018 reports up to 25%%\n'
        printf 'inter-strain dissimilarity in U3 while R stays conserved), and no per-subtype\n'
        printf 'coordinate system exists to substitute, so a window projected onto a non-HXB2\n'
        printf 'reference would be meaningless. Subtype-aware U3 boundaries come from the\n'
        printf 'alignment-anchored liftover downstream, not from this table.\n\n'
        printf "| sample | 5' U3 sequenced | 3' U3 sequenced | del called in 5' | del called in 3' | usable copy |\n"
        printf '|---|---|---|---|---|---|\n'
        for S in ${SAMPLES}; do
            BED="${ARMD}/${S}.uncovered.bed"; VCF="${ARMD}/${S}.vcf.gz"
            [ -s "${BED}" ] || continue
            REFID=$(head -1 "${ARMD}/${S}_reference.txt" 2>/dev/null || echo "?")
            # U3 windows are HXB2 coordinates; they are meaningless against any
            # other reference, so say so rather than printing a wrong number
            case "${REFID}" in
                HXB2|*K03455*) ;;
                *) printf '| %s | n/a | n/a | n/a | n/a | reference is %s, not HXB2 -- U3 window undefined in this frame |\n' "${S}" "${REFID}"; continue ;;
            esac
            L5=$((U5E - U5S + 1)); L3=$((U3E - U3S + 1))
            M5=$(masked_in_window "${BED}" "${U5S}" "${U5E}")
            M3=$(masked_in_window "${BED}" "${U3S}" "${U3E}")
            S5=$((L5 - M5)); S3=$((L3 - M3))
            D5=$(variants_in_window "${VCF}" "${U5S}" "${U5E}" indels)
            D3=$(variants_in_window "${VCF}" "${U3S}" "${U3E}" indels)
            P5=$(awk -v a="${S5}" -v b="${L5}" 'BEGIN{printf "%.0f", 100*a/b}')
            P3=$(awk -v a="${S3}" -v b="${L3}" 'BEGIN{printf "%.0f", 100*a/b}')
            # the extraction downstream prefers whichever copy carries more
            # sequenced base, so state the same verdict here
            if [ "${S5}" -eq 0 ] && [ "${S3}" -eq 0 ]; then
                USABLE="**neither -- exclude from motif analysis**"
            elif [ "${S3}" -gt "${S5}" ]; then
                USABLE="3' copy"
            else
                USABLE="5' copy"
            fi
            printf "| %s | %s/%s bp (%s%%) | %s/%s bp (%s%%) | %s | %s | %s |\n" \
                "${S}" "${S5}" "${L5}" "${P5}" "${S3}" "${L3}" "${P3}" "${D5}" "${D3}" "${USABLE}"
        done
        printf '\n'
    fi
done

# ============================ cross-arm comparison ==========================
printf '%s\n\n## ' '---'; printf 'Does subtype-matching reduce divergence?\n\n'
printf 'Total called variants against each arm reference. Fewer variants means the\n'
printf 'reference is closer to the sample, so less of the consensus is inherited from\n'
printf 'a reference backbone.\n\n'
printf '| sample |'
for ARM in ${ARMS}; do printf ' %s (ref, SNV+indel) |' "${ARM}"; done
printf '\n|---|'
for ARM in ${ARMS}; do printf '%s' '---|'; done
printf '\n'
for S in ${SAMPLES}; do
    printf '| %s |' "${S}"
    for ARM in ${ARMS}; do
        VCF="${ASM_DIR}/${ARM}_out/${S}.vcf.gz"
        if [ -s "${VCF}" ]; then
            RID=$(head -1 "${ASM_DIR}/${ARM}_out/${S}_reference.txt" 2>/dev/null || echo "?")
            NS=$(bcftools view -H -v snps "${VCF}" 2>/dev/null | wc -l)
            NI=$(bcftools view -H -v indels "${VCF}" 2>/dev/null | wc -l)
            printf ' %s: %s |' "${RID}" "$((NS + NI))"
        else
            printf ' n/a |'
        fi
    done
    printf '\n'
done
printf '\n'

# ============================ warnings ==========================
printf '%s\n\n## ' '---'; printf 'Flags requiring attention\n\n'
FLAGGED=0
for ARM in ${ARMS}; do
    ARMD="${ASM_DIR}/${ARM}_out"
    [ -d "${ARMD}" ] || continue
    for S in ${SAMPLES}; do
        BED="${ARMD}/${S}.uncovered.bed"
        [ -s "${BED}" ] || continue
        RL=$(ref_length "${ARMD}/${S}.sorted.bam")
        # interior masking is a pipeline regression, not a data property
        IC=$(awk -v L="${RL}" '$2 != 0 && $3 != L {n++} END{print n+0}' "${BED}")
        if [ "${IC}" -gt 0 ]; then
            printf '%s `%s`/%s: **%s interior masked interval(s)** -- reads span these positions,\n' '-' "${ARM}" "${S}" "${IC}"
            printf '  so labelling them unsequenced is wrong. Check that the mask is still built\n'
            printf '  from read spans and that indel calling is enabled.\n'
            FLAGGED=1
        fi
        if [ -n "${U3_COORDS:-}" ]; then
            REFID=$(head -1 "${ARMD}/${S}_reference.txt" 2>/dev/null || echo "?")
            case "${REFID}" in
                HXB2|*K03455*)
                    M5=$(masked_in_window "${BED}" "${U5S}" "${U5E}")
                    M3=$(masked_in_window "${BED}" "${U3S}" "${U3E}")
                    L5=$((U5E - U5S + 1)); L3=$((U3E - U3S + 1))
                    if [ "${M5}" -eq "${L5}" ] && [ "${M3}" -eq "${L3}" ]; then
                        printf '%s `%s`/%s: **no U3 sequence at either LTR.** Any motif result for this\n' '-' "${ARM}" "${S}"
                        printf '  sample would come entirely from the reference. Exclude it from motif\n'
                        printf '  analysis rather than reporting reference motifs as sample findings.\n'
                        FLAGGED=1
                    fi
                    ;;
            esac
        fi
    done
done
[ "${FLAGGED}" -eq 0 ] && printf 'None.\n'
printf '\n'
printf '%s\n\n' '---'
printf 'Per-run timings and validity: `results/assembly/pacbio/summary.tsv`\n'
} > "${OUT}"

echo "Wrote assembly findings report to ${OUT}"
