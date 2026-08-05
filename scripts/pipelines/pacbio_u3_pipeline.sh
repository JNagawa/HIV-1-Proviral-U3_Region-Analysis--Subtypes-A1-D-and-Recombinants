#!/bin/bash
#SBATCH --job-name=pb_u3_e2e
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

##################################################################################
# End-to-end HIV-1 proviral U3 pipeline for PacBio HIV-SMRTcap (HiFi) data.
#
# This is the PRODUCTION pipeline: it runs only the tools chosen by the
# comparison harness, in one pass, rather than running every candidate. The
# harness (scripts/<stage>/pacbio/compare_*.sh) is what justified these choices;
# this script is what you run once they are settled.
#
#   extraction   minimap2 mask + N-strip
#   assembly     minimap2 map-hifi -> bcftools mpileup/call/consensus
#   MSA          MAFFT
#   U3           HXB2-anchored alignment liftover, per-record LTR choice
#   motifs       FIMO
#   filtering    HIV-Intact (per-sample subtype) + Poplars
#   subtyping    jpHMM + IQ-TREE2   (OPTIONAL -- jpHMM runs >1h per sample)
#
# Two configuration choices below are load-bearing and were established by
# measurement on 2026-08-05; changing them silently breaks correctness:
#
#   1. `bcftools mpileup -X pacbio-ccs`. With generic defaults this step called
#      ZERO indels in all four samples even where every read agreed, because
#      mpileup's indel-candidate detection is tuned for short reads. Switching
#      the caller -c -> -m alone did NOT fix it.
#   2. The N-mask is built from read alignment SPANS, not `samtools depth`. depth
#      counts bases, so a position deleted in every read reports depth 0 and is
#      indistinguishable from a position no read reached. Masking on depth wrote
#      N over the dataset's best-supported deletions.
#
# Usage:
#   sbatch scripts/pipelines/pacbio_u3_pipeline.sh
#   SAMPLES="124_4 203_3" RUN_SUBTYPING=1 sbatch scripts/pipelines/pacbio_u3_pipeline.sh
#
# Environment overrides:
#   SAMPLES         space-separated sample ids   (default: the pacbio subset TSV)
#   REFERENCE       reference FASTA              (default: HXB2 K03455.1)
#   OUTROOT         output directory              (default: results/pipeline/pacbio)
#   RUN_SUBTYPING   1 to include jpHMM/IQ-TREE    (default: 0, it is very slow)
#   THREADS         thread count                  (default: Slurm allocation or 4)
#
# -u errors on unset vars, pipefail fails a pipe if any stage fails. NOT -e:
# a single sample failing must not abort the whole run, so failures are counted
# and reported at the end instead.
set -uo pipefail

##################################################################################
#  Environment
##################################################################################
CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"
if [ -f "${CONDA_SH}" ]; then source "${CONDA_SH}"
else source "$(conda info --base)/etc/profile.d/conda.sh"
fi
conda activate HIV_U3analysis

REPO_ROOT=""
for CAND in "${SLURM_SUBMIT_DIR:-}" "$(cd "$(dirname "$0")/../.." 2>/dev/null && pwd)" "$(pwd)"; do
    if [ -n "${CAND}" ] && [ -f "${CAND}/scripts/common/lib_compare.sh" ]; then
        REPO_ROOT="${CAND}"; break
    fi
done
[ -n "${REPO_ROOT}" ] || { echo "ERROR: cannot locate repo root." >&2; exit 1; }
cd "${REPO_ROOT}" || exit 1
source "${REPO_ROOT}/scripts/common/lib_compare.sh"

THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
REFERENCE="${REFERENCE:-${REPO_ROOT}/data/reference/K03455.1.fasta}"
OUTROOT="${OUTROOT:-${REPO_ROOT}/results/pipeline/pacbio}"
RUN_SUBTYPING="${RUN_SUBTYPING:-0}"
GB_CACHE="${REPO_ROOT}/data/reference/K03455.1.gb"
# reads that have already been host-N-masked upstream
READS_DIR="${REPO_ROOT}/data/raw/pacbio/local_masked"

if [ -z "${SAMPLES:-}" ]; then
    SAMPLES=$(subset_accessions pacbio "${REPO_ROOT}/scripts/common/subset_samples.tsv" | tr '\n' ' ')
fi
[ -n "${SAMPLES}" ] || { echo "ERROR: no samples resolved." >&2; exit 1; }

mkdir -p "${OUTROOT}"/{extraction,assembly,msa,u3,motifs,filtering,subtyping}
STATUS_TSV="${OUTROOT}/pipeline_status.tsv"
printf 'stage\tsample\texit\tdetail\n' > "${STATUS_TSV}"
N_FAIL=0

# Record one stage/sample outcome and keep a running failure count.
note() {
    local stage="$1" sample="$2" code="$3" detail="$4"
    printf '%s\t%s\t%s\t%s\n' "${stage}" "${sample}" "${code}" "${detail}" >> "${STATUS_TSV}"
    [ "${code}" -eq 0 ] || N_FAIL=$((N_FAIL+1))
    printf '  [%s] %-8s exit=%s  %s\n' "${stage}" "${sample}" "${code}" "${detail}"
}

echo "=== PacBio U3 end-to-end pipeline ==="
echo "samples : ${SAMPLES}"
echo "threads : ${THREADS}"
echo "output  : ${OUTROOT}"
echo "subtyping: ${RUN_SUBTYPING} (0 = skipped; jpHMM exceeds an hour per sample)"

##################################################################################
#  Stage 1 -- proviral extraction
#
#  The input reads are host-N-masked, so the proviral core is the longest run of
#  non-N sequence. Reads are also oriented forward, because a proviral core can
#  arrive reverse-complemented and every later stage assumes one orientation.
##################################################################################
echo "=== [1/7] proviral extraction ==="
for S in ${SAMPLES}; do
    IN=$(ls "${READS_DIR}/${S}"*.fasta "${READS_DIR}/${S}"*.fa "${READS_DIR}/${S}"*.fastq.gz 2>/dev/null | head -1)
    OUT="${OUTROOT}/extraction/${S}.provirus.fasta"
    if [ -z "${IN}" ]; then note extraction "${S}" 1 "no reads found under ${READS_DIR}"; continue; fi
    # keep the longest non-N stretch per read, then drop anything under 500bp
    awk -v minlen=500 '
        /^>/ { if (id != "") emit(); id = substr($0,2); seq = ""; next }
        { seq = seq $0 }
        function emit(   n, parts, i, best) {
            n = split(toupper(seq), parts, /N+/); best = ""
            for (i = 1; i <= n; i++) if (length(parts[i]) > length(best)) best = parts[i]
            if (length(best) >= minlen) { print ">" id; print best }
        }
        END { if (id != "") emit() }
    ' "${IN}" > "${OUT}" 2>/dev/null
    NREADS=$(grep -c '^>' "${OUT}" 2>/dev/null || echo 0)
    if [ "${NREADS}" -gt 0 ]; then note extraction "${S}" 0 "${NREADS} proviral cores"
    else note extraction "${S}" 1 "no core >=500bp"; fi
done

##################################################################################
#  Stage 2 -- reference-guided assembly
##################################################################################
echo "=== [2/7] assembly (minimap2 + bcftools) ==="
for S in ${SAMPLES}; do
    READS="${OUTROOT}/extraction/${S}.provirus.fasta"
    D="${OUTROOT}/assembly"; BAM="${D}/${S}.sorted.bam"
    VCF="${D}/${S}.vcf.gz"; CONS="${D}/${S}_consensus.fasta"
    MASK="${D}/${S}.uncovered.bed"; SPANS="${D}/${S}.read_spans.tsv"
    [ -s "${READS}" ] || { note assembly "${S}" 1 "no extracted reads"; continue; }

    minimap2 -a -x map-hifi -t "${THREADS}" "${REFERENCE}" "${READS}" 2>"${D}/${S}_minimap2.log" \
        | samtools view -b - | samtools sort -o "${BAM}" 2>/dev/null
    samtools index "${BAM}" 2>/dev/null

    # -X pacbio-ccs is required for indel calling; see the header note
    bcftools mpileup -X pacbio-ccs -Ou -f "${REFERENCE}" "${BAM}" 2>"${D}/${S}_bcftools.log" \
        | bcftools call -m --ploidy 1 -Oz -o "${VCF}" 2>>"${D}/${S}_bcftools.log"
    tabix -f -p vcf "${VCF}" 2>/dev/null

    # Coverage from read alignment spans, NOT samtools depth; see the header note.
    # Reference span = POS .. POS + (reference-consuming CIGAR ops) - 1, where
    # M/D/N/=/X advance the reference and I/S/H/P do not.
    samtools view -F 0x904 "${BAM}" 2>/dev/null | awk -F'\t' '
        { cig = $6; if (cig == "*") next
          span = 0
          while (match(cig, /^[0-9]+[MIDNSHP=X]/)) {
              tok = substr(cig, RSTART, RLENGTH)
              n = substr(tok, 1, length(tok)-1) + 0
              op = substr(tok, length(tok))
              if (op ~ /^[MDN=X]$/) span += n
              cig = substr(cig, RSTART + RLENGTH)
          }
          if (span > 0) print $4, $4 + span - 1 }' > "${SPANS}"
    REFID=$(seqkit fx2tab -n -i "${REFERENCE}" | head -1 | cut -f1)
    REFLEN=$(seqkit fx2tab -n -l -i "${REFERENCE}" | head -1 | cut -f2)
    # complement the spanned intervals, so an empty BAM still yields a full mask
    awk -v L="${REFLEN}" -v CHR="${REFID}" -v OFS='\t' '
        { for (p = $1; p <= $2; p++) cov[p] = 1 }
        END { inrun = 0
              for (i = 1; i <= L; i++) {
                  if (!(i in cov)) { if (!inrun) { s = i; inrun = 1 } }
                  else if (inrun) { print CHR, s-1, i-1; inrun = 0 }
              }
              if (inrun) print CHR, s-1, L }' "${SPANS}" > "${MASK}"

    bcftools consensus -m "${MASK}" -f "${REFERENCE}" "${VCF}" > "${CONS}" 2>>"${D}/${S}_bcftools.log"
    sed -i "1s/.*/>${S}/" "${CONS}" 2>/dev/null

    if [ -s "${CONS}" ]; then
        read -r LEN NN <<EOF
$(awk '!/^>/{s = s $0} END{n = gsub(/[Nn]/, "N", s); print length(s), n+0}' "${CONS}")
EOF
        CALLED=$((LEN - NN))
        SNV=$(bcftools view -H -v snps "${VCF}" 2>/dev/null | wc -l)
        IND=$(bcftools view -H -v indels "${VCF}" 2>/dev/null | wc -l)
        # an interior masked interval means a covered position was called unsequenced
        INT=$(awk -v L="${REFLEN}" '$2 != 0 && $3 != L {n++} END{print n+0}' "${MASK}")
        note assembly "${S}" 0 "${LEN}bp, ${CALLED} called, ${SNV} SNV, ${IND} indel, interior-mask=${INT}"
        [ "${INT}" -eq 0 ] || echo "    WARNING: ${S} has ${INT} interior masked interval(s) -- mask may be misbuilt" >&2
    else
        note assembly "${S}" 1 "no consensus produced"
    fi
done

##################################################################################
#  Stage 3 -- MSA (MAFFT)
#
#  HXB2 is included as a record on purpose: its GenBank annotation is what the U3
#  liftover projects onto the samples, so it has to share the coordinate frame.
##################################################################################
echo "=== [3/7] MSA (MAFFT) ==="
COMBINED="${OUTROOT}/msa/combined_input.fasta"
ALIGNED="${OUTROOT}/msa/mafft_aligned.fasta"
cat "${REFERENCE}" "${OUTROOT}/assembly"/*_consensus.fasta > "${COMBINED}" 2>/dev/null
NREC=$(grep -c '^>' "${COMBINED}" 2>/dev/null || echo 0)
if [ "${NREC}" -ge 2 ]; then
    mafft --auto --thread "${THREADS}" "${COMBINED}" > "${ALIGNED}" 2>"${OUTROOT}/msa/mafft.log"
    COLS=$(awk '/^>/{if(n++)exit;next}{c += length($0)} END{print c+0}' "${ALIGNED}")
    note msa all "$?" "${NREC} records, ${COLS} columns"
else
    note msa all 1 "fewer than 2 records to align"
fi

##################################################################################
#  Stage 4 -- U3 extraction
##################################################################################
echo "=== [4/7] U3 extraction ==="
U3="${OUTROOT}/u3/U3_extracted.fasta"
if [ -s "${ALIGNED}" ]; then
    bash "${REPO_ROOT}/scripts/utils/extract_u3_by_hxb2_anchor.sh" \
        --alignment "${ALIGNED}" --hxb2-id K03455.1 --gb-cache "${GB_CACHE}" \
        --out-gapped "${OUTROOT}/u3/U3_aligned.fasta" --out "${U3}" \
        --warnings-log "${OUTROOT}/u3/warnings.log" > "${OUTROOT}/u3/extract.log" 2>&1
    EC=$?
    NU3=$(grep -c '^>' "${U3}" 2>/dev/null || echo 0)
    note u3 all "${EC}" "${NU3} U3 sequences; $(wc -l < "${OUTROOT}/u3/warnings.log" 2>/dev/null || echo 0) flagged"
    # a record with no sequenced base in either LTR carries no U3 evidence at all
    [ -s "${OUTROOT}/u3/warnings.log" ] && sed 's/^/    FLAG: /' "${OUTROOT}/u3/warnings.log"
else
    note u3 all 1 "no alignment"
fi

##################################################################################
#  Stage 5 -- motif mapping (FIMO)
##################################################################################
echo "=== [5/7] motif mapping (FIMO) ==="
MEME=$(ls "${REPO_ROOT}"/data/reference/jaspar/*.meme 2>/dev/null | head -1)
if [ -s "${U3}" ] && [ -n "${MEME}" ]; then
    fimo --oc "${OUTROOT}/motifs/fimo_out" --verbosity 1 "${MEME}" "${U3}" \
        > "${OUTROOT}/motifs/fimo.log" 2>&1
    EC=$?
    HITS=$(awk -F'\t' 'NR>1 && $3!="" && $1!~/^#/' "${OUTROOT}/motifs/fimo_out/fimo.tsv" 2>/dev/null | wc -l)
    note motifs all "${EC}" "${HITS} hits"
    # per-sample motif table, the actual deliverable of this stage
    awk -F'\t' 'NR>1 && $3!="" && $1!~/^#/ {printf "%s\t%s\t%s-%s\n",$3,$2,$4,$5}' \
        "${OUTROOT}/motifs/fimo_out/fimo.tsv" 2>/dev/null | sort -u > "${OUTROOT}/motifs/hits_by_sequence.tsv"
else
    note motifs all 1 "missing U3 FASTA or JASPAR .meme"
fi

##################################################################################
#  Stage 6 -- biological filtering
#
#  HIV-Intact scores against a subtype-specific reference, so each sequence is run
#  on its own subtype. Scoring a whole batch under one subtype reports HXB2 itself
#  as defective, which is how the previous configuration went unnoticed.
#  Poplars needs an ALIGNED FASTA -- it asserts equal record lengths.
##################################################################################
echo "=== [6/7] biological filtering ==="
FD="${OUTROOT}/filtering"
RUN_BF="${REPO_ROOT}/scripts/biological_filtering/illumina"
BESTREF_DIR="${REPO_ROOT}/results/assembly/pacbio/minimap2_bestref_out"
SUPPORTED="A1 A2 B C D F1 F2 G H HXB2"
printf 'sequence\tsubtype\tsource\n' > "${FD}/subtype_used.tsv"
for S in ${SAMPLES}; do
    CONS="${OUTROOT}/assembly/${S}_consensus.fasta"
    [ -s "${CONS}" ] || { note filtering "${S}" 1 "no consensus"; continue; }
    # subtype estimate comes from the harness's best-reference selection if present
    RAW=$(head -1 "${BESTREF_DIR}/${S}_reference.txt" 2>/dev/null)
    CAND="${RAW%%:*}"
    case " ${SUPPORTED} " in
        *" ${CAND} "*) SUB="${CAND}"; SRC="best-reference ${RAW}" ;;
        *) SUB=HXB2; SRC="no estimate (${RAW:-none}) -- defaulted to HXB2" ;;
    esac
    printf '%s\t%s\t%s\n' "${S}" "${SUB}" "${SRC}" >> "${FD}/subtype_used.tsv"
    bash "${RUN_BF}/run_hivintact.sh" "${CONS}" "${FD}/${S}" "${SUB}" >> "${FD}/hivintact.log" 2>&1
    NI=$(grep -c '^>' "${FD}/${S}/intact.fasta" 2>/dev/null || echo 0)
    ERRS=$(grep -oE '"error": *"[^"]*"' "${FD}/${S}/errors.json" 2>/dev/null \
           | sed 's/.*": *"//; s/"$//' | sort -u | tr '\n' ',' | sed 's/,$//')
    note filtering "${S}" 0 "subtype=${SUB}, $([ "${NI}" -gt 0 ] && echo intact || echo defective)${ERRS:+ [${ERRS}]}"
done
# HXB2 as positive control: if it is not called intact, this stage is untrustworthy
bash "${RUN_BF}/run_hivintact.sh" "${REFERENCE}" "${FD}/HXB2_control" HXB2 >> "${FD}/hivintact.log" 2>&1
CTRL=$(grep -c '^>' "${FD}/HXB2_control/intact.fasta" 2>/dev/null || echo 0)
if [ "${CTRL}" -gt 0 ]; then note filtering HXB2_ctrl 0 "control INTACT as expected"
else note filtering HXB2_ctrl 1 "control NOT intact -- stage output is not trustworthy"; fi
# Poplars over the alignment
if [ -s "${ALIGNED}" ]; then
    bash "${RUN_BF}/run_poplars.sh" "${ALIGNED}" "${FD}/poplars_out.tsv" > "${FD}/poplars.log" 2>&1
    note filtering poplars "$?" "$(wc -l < "${FD}/poplars_out.tsv" 2>/dev/null || echo 0) rows"
fi

##################################################################################
#  Stage 7 -- subtyping (optional)
##################################################################################
echo "=== [7/7] subtyping ==="
if [ "${RUN_SUBTYPING}" = "1" ]; then
    RUN_ST="${REPO_ROOT}/scripts/subtyping/illumina"
    for S in ${SAMPLES}; do
        CONS="${OUTROOT}/assembly/${S}_consensus.fasta"
        [ -s "${CONS}" ] || continue
        bash "${RUN_ST}/run_jphmm.sh" "${CONS}" "${OUTROOT}/subtyping/${S}" \
            > "${OUTROOT}/subtyping/jphmm_${S}.log" 2>&1
        note subtyping "${S}" "$?" "jpHMM complete"
    done
    if [ -s "${ALIGNED}" ]; then
        bash "${RUN_ST}/run_iqtree.sh" "${ALIGNED}" "${OUTROOT}/subtyping/iqtree_out" \
            > "${OUTROOT}/subtyping/iqtree.log" 2>&1
        note subtyping iqtree "$?" "tree written"
    fi
else
    echo "  skipped (set RUN_SUBTYPING=1 to include; jpHMM exceeds an hour per sample)"
    note subtyping all 0 "skipped by configuration"
fi

##################################################################################
#  Summary
##################################################################################
echo
echo "=== pipeline status ==="
column -t -s$'\t' "${STATUS_TSV}"
echo
if [ -s "${OUTROOT}/motifs/hits_by_sequence.tsv" ]; then
    echo "=== motif hits by sequence ==="
    column -t -s$'\t' "${OUTROOT}/motifs/hits_by_sequence.tsv"
    echo
fi
echo "Outputs under ${OUTROOT}"
echo "Status table: ${STATUS_TSV}"
if [ "${N_FAIL}" -gt 0 ]; then
    echo "Completed with ${N_FAIL} failed step(s) -- see the table above." >&2
    exit 1
fi
echo "All steps completed."
