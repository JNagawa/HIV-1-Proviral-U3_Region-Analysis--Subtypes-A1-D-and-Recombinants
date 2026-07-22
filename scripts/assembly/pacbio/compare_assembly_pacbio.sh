#!/bin/bash
# Assembly step of the PacBio (HIV-SMRTcap) harness: reference-guided
# minimap2->HXB2 + bcftools consensus  vs  de novo hifiasm, on the
# per-sample proviral-core reads from the extraction step. This mirrors the
# Illumina assembly comparison (bwa_consensus vs SPAdes): a reference-guided
# consensus that is robust at low coverage, against a reference-free de novo
# assembler that can capture divergent/structural variation the reference
# would mask. The tools review names minimap2 as the long-read mapper but no
# long-read de novo assembler, so hifiasm (the PacBio-HiFi standard) is the
# de novo counterpart chosen here.
# Usage: ./compare_assembly_pacbio.sh   (via sbatch scripts/utils/run_comparison_step.slurm.sh)
set -uo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
source "${REPO_ROOT}/scripts/common/lib_compare.sh"

REF_FASTA="${REPO_ROOT}/data/reference/K03455.1.fasta"
PROVIRUS_DIR="${REPO_ROOT}/results/proviral_extraction/pacbio"
RESULTS_DIR="${REPO_ROOT}/results/assembly/pacbio"
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
mkdir -p "${RESULTS_DIR}"
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

if [ ! -s "${REF_FASTA}" ]; then echo "ERROR: HXB2 reference ${REF_FASTA} not found." >&2; exit 1; fi
THREADS="${THREADS:-4}"
export THREADS

# Reference-guided: map proviral-core reads to HXB2 (HiFi preset) and call a
# consensus. Same shape as the Illumina bwa_consensus function.
run_minimap2_consensus() {
    local srr="$1" reads="$2" ref="$3" outdir="$4"
    mkdir -p "${outdir}"
    local bam="${outdir}/${srr}.sorted.bam" vcf="${outdir}/${srr}.vcf.gz" consensus="${outdir}/${srr}_consensus.fasta"
    minimap2 -a -x map-hifi -t "${THREADS:-4}" "${ref}" "${reads}" 2>"${outdir}/${srr}_minimap2.log" \
        | samtools view -b - | samtools sort -o "${bam}" || exit 1
    samtools index "${bam}" || exit 1
    bcftools mpileup -Ou -f "${ref}" "${bam}" 2>"${outdir}/${srr}_bcftools.log" \
        | bcftools call -c --ploidy 1 -Oz -o "${vcf}" 2>>"${outdir}/${srr}_bcftools.log" || exit 1
    tabix -p vcf "${vcf}" || exit 1
    cat "${ref}" | bcftools consensus "${vcf}" > "${consensus}" 2>>"${outdir}/${srr}_bcftools.log" || exit 1
    sed -i "1s/.*/>${srr}/" "${consensus}"
}
export -f run_minimap2_consensus

# De novo: hifiasm assembles the proviral-core reads; take the primary contig
# best matching HXB2 (by BLAST) as the assembly, else the longest -- same
# best-contig selection as the Illumina SPAdes function.
run_hifiasm() {
    local srr="$1" reads="$2" ref="$3" outdir="$4"
    mkdir -p "${outdir}"
    local prefix="${outdir}/${srr}"
    hifiasm -o "${prefix}" -t "${THREADS:-4}" "${reads}" > "${outdir}/${srr}_hifiasm.log" 2>&1 || exit 1
    local gfa="${prefix}.bp.p_ctg.gfa"
    [ -s "${gfa}" ] || gfa="${prefix}.p_ctg.gfa"
    [ -s "${gfa}" ] || { echo "ERROR: hifiasm produced no primary-contig GFA for ${srr}" >&2; exit 1; }
    local contigs="${outdir}/${srr}_contigs.fasta"
    awk '/^S/{print ">"$2"\n"$3}' "${gfa}" > "${contigs}"
    [ -s "${contigs}" ] || { echo "ERROR: no contigs parsed from ${gfa}" >&2; exit 1; }

    if command -v makeblastdb >/dev/null 2>&1 && command -v blastn >/dev/null 2>&1; then
        makeblastdb -in "${ref}" -dbtype nucl -out "${outdir}/${srr}_hxb2db" >/dev/null 2>&1
        local best
        best=$(blastn -query "${contigs}" -db "${outdir}/${srr}_hxb2db" -outfmt "6 qseqid length bitscore" 2>/dev/null \
            | sort -k3,3 -rn | head -1 | cut -f1)
        if [ -n "${best}" ]; then
            seqkit grep -n -p "${best}" "${contigs}" > "${outdir}/${srr}_consensus.fasta"
            sed -i "1s/.*/>${srr}/" "${outdir}/${srr}_consensus.fasta"
            return 0
        fi
    fi
    seqkit sort -l -r "${contigs}" 2>/dev/null | seqkit head -n 1 > "${outdir}/${srr}_consensus.fasta"
    sed -i "1s/.*/>${srr}/" "${outdir}/${srr}_consensus.fasta"
}
export -f run_hifiasm

for SRR in $(subset_accessions pacbio "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do
    READS="${PROVIRUS_DIR}/${SRR}.provirus.fasta"
    if [ ! -s "${READS}" ]; then
        echo "WARNING: no proviral cores for ${SRR} (run extract_provirus_pacbio.sh first), skipping." >&2
        continue
    fi

    for TOOL in minimap2_consensus hifiasm; do
        if [ "${TOOL}" = "hifiasm" ] && ! command -v hifiasm >/dev/null 2>&1; then
            echo "NOTE: hifiasm not installed, skipping (see HIV_U3analysis_env.yml)." >&2
            continue
        fi
        OUTDIR="${RESULTS_DIR}/${TOOL}_out"
        TIMELOG="${RESULTS_DIR}/${TOOL}_${SRR}.time"
        LOG="${RESULTS_DIR}/${TOOL}_${SRR}.log"
        CONSENSUS="${OUTDIR}/${SRR}_consensus.fasta"

        echo "=== ${TOOL} on ${SRR} ==="
        measure_and_run "${TIMELOG}" -- \
            bash -c 'run_'"${TOOL}"' "$@"' _ "${SRR}" "${READS}" "${REF_FASTA}" "${OUTDIR}" > "${LOG}" 2>&1
        EXIT_CODE=$?
        parse_time_metrics "${TIMELOG}"

        VALID=0; METRIC="n/a"
        if [ -s "${CONSENSUS}" ]; then
            LEN=$(seqkit stats -T "${CONSENSUS}" 2>/dev/null | tail -1 | cut -f5)
            N_PCT=$(seqkit fx2tab -n -g -B N "${CONSENSUS}" 2>/dev/null | awk -F'\t' '{print $NF}' | tail -1)
            if [ -n "${LEN}" ] && [ "${LEN}" -ge 8000 ] && [ "${LEN}" -le 10000 ]; then VALID=1; fi
            METRIC="${LEN}bp, ${N_PCT:-?}% N"
        fi
        append_summary_row "assembly_pacbio" "${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"
    done
done

echo "Done. See ${SUMMARY_TSV}"
