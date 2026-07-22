#!/bin/bash
# Proviral extraction step of the PacBio (HIV-SMRTcap) harness -- the special
# step: recover the integrated provirus from HiFi reads that carry host-genome
# flanks, where host is N-masked and the provirus is in ACGT, e.g.
#     NNNNNACGGTTCTAGGGTTTCCACTANNNNN  ->  ACGGTTCTAGGGTTTCCACTA
# The virus integrates at a different host position in every read, so the
# flanks vary in length and sit at both ends. Extraction = strip the leading
# and trailing N-runs (scripts/utils/extract_provirus_strip_hostN.sh).
#
# Two input paths:
#   * Path A (PROVIRUS_INPUT_MASKED=1): input reads are ALREADY host-N-masked
#     (e.g. output of the SMRTCap mask step / supplied pre-masked). Just strip.
#   * Path B (default): input is raw clean HiFi reads (host+provirus, no mask,
#     as they come from SRA). We first PRODUCE the N-masked form ourselves:
#     map each read to HXB2 with minimap2 (HiFi preset), and N-mask the
#     soft-clipped read ends (the host flanks that don't align to HIV). This
#     yields exactly the NNN...ACGT...NNN format above, which the strip tool
#     then reduces to the proviral core. Reads that don't map to HIV at all
#     are dropped (no provirus). minimap2 + samtools + awk + seqkit only.
# Usage: ./extract_provirus_pacbio.sh   (via sbatch scripts/utils/run_comparison_step.slurm.sh)
set -uo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
source "${REPO_ROOT}/scripts/common/lib_compare.sh"

REF_FASTA="${REPO_ROOT}/data/reference/K03455.1.fasta"
QC_DIR="${REPO_ROOT}/results/download_qc/pacbio"
RAW_DIR="${REPO_ROOT}/data/raw/pacbio"
MASK_DIR="${REPO_ROOT}/data/processed/pacbio/masked"
RESULTS_DIR="${REPO_ROOT}/results/proviral_extraction/pacbio"
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
STRIP_TOOL="${REPO_ROOT}/scripts/utils/extract_provirus_strip_hostN.sh"
mkdir -p "${RESULTS_DIR}" "${MASK_DIR}"
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

if [ ! -s "${REF_FASTA}" ]; then
    echo "ERROR: HXB2 reference ${REF_FASTA} not found." >&2
    exit 1
fi
THREADS="${THREADS:-4}"
export THREADS
PROVIRUS_INPUT_MASKED="${PROVIRUS_INPUT_MASKED:-0}"

# Resolve the cleanest available reads for a sample: prefer deduped, then
# Kraken2-cleaned, then filtered, then raw.
resolve_reads() {
    local srr="$1"
    local lm="${RAW_DIR}/local_masked"
    # Highest priority: user-provided local host-N-masked HiFi (Path A). The
    # SMRTcap naming is "<sample>.fastq.hiv.unmasked.fa" -- "unmasked" = the HIV
    # provirus is the unmasked (ACGT) part, host is N-masked. Search RECURSIVELY
    # under local_masked/ so an extra folder level from scp (e.g. the copied
    # raw_smrtcap/ wrapper) doesn't hide the files. Prefer the exact SMRTcap
    # name, then any .fa/.fasta whose basename starts with the sample id.
    if [ -d "${lm}" ]; then
        local hit
        hit=$(find "${lm}" -type f -iname "${srr}.fastq.hiv.unmasked.fa" 2>/dev/null | head -1)
        [ -z "${hit}" ] && hit=$(find "${lm}" -type f \( -iname "${srr}*.fa" -o -iname "${srr}*.fasta" \) 2>/dev/null | head -1)
        [ -n "${hit}" ] && [ -s "${hit}" ] && { echo "${hit}"; return; }
    fi
    for cand in \
        "${QC_DIR}/dedup_fastp_out/${srr}.dedup.fastq.gz" \
        "${QC_DIR}/kraken2_fastp_out/${srr}.kraken_filtered.fastq.gz" \
        "${QC_DIR}/fastp_out/${srr}.filtered.fastq.gz" \
        "${RAW_DIR}/${srr}.fastq.gz"; do
        [ -s "${cand}" ] && { echo "${cand}"; return; }
    done
    echo ""
}

for SRR in $(subset_accessions pacbio "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do
    READS="$(resolve_reads "${SRR}")"
    if [ -z "${READS}" ]; then
        echo "WARNING: no reads found for ${SRR} (run download + download_qc first), skipping." >&2
        continue
    fi
    echo "=== proviral extraction for ${SRR} (input: ${READS##*/}) ==="

    MASKED="${MASK_DIR}/${SRR}.masked.fasta"
    TIMELOG="${RESULTS_DIR}/extract_${SRR}.time"
    LOG="${RESULTS_DIR}/extract_${SRR}.log"
    PROVIRUS="${RESULTS_DIR}/${SRR}.provirus.fasta"
    COORDS="${RESULTS_DIR}/${SRR}.provirus_coords.tsv"

    if [ "${PROVIRUS_INPUT_MASKED}" = "1" ]; then
        # Path A: input already N-masked; convert fastq->fasta if needed.
        case "${READS}" in
            *.fastq.gz|*.fq.gz) seqkit fq2fa "${READS}" -o "${MASKED}" 2>/dev/null ;;
            *) cp "${READS}" "${MASKED}" ;;
        esac
    else
        # Path B: build the N-masked form by mapping to HXB2 and masking the
        # soft-clipped (host) read ends. -F 0x904 keeps only primary mapped
        # alignments (drops unmapped 0x4, secondary 0x100, supplementary
        # 0x800), so each surviving read appears once. Leading/trailing S in
        # the CIGAR give the host-flank lengths to mask.
        minimap2 -a -x map-hifi --secondary=no -t "${THREADS}" "${REF_FASTA}" "${READS}" 2>"${LOG}" \
            | samtools view -F 0x904 - 2>>"${LOG}" \
            | awk '
                {
                    qname=$1; cigar=$6; seq=$10;
                    lead=0; trail=0;
                    if (match(cigar, /^[0-9]+S/)) lead = substr(cigar, RSTART, RLENGTH-1) + 0;
                    if (match(cigar, /[0-9]+S$/)) trail = substr(cigar, RSTART, RLENGTH-1) + 0;
                    L = length(seq);
                    mid = substr(seq, lead+1, L-lead-trail);
                    np = sprintf("%*s", lead, "");  gsub(/ /, "N", np);
                    ns = sprintf("%*s", trail, ""); gsub(/ /, "N", ns);
                    print ">" qname;
                    print np mid ns;
                }
            ' > "${MASKED}" 2>>"${LOG}"
    fi

    N_MASKED=$(grep -c "^>" "${MASKED}" 2>/dev/null || echo 0)
    if [ "${N_MASKED}" -eq 0 ]; then
        echo "WARNING: no HIV-mapping reads for ${SRR}; no provirus extracted." >&2
        append_summary_row "proviral_extraction_pacbio" "minimap2mask+strip" "${SRR}" "" "" "1" "0" "0 reads mapped to HXB2"
        continue
    fi

    # The requested step: strip host N-flanks -> proviral cores.
    measure_and_run "${TIMELOG}" -- \
        "${STRIP_TOOL}" "${MASKED}" "${PROVIRUS}" "${COORDS}" >> "${LOG}" 2>&1
    EXIT_CODE=$?
    parse_time_metrics "${TIMELOG}"

    N_PROV=$(grep -c "^>" "${PROVIRUS}" 2>/dev/null || echo 0)
    VALID=0; METRIC="n/a"
    if [ "${N_PROV}" -gt 0 ]; then
        VALID=1
        MEDLEN=$(awk -F'\t' 'NR>1 && $8=="kept"{print $5}' "${COORDS}" 2>/dev/null | sort -n | awk '{a[NR]=$1} END{if(NR)print (NR%2)?a[(NR+1)/2]:int((a[NR/2]+a[NR/2+1])/2)}')
        METRIC="${N_PROV}/${N_MASKED} reads yielded a proviral core; median core ${MEDLEN:-?}bp"
    fi
    append_summary_row "proviral_extraction_pacbio" "minimap2mask+strip" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"
done

echo "Done. Proviral cores in ${RESULTS_DIR}/<SRR>.provirus.fasta ; see ${SUMMARY_TSV}"
