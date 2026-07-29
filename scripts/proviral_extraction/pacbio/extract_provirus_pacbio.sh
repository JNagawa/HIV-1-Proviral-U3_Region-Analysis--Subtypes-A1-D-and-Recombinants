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
# -u errors on unset vars, pipefail catches a failing pipe stage; no -e so one bad sample doesn't
# abort the run
set -uo pipefail

# absolute path of this script's own dir, so paths work from any launch location
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
# repo top-level (three dirs up), the base for every other path below
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
# load shared helpers: measure_and_run, parse_time_metrics, append_summary_row, subset_accessions
source "${REPO_ROOT}/scripts/common/lib_compare.sh"

# HXB2 reference genome reads are mapped against to find the provirus
REF_FASTA="${REPO_ROOT}/data/reference/K03455.1.fasta"
# where the QC step wrote filtered/cleaned reads
QC_DIR="${REPO_ROOT}/results/download_qc/pacbio"
# downloaded raw HiFi FASTQs (and any local pre-masked inputs)
RAW_DIR="${REPO_ROOT}/data/raw/pacbio"
# intermediate host-N-masked FASTAs produced here
MASK_DIR="${REPO_ROOT}/data/processed/pacbio/masked"
# output: proviral cores, coords, timing, summary
RESULTS_DIR="${REPO_ROOT}/results/proviral_extraction/pacbio"
# the TSV each sample appends a timing/validity row to
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
# the actual N-flank stripper that yields the proviral core
STRIP_TOOL="${REPO_ROOT}/scripts/utils/extract_provirus_strip_hostN.sh"
mkdir -p "${RESULTS_DIR}" "${MASK_DIR}"              # ensure output and intermediate dirs exist
# seed the manual notes file from the template on first run
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" "${RESULTS_DIR}/ease_of_use_notes.md"

# can't map to HIV without the reference, so fail early...
if [ ! -s "${REF_FASTA}" ]; then
    echo "ERROR: HXB2 reference ${REF_FASTA} not found." >&2  # ...with a clear error on stderr...
    exit 1                                           # ...and a non-zero exit
fi
# honour an externally-set THREADS, otherwise default to 4
THREADS="${THREADS:-4}"
export THREADS                                       # export so child tools/scripts inherit it
# 0 = raw reads (Path B, mask then strip); 1 = input already N-masked (Path A, just strip)
PROVIRUS_INPUT_MASKED="${PROVIRUS_INPUT_MASKED:-0}"

# Resolve the cleanest available reads for a sample: prefer deduped, then
# Kraken2-cleaned, then filtered, then raw.
# pick the best available reads for a sample, echoing the chosen path
resolve_reads() {
    local srr="$1"                                   # arg 1 = accession/sample id
    # dir holding user-supplied pre-masked SMRTcap FASTAs (Path A)
    local lm="${RAW_DIR}/local_masked"
    # Highest priority: user-provided local host-N-masked HiFi (Path A). The
    # SMRTcap naming is "<sample>.fastq.hiv.unmasked.fa" -- "unmasked" = the HIV
    # provirus is the unmasked (ACGT) part, host is N-masked. Search RECURSIVELY
    # under local_masked/ so an extra folder level from scp (e.g. the copied
    # raw_smrtcap/ wrapper) doesn't hide the files. Prefer the exact SMRTcap
    # name, then any .fa/.fasta whose basename starts with the sample id.
    if [ -d "${lm}" ]; then                          # only search if a local_masked dir exists
        local hit                                    # will hold the matched pre-masked file, if any
        # prefer the exact SMRTcap name (recursive, case-insensitive)
        hit=$(find "${lm}" -type f -iname "${srr}.fastq.hiv.unmasked.fa" 2>/dev/null | head -1)
        # else any .fa/.fasta whose basename starts with the sample id
        [ -z "${hit}" ] && hit=$(find "${lm}" -type f \( -iname "${srr}*.fa" -o -iname "${srr}*.fasta" \) 2>/dev/null | head -1)
        # if a non-empty match was found, use it (Path A) and stop
        [ -n "${hit}" ] && [ -s "${hit}" ] && { echo "${hit}"; return; }
    fi
    # otherwise try QC outputs from cleanest to rawest, then raw reads
    for cand in \
        "${QC_DIR}/dedup_fastp_out/${srr}.dedup.fastq.gz" \
        "${QC_DIR}/kraken2_fastp_out/${srr}.kraken_filtered.fastq.gz" \
        "${QC_DIR}/fastp_out/${srr}.filtered.fastq.gz" \
        "${RAW_DIR}/${srr}.fastq.gz"; do
        # return the first candidate that exists and is non-empty
        [ -s "${cand}" ] && { echo "${cand}"; return; }
    done
    # nothing usable found; echo empty so the caller can detect it
    echo ""
}

# loop over just the PacBio accessions chosen for this comparison
for SRR in $(subset_accessions pacbio "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do
    READS="$(resolve_reads "${SRR}")"                # locate the best input reads for this sample
    if [ -z "${READS}" ]; then                       # if none were found...
        # ...warn...
        echo "WARNING: no reads found for ${SRR} (run download + download_qc first), skipping." >&2
        continue                                     # ...and move on to the next sample
    fi
    # progress marker showing which input was picked
    echo "=== proviral extraction for ${SRR} (input: ${READS##*/}) ==="

    # host-N-masked FASTA (produced or converted below)
    MASKED="${MASK_DIR}/${SRR}.masked.fasta"
    # file where measure_and_run records wallclock/RSS for the strip step
    TIMELOG="${RESULTS_DIR}/extract_${SRR}.time"
    LOG="${RESULTS_DIR}/extract_${SRR}.log"          # captured stdout+stderr for this sample
    PROVIRUS="${RESULTS_DIR}/${SRR}.provirus.fasta"  # final proviral cores output
    # per-read provenance (lead/trail N, core length, status)
    COORDS="${RESULTS_DIR}/${SRR}.provirus_coords.tsv"

    if [ "${PROVIRUS_INPUT_MASKED}" = "1" ]; then    # Path A: input is already host-N-masked
        # Path A: input already N-masked; convert fastq->fasta if needed.
        case "${READS}" in                           # normalise the input to a FASTA at ${MASKED}
            # fastq input: convert to fasta
            *.fastq.gz|*.fq.gz) seqkit fq2fa "${READS}" -o "${MASKED}" 2>/dev/null ;;
            *) cp "${READS}" "${MASKED}" ;;          # already fasta: just copy it into place
        esac
    # Path B: raw reads, build the masked form by mapping to HXB2
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
                    qname=$1; cigar=$6; seq=$10;  # SAM read name, CIGAR string, and read sequence
                    lead=0; trail=0;              # host-flank lengths to N-mask at each end (default none)
                    if (match(cigar, /^[0-9]+S/)) lead = substr(cigar, RSTART, RLENGTH-1) + 0;   # leading soft-clip = 5-prime host flank length
                    if (match(cigar, /[0-9]+S$/)) trail = substr(cigar, RSTART, RLENGTH-1) + 0;  # trailing soft-clip = 3-prime host flank length
                    L = length(seq);              # total read length
                    mid = substr(seq, lead+1, L-lead-trail);  # the HIV-aligned core between the two clips
                    np = sprintf("%*s", lead, "");  gsub(/ /, "N", np);   # build a run of N of length lead (5-prime mask)
                    ns = sprintf("%*s", trail, ""); gsub(/ /, "N", ns);   # build a run of N of length trail (3-prime mask)
                    print ">" qname;             # emit FASTA header
                    print np mid ns;             # emit N-masked-flanks + ACGT core, the strip tool input format
                }
            ' > "${MASKED}" 2>>"${LOG}"          # minimap2 (HiFi preset) -> keep only primary mapped reads (-F 0x904) -> N-mask soft-clipped host ends
    fi

    # number of reads that mapped to HIV and made it into the masked FASTA
    N_MASKED=$(grep -c "^>" "${MASKED}" 2>/dev/null || echo 0)
    if [ "${N_MASKED}" -eq 0 ]; then                 # no HIV-mapping reads means nothing to extract
        # warn on stderr
        echo "WARNING: no HIV-mapping reads for ${SRR}; no provirus extracted." >&2
        # record a zero-result row and move on
        append_summary_row "proviral_extraction_pacbio" "minimap2mask+strip" "${SRR}" "" "" "1" "0" "0 reads mapped to HXB2"
        continue                                     # skip the strip step for this sample
    fi

    # The requested step: strip host N-flanks -> proviral cores.
    # time+run the N-flank stripper to produce the proviral cores and coords
    measure_and_run "${TIMELOG}" -- \
        "${STRIP_TOOL}" "${MASKED}" "${PROVIRUS}" "${COORDS}" >> "${LOG}" 2>&1
    # capture the strip tool's exit status before $? is overwritten
    EXIT_CODE=$?
    # set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file
    parse_time_metrics "${TIMELOG}"

    # number of reads that yielded a non-empty proviral core
    N_PROV=$(grep -c "^>" "${PROVIRUS}" 2>/dev/null || echo 0)
    VALID=0; METRIC="n/a"                            # assume invalid until proven otherwise
    if [ "${N_PROV}" -gt 0 ]; then                   # at least one core extracted = success
        VALID=1                                      # mark valid
        # median core length over kept reads
        MEDLEN=$(awk -F'\t' 'NR>1 && $8=="kept"{print $5}' "${COORDS}" 2>/dev/null | sort -n | awk '{a[NR]=$1} END{if(NR)print (NR%2)?a[(NR+1)/2]:int((a[NR/2]+a[NR/2+1])/2)}')
        # human-readable key metric
        METRIC="${N_PROV}/${N_MASKED} reads yielded a proviral core; median core ${MEDLEN:-?}bp"
    fi
    # write this sample's row to summary.tsv
    append_summary_row "proviral_extraction_pacbio" "minimap2mask+strip" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"
done

# final confirmation pointing at the outputs
echo "Done. Proviral cores in ${RESULTS_DIR}/<SRR>.provirus.fasta ; see ${SUMMARY_TSV}"
