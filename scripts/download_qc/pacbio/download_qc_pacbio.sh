#!/bin/bash
#SBATCH --job-name=pb_download_qc
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
# 64G because Kraken2 loads the whole 16GB standard DB into RAM and is run
# once after each filter; 24h is sized for the SRA FASTQ arm (five raw HiFi
# runs, ~12GB gzipped) where Kraken2 on 10-15kb reads is the long pole --
# merely downloading that data took 4h40m (job 111580). The default local
# FASTA arm is four tiny samples (19-75 reads) and finishes in minutes; the
# allocation is left at the worst case so one header covers both arms.

##################################################################################
# QC analysis of single-end PacBio HIV-SMRTcap data. Accepts BOTH FASTA and
# FASTQ input, decided per file from the DATA rather than the filename: the
# SMRTcap capture files are named "<sample>.fastq.hiv.unmasked.fa" but hold
# FASTA with no quality scores, so the extension cannot be trusted.
#
# STEP 0 runs before any QC: strip the host-N flanks.
#   The SMRTcap assay captures the integrated provirus together with its
#   flanking host genome, and the host portion arrives masked to N while the
#   provirus stays in ACGT. The virus integrates at a different host position
#   in every read, so each read looks like
#       NNNN...NNN  ACGT...proviral...ACGT  NNN...NNNN
#   and the sequence of interest is the ACGT core. Measured across the four
#   local samples the flanks are 43-78% of each read and there are ZERO
#   internal Ns, so stripping the two anchored N-runs removes the entire mask
#   and leaves pure-ACGT provirus. Every QC step below therefore runs on the
#   stripped core -- filtering before stripping would measure the host mask,
#   not the provirus.
#
# STEP 0b then puts every read on the same strand.
#   PacBio CCS names carry the strand as the FINAL path element: .../ccs/0 is
#   forward, .../ccs/1 is reverse. Confirmed by mapping all four local samples
#   to HXB2 (K03455.1) with minimap2 -- 77/77 reads ending in /0 map to '+' and
#   87/87 ending in /1 map to '-', and after this step every mapped read is on
#   '+'. Some names ALSO carry a fwd/rev word; that word contradicts the
#   alignment (.../ccs/fwd/1 maps to '-') and is deliberately ignored.
#   Orienting here means dedup, MSA and motif mapping never have to reason
#   about strand.
#
# The four QC steps, and which tools apply to which format:
#   1. NanoPlot, run TWICE (pre-strip and post-strip) so the read-length
#      distributions show exactly what the mask removal did.
#      FASTQ -> --fastq (length and quality); FASTA -> --fasta (length only).
#   2. length/quality FILTER, comparing:
#      FASTQ -> NanoFilt vs fastp vs chopper, all held to the same LEN_MIN/Q_MIN
#               so the comparison reflects the tools, not mismatched settings.
#      FASTA -> seqkit vs awk, length only. NanoFilt and chopper ABORT on FASTA,
#               and fastp exits 0 while writing an EMPTY file, so all three get
#               explicit not_applicable rows rather than being silently skipped
#               or, worse, recorded as a legitimate zero-read result.
#   3. Kraken2 host/bacterial removal (single-end) after each filter. On the
#      already-stripped SMRTcap reads this is a CHECK, not a cleanup: the mask
#      should have removed the host already, so a large Kraken2 loss here means
#      the mask under-called host sequence.
#   4. DEDUP. seqkit and awk are each run under BOTH matching rules, because
#      they are different questions and the tools do not default to the same
#      one: *_exact is byte-identical only (seqkit -s -P), *_bothstrands also
#      collapses a sequence against its reverse complement (seqkit -s, the
#      default). After Step 0b the two should agree, so their agreement doubles
#      as a check that orientation worked. fastp --dedup joins in on FASTQ only.
#
# LENGTH BAR: LEN_MIN is left unset by default and chosen per format (500 for
# the SMRTcap FASTA arm, 1000 for raw HiFi FASTQ) -- see the LEN_MIN block below
# for the measurements behind those numbers.
#
# INPUT RESOLUTION is local-first: each sample is looked up in LOCAL_INPUT_DIR
# (nested <dir>/<sample>/ layout first, then flat <dir>/), and only when nothing
# is found there, and the name looks like an SRA accession, is it downloaded.
#
# Usage: sbatch scripts/download_qc/pacbio/download_qc_pacbio.sh
#   Submit from the repo root -- the --output/--error paths above are relative to
#   the submitting directory, so logs/ must exist where you run sbatch.
#   Still works unchanged as `./download_qc_pacbio.sh` for a quick direct run, and
#   via `sbatch scripts/utils/run_comparison_step.slurm.sh <this script>`; the
#   #SBATCH lines are ordinary comments in both of those cases.
#
# Environment overrides:
#   SAMPLE_SET       subset_samples.tsv rows to loop over (default: pacbio)
#   LOCAL_INPUT_DIR  where to look for input before downloading
#   FORCE_FORMAT     fasta|fastq, overrides per-file format detection
#   SKIP_STRIP       1 to process input as-is, without removing host-N flanks
#   LEN_MIN, Q_MIN, THREADS
#   e.g. the SRA FASTQ arm:
#     SAMPLE_SET=pacbio_sra sbatch scripts/download_qc/pacbio/download_qc_pacbio.sh
#
# -u errors on unset vars, pipefail fails a pipe if any stage fails; no -e so one tool failing
# doesn't kill the whole comparison
#######################################################################################
set -uo pipefail

# Activate the tool env ourselves so this is submittable on its own, not only
# through run_comparison_step.slurm.sh (which does the same thing). Skipped when
# the env is already active, so a direct run in an activated shell is untouched.
if [ "${CONDA_DEFAULT_ENV:-}" != "HIV_U3analysis" ]; then
    CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"  # usual miniconda hook location
    # source it, or fall back to whatever conda is on PATH
    if [ -f "$CONDA_SH" ]; then source "$CONDA_SH"; else source "$(conda info --base)/etc/profile.d/conda.sh"; fi
    conda activate HIV_U3analysis                    # env holding NanoPlot/fastp/chopper/Kraken2/seqkit
fi

# Locate the repo root. A plain `$(dirname "$0")/../../..` walk is NOT enough:
# under `sbatch <this script>` Slurm copies the script into a spool directory
# (/var/spool/slurmd/job<N>/slurm_script), so $0 no longer sits inside the repo
# and the walk lands on /var. Try each candidate in turn and accept the first
# one that actually contains the shared helper library:
#   1. SLURM_SUBMIT_DIR -- the directory sbatch was invoked from (set only under Slurm)
#   2. the $0-relative walk -- correct for a direct run or via run_comparison_step.slurm.sh
#   3. the current working directory -- last resort
REPO_ROOT=""                                         # filled in by the loop below
for CAND in "${SLURM_SUBMIT_DIR:-}" \
            "$(cd "$(dirname "$0")/../../.." 2>/dev/null && pwd)" \
            "$(pwd)"; do
    # a repo root is only a repo root if the helper library is under it
    if [ -n "${CAND}" ] && [ -f "${CAND}/scripts/common/lib_compare.sh" ]; then
        REPO_ROOT="${CAND}"                          # first match wins
        break
    fi
done
if [ -z "${REPO_ROOT}" ]; then                       # none of the candidates panned out...
    # ...so fail immediately rather than writing results into the wrong tree
    echo "ERROR: cannot locate the repo root (no scripts/common/lib_compare.sh found)." >&2
    exit 1
fi
# load shared helpers: measure_and_run, parse_time_metrics, append_summary_row, subset_accessions
source "${REPO_ROOT}/scripts/common/lib_compare.sh"

##################################################################################
#  Configuration
##################################################################################
# which subset_samples.tsv rows to process; `pacbio` is the four local SMRTcap
# samples, `pacbio_sra` the five raw SRA HiFi runs
SAMPLE_SET="${SAMPLE_SET:-pacbio}"
# searched first for every sample, before any download is attempted
LOCAL_INPUT_DIR="${LOCAL_INPUT_DIR:-${REPO_ROOT}/data/raw/pacbio/local_masked/raw_smrtcap}"
# where SRA downloads land, and where an earlier download would already be
DOWNLOAD_DIR="${REPO_ROOT}/data/raw/pacbio"
# scratch dir for the large prefetched .sra files
PREFETCH_DIR="${DOWNLOAD_DIR}/.sra"
# output: all QC results + summary go here
RESULTS_DIR="${REPO_ROOT}/results/download_qc/pacbio"
# the single TSV every tool appends a timing/validity row to
SUMMARY_TSV="${RESULTS_DIR}/summary.tsv"
# Kraken2 database used for host/bacterial read removal
KRAKEN2_DB="${REPO_ROOT}/data/reference/kraken2_standard_16gb_db"
# where the N-stripped proviral cores and their coordinate logs are written
PROV_DIR="${RESULTS_DIR}/provirus_stripped"

# create the results dirs (and parents) if they don't exist
mkdir -p "${RESULTS_DIR}" "${PROV_DIR}"
# seed the manual notes file from the template only on the first run
[ -f "${RESULTS_DIR}/ease_of_use_notes.md" ] \
    || cp "${REPO_ROOT}/scripts/common/ease_of_use_template.md" \
          "${RESULTS_DIR}/ease_of_use_notes.md"

# Minimum read length. Left EMPTY here on purpose: the sensible bar differs by
# arm, so the effective value (LEN_CUT) is chosen per sample once the format is
# known. Setting LEN_MIN in the environment overrides both defaults.
#
#   FASTA / SMRTcap arm -> 500.  Measured over the four local samples, the
#     shortest proviral core after N-stripping is 763bp and retention is FLAT at
#     165/165 reads for every cutoff from 0 to 763 -- there is no short junk to
#     remove. The first thing a higher bar deletes is real biology: at 800+,
#     sample 124_4 drops 41 reads to 16, losing all 24 reads of its ~765bp
#     3'LTR population (HXB2 8985-9719), which is exactly the U3-bearing
#     fragment this project exists to analyse. 500 keeps everything while still
#     being a meaningful floor, since a read below ~450bp cannot carry a
#     complete U3 and is useless here anyway.
#   FASTQ / SRA arm -> 1000. Those are raw unmasked 10-15kb HiFi reads with no
#     host mask and no deletion structure, so the old bar still applies.
LEN_MIN="${LEN_MIN:-}"                               # empty = pick the per-format default below
# minimum mean read quality (Phred) kept by every FASTQ filter
Q_MIN="${Q_MIN:-20}"
# set to 1 to process input as-is, skipping the host-N flank removal
SKIP_STRIP="${SKIP_STRIP:-0}"
# thread count; honour an externally-set THREADS, else Slurm's allocation, else 4
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
# export so child scripts/tools (e.g. the Kraken2 wrapper) inherit it
export THREADS

##################################################################################
#  Pipeline steps (inlined)
#
#  Every processing step from the host-N strip onward lives in this file, so the
#  whole QC pipeline is one self-contained script with no sibling helper scripts
#  to keep in sync. Each is exported with `export -f` because measure_and_run
#  runs its command in a child process (`bash -c 'fn "$@"' _ ...`), and only
#  exported functions survive that boundary.
#
#  All five accept FASTA or FASTQ, gzipped or plain, and detect which from the
#  DATA rather than the filename -- the SMRTcap inputs are called
#  "<sample>.fastq.hiv.unmasked.fa" but hold FASTA, so extensions cannot be
#  trusted anywhere in this pipeline.
##################################################################################

# --- Step 0: strip the host-N flanks, keeping the ACGT proviral core. ---------
# The SMRTcap convention is that host genome is N-masked and the integrated
# provirus stays in ACGT. The virus integrates at a different host position in
# every read, so the flanks vary in length and sit at BOTH ends:
#     NNNNNACGGTTCTAGGGTTTCCACTANNNNN  ->  ACGGTTCTAGGGTTTCCACTA
# Only the anchored end-runs are removed; any internal base is preserved, so a
# genuine internal N never truncates the genome. For FASTQ the quality string is
# sliced with the same coordinates, keeping bases and qualities in register.
#
# Usage: strip_host_n <in> <out[.gz]> [coords.tsv]
#   coords.tsv gets one row per INPUT record:
#     name orig_len lead_N trail_N provirus_len prov_start prov_end status
#   (1-based inclusive coords; status = kept | dropped_all_N | dropped_empty)
strip_host_n() {
    local in="$1" out="$2" coords="${3:-}"           # input, output, optional provenance log
    [ -s "${in}" ] || { echo "ERROR: strip_host_n: '${in}' missing or empty." >&2; return 1; }
    mkdir -p "$(dirname "${out}")"                   # ensure the output dir exists
    # truncate any coords file from a previous run: the awk below APPENDS rows
    [ -n "${coords}" ] && { mkdir -p "$(dirname "${coords}")"; : > "${coords}"; }

    local fmt fx2tab_opts out_cmd                    # format, fx2tab flags, output filter
    fmt=$(detect_format "${in}")                     # '>' = fasta, '@' = fastq
    # anything that is not sequence data cannot be stripped
    [ "${fmt}" = "unknown" ] && { echo "ERROR: strip_host_n: '${in}' is not FASTA or FASTQ." >&2; return 1; }
    # -q makes fx2tab carry the quality string as a third column
    [ "${fmt}" = "fastq" ] && fx2tab_opts="-q" || fx2tab_opts=""
    # gzip the output when the caller asked for a .gz path
    [ "${out%.gz}" != "${out}" ] && out_cmd="gzip -c" || out_cmd="cat"

    # shellcheck disable=SC2086
    zcat -f "${in}" 2>/dev/null | seqkit fx2tab ${fx2tab_opts} 2>/dev/null \
      | awk -F'\t' -v coords="${coords}" -v fmt="${fmt}" '
        {
            name = $1                    # record name (first tab field)
            seq  = $2                    # linearised sequence (second tab field)
            qual = (fmt == "fastq" ? $3 : "")  # quality string, FASTQ only
            orig = length(seq)           # original read length before stripping

            lead = 0                     # count of leading N bases (the 5-prime host flank)
            while (lead < orig && substr(seq, lead+1, 1) ~ /[Nn]/) lead++  # advance past every leading N

            if (lead == orig) {            # sequence is entirely N -> no provirus
                if (coords != "")
                    printf "%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n", name, orig, orig, 0, 0, "NA", "NA", "dropped_all_N" >> coords
                next                     # emit no core, move to the next record
            }

            trail = 0                    # count of trailing N bases (the 3-prime host flank)
            while (trail < orig && substr(seq, orig-trail, 1) ~ /[Nn]/) trail++  # advance inward past every trailing N

            core_len   = orig - lead - trail  # length of the ACGT core left after removing both flanks
            prov_start = lead + 1             # 1-based start of the core in original coordinates
            prov_end   = orig - trail         # 1-based inclusive end of the core
            core       = substr(seq, prov_start, core_len)  # the proviral core itself

            if (core_len <= 0) {           # defensive: nothing left after stripping
                if (coords != "")
                    printf "%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n", name, orig, lead, trail, 0, "NA", "NA", "dropped_empty" >> coords
                next                     # emit no core
            }

            # the quality string is cut with the SAME start/length as the sequence
            if (fmt == "fastq")
                printf "%s\t%s\t%s\n", name, core, substr(qual, prov_start, core_len)
            else
                printf "%s\t%s\n", name, core
            if (coords != "")
                printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", name, orig, lead, trail, core_len, prov_start, prov_end, "kept" >> coords
        }
      ' | seqkit tab2fx -w0 2>/dev/null | ${out_cmd} > "${out}"

    # prepend the coords header (awk appended rows without one)
    if [ -n "${coords}" ] && [ -f "${coords}" ]; then
        local hdr tmp_c                              # header text and scratch file
        hdr="name\torig_len\tlead_N\ttrail_N\tprovirus_len\tprov_start\tprov_end\tstatus"
        tmp_c="$(mktemp)"                            # scratch file to prepend the header
        { printf "%b\n" "${hdr}"; cat "${coords}"; } > "${tmp_c}" && mv "${tmp_c}" "${coords}"
    fi
    # summary to stderr so it lands in the caller's log
    echo "strip_host_n: $(count_records "${out}")/$(count_records "${in}") records kept a proviral core (${fmt})" >&2
}

# --- Step 0b: put every read on the forward strand. ---------------------------
# PacBio CCS names carry the strand as the FINAL path element:
#     .../ccs/0   forward (+) -> keep as is
#     .../ccs/1   reverse (-) -> reverse complement it
# Verified by mapping all four local SMRTcap samples to HXB2 (K03455.1) with
# minimap2: 77/77 reads ending in /0 map to '+' and 87/87 ending in /1 map to
# '-', with zero exceptions, and every read maps to '+' afterwards. Some names
# ALSO carry a fwd/rev word; that word CONTRADICTS the alignment (.../ccs/fwd/1
# maps to '-'), so it is deliberately ignored and only the digit is used.
# Flipped records get ":oriented_fwd" appended so the change is never invisible.
# For FASTQ the quality string is REVERSED (not complemented) alongside the bases.
#
# Usage: orient_forward <in> <out[.gz]> [report.tsv]
#   report.tsv gets one row per record: name  strand_tag  action
orient_forward() {
    local in="$1" out="$2" report="${3:-}"           # input, output, optional per-read log
    [ -s "${in}" ] || { echo "ERROR: orient_forward: '${in}' missing or empty." >&2; return 1; }
    mkdir -p "$(dirname "${out}")"                   # ensure the output dir exists
    # truncate any report from a previous run: the awk below APPENDS rows
    [ -n "${report}" ] && { mkdir -p "$(dirname "${report}")"; : > "${report}"; }

    local fmt fx2tab_opts out_cmd                    # format, fx2tab flags, output filter
    fmt=$(detect_format "${in}")                     # '>' = fasta, '@' = fastq
    [ "${fmt}" = "unknown" ] && { echo "ERROR: orient_forward: '${in}' is not FASTA or FASTQ." >&2; return 1; }
    [ "${fmt}" = "fastq" ] && fx2tab_opts="-q" || fx2tab_opts=""   # carry quality through for FASTQ
    [ "${out%.gz}" != "${out}" ] && out_cmd="gzip -c" || out_cmd="cat"  # gzip if asked

    # shellcheck disable=SC2086
    zcat -f "${in}" 2>/dev/null | seqkit fx2tab ${fx2tab_opts} 2>/dev/null \
      | awk -F'\t' -v report="${report}" -v fmt="${fmt}" '
        function revcomp(s,   i, n, out, c) {
            n = length(s); out = ""
            for (i = n; i >= 1; i--) {       # walk the sequence backwards...
                c = substr(s, i, 1)
                out = out ((c in comp) ? comp[c] : c)  # ...complementing each base
            }
            return out
        }
        function reverse(s,   i, n, out) {    # plain reversal, for the quality string
            n = length(s); out = ""
            for (i = n; i >= 1; i--) out = out substr(s, i, 1)
            return out
        }
        BEGIN {
            comp["A"]="T"; comp["C"]="G"; comp["G"]="C"; comp["T"]="A"; comp["U"]="A"; comp["N"]="N"
            comp["a"]="t"; comp["c"]="g"; comp["g"]="c"; comp["t"]="a"; comp["u"]="a"; comp["n"]="n"
        }
        {
            name = $1                        # record name (first tab field)
            seq  = $2                        # linearised sequence
            qual = (fmt == "fastq" ? $3 : "")  # quality string, FASTQ only

            id = name                        # the id ends at the first space
            sub(/[ \t].*$/, "", id)          # drop any trailing description
            tag = id
            sub(/^.*\//, "", tag)            # keep only what follows the final slash

            if (tag == "1") {                  # reverse-strand read -> flip it
                seq = revcomp(seq)             # reverse complement the bases
                if (fmt == "fastq") qual = reverse(qual)  # qualities reverse but do NOT complement
                name = name ":oriented_fwd"    # mark the record as transformed
                action = "revcomp"; nrev++
            } else if (tag == "0") {           # forward-strand read -> already correct
                action = "kept"; nfwd++
            } else {                           # name does not end in /0 or /1
                action = "kept_untagged"; nunk++   # leave it alone rather than guess
            }

            if (report != "")
                printf "%s\t%s\t%s\n", id, tag, action >> report

            if (fmt == "fastq")
                printf "%s\t%s\t%s\n", name, seq, qual
            else
                printf "%s\t%s\n", name, seq
        }
        END {
            printf "orient_forward: %d forward kept, %d reverse-complemented, %d untagged (%s)\n", nfwd, nrev, nunk, fmt > "/dev/stderr"
        }
      ' | seqkit tab2fx -w0 2>/dev/null | ${out_cmd} > "${out}"

    # prepend the report header (awk appended rows without one)
    if [ -n "${report}" ] && [ -f "${report}" ]; then
        local tmp_r                                  # scratch file to prepend the header
        tmp_r="$(mktemp)"
        { printf "%b\n" "name\tstrand_tag\taction"; cat "${report}"; } > "${tmp_r}" && mv "${tmp_r}" "${report}"
    fi
}

# --- Filter: pure-awk minimum-length filter. ----------------------------------
# The independent second implementation the FASTA arm needs to hold seqkit
# against: on FASTA, NanoFilt and chopper abort outright and fastp exits 0 while
# writing an EMPTY file, which would otherwise leave seqkit unopposed. Verified
# to agree with `seqkit seq -m` on identical inputs.
#
# Usage: filter_len_awk <in> <out.gz> <min_length>
filter_len_awk() {
    local in="$1" out="$2" minlen="$3"               # input, gzipped output, shortest record to keep
    mkdir -p "$(dirname "${out}")"                   # ensure the output dir exists
    if [ "$(detect_format "${in}")" = "fasta" ]; then
        # FASTA: records may wrap over many lines, so accumulate until the next '>'
        zcat -f "${in}" 2>/dev/null | awk -v m="${minlen}" '
            /^>/ {
                if (n && length(s) >= m) print h ORS s   # flush the previous record if long enough
                h = $0; s = ""; n = 1                    # start collecting the new record
                next
            }
            { s = s $0 }                                 # append this sequence line
            END { if (n && length(s) >= m) print h ORS s }  # flush the final record
        ' | gzip > "${out}"
    else
        # FASTQ: fixed 4-line records (header, sequence, separator, quality)
        zcat -f "${in}" 2>/dev/null | awk -v m="${minlen}" '
            NR % 4 == 1 { h = $0 }                       # line 1 = @header
            NR % 4 == 2 { s = $0 }                       # line 2 = sequence
            NR % 4 == 3 { p = $0 }                       # line 3 = + separator
            NR % 4 == 0 { if (length(s) >= m) print h ORS s ORS p ORS $0 }  # line 4 = quality
        ' | gzip > "${out}"
    fi
}

# --- Dedup: pure-awk duplicate removal, two matching rules. -------------------
# `seqkit rmdup -s` and a naive exact string match are NOT the same rule, and on
# real SMRTcap data they disagree (sample 203_3: 12 vs 15 records from the same
# 19 inputs). Both rules are provided so the seqkit-vs-awk comparison is
# like-for-like rather than comparing two different questions:
#
#   exact        byte-identical sequences only.        Equivalent to seqkit rmdup -s -P
#   bothstrands  a sequence and its REVERSE COMPLEMENT. Equivalent to seqkit rmdup -s
#                (seqkit's default, which compares both strands unless -P is given)
#
# After orient_forward the two should agree, since a read and its former reverse
# complement become byte-identical -- so their agreement is a useful check that
# the orientation step did what it claims.
#
# Sequences are held in an awk array, so peak memory scales with the number of
# UNIQUE bases. Trivial for the SMRTcap samples (19-75 reads), but the reason
# this is offered alongside seqkit rather than instead of it.
#
# Usage: rmdup_awk <in> <out.gz> [exact|bothstrands]
rmdup_awk() {
    local in="$1" out="$2" mode="${3:-exact}"        # input, gzipped output, matching rule
    case "${mode}" in
        exact|bothstrands) ;;                        # the two supported rules
        *)
            # a typo here would silently corrupt the comparison, so fail loudly
            echo "ERROR: rmdup_awk: mode must be 'exact' or 'bothstrands', got '${mode}'." >&2
            return 1 ;;
    esac
    mkdir -p "$(dirname "${out}")"                   # ensure the output dir exists

    # the key builder shared by both format branches: in exact mode the key IS
    # the sequence; in bothstrands mode it is whichever of the sequence and its
    # reverse complement sorts first, so the two collapse onto one key
    local keyfunc='
        function revcomp(s,   i, n, out, c) {
            n = length(s); out = ""
            for (i = n; i >= 1; i--) { c = substr(s, i, 1); out = out ((c in comp) ? comp[c] : c) }
            return out
        }
        function key(s,   r) {
            if (mode != "bothstrands") return s
            r = revcomp(s)
            return (s < r) ? s : r
        }
        BEGIN {
            comp["A"]="T"; comp["C"]="G"; comp["G"]="C"; comp["T"]="A"; comp["U"]="A"; comp["N"]="N"
            comp["a"]="t"; comp["c"]="g"; comp["g"]="c"; comp["t"]="a"; comp["u"]="a"; comp["n"]="n"
        }
    '
    if [ "$(detect_format "${in}")" = "fasta" ]; then
        # FASTA: accumulate wrapped sequence lines, then test the joined sequence
        zcat -f "${in}" 2>/dev/null | awk -v mode="${mode}" "${keyfunc}"'
            /^>/ {
                if (n) { k = key(s); if (!(k in seen)) { seen[k]; print h ORS s } }  # first sighting -> keep
                h = $0; s = ""; n = 1                    # start collecting the new record
                next
            }
            { s = s $0 }                                 # append this sequence line
            END { if (n) { k = key(s); if (!(k in seen)) { seen[k]; print h ORS s } } }  # flush the final record
        ' | gzip > "${out}"
    else
        # FASTQ: fixed 4-line records; the sequence line supplies the key
        zcat -f "${in}" 2>/dev/null | awk -v mode="${mode}" "${keyfunc}"'
            NR % 4 == 1 { h = $0 }                       # line 1 = @header
            NR % 4 == 2 { s = $0 }                       # line 2 = sequence, the key source
            NR % 4 == 3 { p = $0 }                       # line 3 = + separator
            NR % 4 == 0 { k = key(s); if (!(k in seen)) { seen[k]; print h ORS s ORS p ORS $0 } }  # keep first per key
        ' | gzip > "${out}"
    fi
}

# --- Kraken2 single-end host/bacterial removal. -------------------------------
# Classify every read against a human+bacterial+viral database, then drop any
# read whose assigned taxID descends from Homo sapiens (9606) or Bacteria (2),
# keeping unclassified and viral reads (the real HIV-1 signal). Single-end
# because HiFi/CCS reads are one record with no mate.
#
# On the already-stripped SMRTcap reads this is a CHECK rather than a cleanup:
# the N-mask should have removed the host already, so a large loss here would
# mean the mask under-called host sequence.
#
# Usage: kraken2_filter_se <reads> <outdir> <sample> <kraken2_db_dir>
# Produces: <outdir>/<sample>.kraken_filtered.<fasta|fastq>.gz and .kreport
kraken2_filter_se() {
    local reads="$1" outdir="$2" sample="$3" db="$4"  # input, output dir, sample name, DB dir
    local nodes="${db}/nodes.dmp"                     # taxonomy tree, used for the ancestor walk
    # bail if the DB isn't a real extracted Kraken2 database
    [ -s "${nodes}" ] || { echo "ERROR: ${nodes} not found -- is '${db}' an extracted Kraken2 database?" >&2; return 1; }
    mkdir -p "${outdir}"                              # ensure the output dir exists

    local gz_opt fmt                                  # gzip flag for kraken2, and output format
    # Kraken2 must be TOLD its input is gzipped; probe the magic bytes rather
    # than the filename, since the SMRTcap files carry misleading extensions
    gzip -t "${reads}" 2>/dev/null && gz_opt="--gzip-compressed" || gz_opt=""
    fmt=$(detect_format "${reads}")                   # decides only the output filename
    [ "${fmt}" = "unknown" ] && fmt="fastq"           # fall back to the common case

    # classify each read; per-read output plus a hierarchical report
    # shellcheck disable=SC2086
    kraken2 --db "${db}" ${gz_opt} --threads "${THREADS:-4}" \
        --output "${outdir}/${sample}.kraken" \
        --report "${outdir}/${sample}.kreport" \
        "${reads}" || return 1

    awk '
        NR==FNR {                                     # first file = nodes.dmp: build taxid->parent
            split($0, f, "\t\\|\t")                   # nodes.dmp fields are separated by "\t|\t"
            parent[f[1]] = f[2]                       # parent[taxid] = parent_taxid
            next
        }
        {
            split($0, f, "\t")                        # second file = kraken output, plain-tab separated
            read_id = f[2]                            # field 2 = read id
            t = f[3]                                  # field 3 = assigned taxid
            keep = 1                                   # keep unless we hit human/bacteria
            depth = 0                                  # loop guard against cycles/broken trees
            while (t != "1" && t != "" && depth < 50) {   # climb ancestors until root (1) or missing
                if (t == "9606" || t == "2") { keep = 0; break }  # human (9606) or Bacteria (2) -> drop
                if (!(t in parent)) break             # no parent (e.g. unclassified taxid 0) -> stop, keep
                t = parent[t]                         # step up to the parent taxid
                depth++
            }
            if (keep) print read_id                   # emit ids of reads to retain
        }
    ' "${nodes}" "${outdir}/${sample}.kraken" > "${outdir}/${sample}.keep_read_ids.txt"

    local n_total n_keep                              # classified vs retained read counts
    n_total=$(wc -l < "${outdir}/${sample}.kraken")   # total classified reads
    n_keep=$(wc -l < "${outdir}/${sample}.keep_read_ids.txt")  # reads kept after host removal
    # progress summary to stderr
    echo "kraken2_filter_se ${sample}: ${n_keep}/${n_total} reads retained (human/bacterial discarded)." >&2

    # subset the reads to just the kept ids, writing back in the input's own format
    seqkit grep -f "${outdir}/${sample}.keep_read_ids.txt" "${reads}" \
        -o "${outdir}/${sample}.kraken_filtered.${fmt}.gz"
}

##################################################################################
#  Helpers
##################################################################################

# Report the record format of a sequence file from its first non-blank
# character ('>' = FASTA, '@' = FASTQ). Never trust the extension here: the
# SMRTcap inputs are named "*.fastq.hiv.unmasked.fa" and hold FASTA.
detect_format() {
    local f="$1"                                     # arg 1 = file to inspect
    if [ -n "${FORCE_FORMAT:-}" ]; then              # an explicit override wins
        echo "${FORCE_FORMAT}"; return               # ...so return it unchanged
    fi
    local c                                          # first non-blank character
    # zcat -f reads plain and gzipped input alike
    c=$(zcat -f "${f}" 2>/dev/null | awk 'NF{print substr($0,1,1); exit}')
    case "${c}" in
        '>') echo "fasta" ;;                         # FASTA record header
        '@') echo "fastq" ;;                         # FASTQ record header
        *)   echo "unknown" ;;                       # not sequence data we can handle
    esac
}

# Count sequence records in a FASTA/FASTQ file, gzipped or plain. seqkit covers
# every combination, so one line serves both formats; a missing, empty or
# unparseable file reports 0 rather than erroring.
count_records() {
    local f="$1"                                     # arg 1 = file to count
    # a missing or empty file has no records
    [ -s "${f}" ] || { echo 0; return; }
    local n                                          # record count from seqkit
    # column 4 of `seqkit stats -T` is num_seqs
    n=$(seqkit stats -T "${f}" 2>/dev/null | awk -F'\t' 'NR==2{print $4+0}')
    echo "${n:-0}"                                   # default to 0 if the parse failed
}

# Find this sample's input, preferring anything already on disk over a
# download. Search order: LOCAL_INPUT_DIR/<sample>/ (the layout the SMRTcap
# data actually uses), then LOCAL_INPUT_DIR/ flat, then a previous download in
# DOWNLOAD_DIR. Prints the path and returns 0, or returns 1 if nothing exists.
resolve_input() {
    local s="$1" hit                                 # arg 1 = sample name; hit = first match found
    # every sequence extension worth looking for, plain or gzipped
    local -a pats=( '*.fa' '*.fasta' '*.fna' '*.fq' '*.fastq'
                    '*.fa.gz' '*.fasta.gz' '*.fna.gz' '*.fq.gz' '*.fastq.gz' )
    local -a findargs=()                             # assembled "-name X -o -name Y ..." group
    local p                                          # loop variable over the patterns
    for p in "${pats[@]}"; do                        # build the OR-group once, reused below
        findargs+=( -name "${p}" -o )                # append this pattern plus an OR
    done
    unset 'findargs[${#findargs[@]}-1]'              # drop the trailing -o so the group is valid

    # 1. nested layout: <LOCAL_INPUT_DIR>/<sample>/<sample>.<ext>
    if [ -d "${LOCAL_INPUT_DIR}/${s}" ]; then
        hit=$(find "${LOCAL_INPUT_DIR}/${s}" -maxdepth 1 \( "${findargs[@]}" \) 2>/dev/null | sort | head -1)
        [ -n "${hit}" ] && { echo "${hit}"; return 0; }
    fi
    # 2. flat layout: <LOCAL_INPUT_DIR>/<sample>.<ext>
    if [ -d "${LOCAL_INPUT_DIR}" ]; then
        hit=$(find "${LOCAL_INPUT_DIR}" -maxdepth 1 -name "${s}.*" \( "${findargs[@]}" \) 2>/dev/null | sort | head -1)
        [ -n "${hit}" ] && { echo "${hit}"; return 0; }
    fi
    # 3. an earlier SRA download
    [ -s "${DOWNLOAD_DIR}/${s}.fastq.gz" ] && { echo "${DOWNLOAD_DIR}/${s}.fastq.gz"; return 0; }
    return 1                                         # nothing on disk for this sample
}

# Fetch a sample from SRA into DOWNLOAD_DIR. Only called when resolve_input
# found nothing locally. Refuses anything that is not an SRA-style accession,
# so a mistyped local sample name fails loudly instead of hitting the network.
download_sample() {
    local s="$1"                                     # arg 1 = accession to fetch
    if ! printf '%s' "${s}" | grep -qE '^[SED]RR[0-9]+$'; then
        # a local-only sample name has no SRA accession to fall back to
        echo "ERROR: '${s}' not found under ${LOCAL_INPUT_DIR} and is not an SRA accession, skipping." >&2
        return 1
    fi
    mkdir -p "${DOWNLOAD_DIR}" "${PREFETCH_DIR}"     # both dirs must exist before prefetch runs
    echo "=== prefetch ${s} (no local copy found) ==="  # progress marker
    # prefetch's default cap is 20G; raise it so a large HiFi run is never silently skipped
    prefetch --max-size 100G --output-directory "${PREFETCH_DIR}" "${s}" || return 1
    echo "=== fasterq-dump ${s} (single-end HiFi) ==="  # progress marker
    # HiFi runs are single-end: --concatenate-reads keeps the run as one stream
    fasterq-dump --threads "${THREADS}" --outdir "${DOWNLOAD_DIR}" \
        --concatenate-reads "${PREFETCH_DIR}/${s}/${s}.sra" || return 1
    # compress the result in place (-f overwrites any stale .gz)
    [ -s "${DOWNLOAD_DIR}/${s}.fastq" ] && gzip -f "${DOWNLOAD_DIR}/${s}.fastq"
    # succeed only if a non-empty FASTQ actually landed
    [ -s "${DOWNLOAD_DIR}/${s}.fastq.gz" ]
}

# Record a step that cannot run on this input format, so summary.tsv states WHY
# a tool is absent instead of leaving a gap. This matters most for fastp on
# FASTA, which exits 0 and writes an empty file -- without an explicit row that
# would be indistinguishable from a genuine zero-read result.
append_na_row() {
    local tool="$1" sample="$2" reason="$3"          # tool, sample, and the human-readable reason
    # NA timings/exit code, output_valid=NA, reason carried in the metric column
    append_summary_row "download_qc_pacbio" "${tool}" "${sample}" "NA" "NA" "NA" "NA" "not_applicable: ${reason}"
}

# Run NanoPlot on one file, choosing --fastq or --fasta from the format. On
# FASTA the quality panels are simply absent from the report (there are no
# quality scores to plot); the read-length distribution, which is the whole
# point of the pre/post-strip pair, is produced for both.
run_nanoplot() {
    local f="$1" outdir="$2" prefix="$3" fmt="$4"    # input, output dir, filename prefix, format
    # already done: keep reruns idempotent
    [ -s "${outdir}/NanoPlot-report.html" ] && return 0
    mkdir -p "${outdir}"                             # ensure the per-sample output dir exists
    local flag="--fastq"                             # FASTQ is the default input flag
    [ "${fmt}" = "fasta" ] && flag="--fasta"         # FASTA needs the other one
    # --tsv_stats gives machine-readable stats alongside the HTML report
    NanoPlot "${flag}" "${f}" --outdir "${outdir}" --prefix "${prefix}" \
        --threads "${THREADS}" --tsv_stats > "${RESULTS_DIR}/nanoplot_${prefix}.log" 2>&1
}

# measure_and_run executes its command in a child `bash -c`, and only exported
# functions cross that boundary -- without this the pipeline steps would be
# "command not found" the moment they are timed. Exported here, after every
# definition above, because `export -f` on a not-yet-defined function fails.
# detect_format and count_records ride along because the steps call them.
export -f detect_format count_records
export -f strip_host_n orient_forward filter_len_awk rmdup_awk kraken2_filter_se

##################################################################################
#  Main loop
##################################################################################
echo "=== QC run: SAMPLE_SET=${SAMPLE_SET} LOCAL_INPUT_DIR=${LOCAL_INPUT_DIR} ==="

for SRR in $(subset_accessions "${SAMPLE_SET}" "${REPO_ROOT}/scripts/common/subset_samples.tsv"); do
    echo ""                                          # blank line between samples in the log
    echo "########## ${SRR} ##########"              # sample banner

    # --- Input resolution: local first, download only as a fallback. ---
    READS=$(resolve_input "${SRR}")                  # look on disk before touching the network
    if [ -z "${READS}" ]; then                       # nothing local...
        # ...so try to fetch it, and skip the sample if that is not possible
        download_sample "${SRR}" || continue
        READS=$(resolve_input "${SRR}")              # re-resolve to pick up what was just downloaded
    fi
    if [ -z "${READS}" ] || [ ! -s "${READS}" ]; then  # still nothing usable...
        # ...warn and move on rather than crashing the whole comparison
        echo "WARNING: no usable input for ${SRR}, skipping." >&2
        continue
    fi

    FORMAT=$(detect_format "${READS}")               # decide the branch from the data itself
    if [ "${FORMAT}" = "unknown" ]; then             # not FASTA and not FASTQ...
        # ...so none of the QC steps below can be applied
        echo "WARNING: ${READS} is neither FASTA nor FASTQ, skipping ${SRR}." >&2
        continue
    fi
    EXT="fa"                                         # extension for this sample's intermediates
    [ "${FORMAT}" = "fastq" ] && EXT="fq"            # FASTQ intermediates get .fq
    N_RAW=$(count_records "${READS}")                # input record count, for the strip metric
    echo "input  : ${READS}"                         # which file was resolved
    echo "format : ${FORMAT} (${N_RAW} records)"     # and what it turned out to be

    # --- Step 0: strip the host-N flanks, leaving the ACGT proviral core. ---
    # Everything downstream runs on this, so lengths and filters describe the
    # provirus rather than the host mask wrapped around it.
    if [ "${SKIP_STRIP}" = "1" ]; then               # explicitly disabled...
        STRIPPED="${READS}"                          # ...so carry the raw input straight through
        echo "=== host-N strip SKIPPED (SKIP_STRIP=1) for ${SRR} ==="  # make the choice visible in the log
        append_na_row "strip_hostN" "${SRR}" "disabled by SKIP_STRIP=1"
    else
        STRIPPED="${PROV_DIR}/${SRR}.provirus.${EXT}.gz"   # the proviral cores for this sample
        COORDS="${PROV_DIR}/${SRR}.strip_coords.tsv"       # per-read flank coordinates + status
        STIME="${RESULTS_DIR}/strip_hostN_${SRR}.time"     # timing file for the strip step
        SLOG="${RESULTS_DIR}/strip_hostN_${SRR}.log"       # log for the strip step
        echo "=== strip host-N flanks on ${SRR} ==="       # progress marker
        # strip_host_n gzips its own output when the path ends in .gz
        measure_and_run "${STIME}" -- \
            bash -c 'strip_host_n "$@"' _ \
            "${READS}" "${STRIPPED}" "${COORDS}" > "${SLOG}" 2>&1
        SEXIT=$?                                     # capture the stripper's exit status
        # set WALLCLOCK_SEC / PEAK_RSS_MB from the .time file
        parse_time_metrics "${STIME}"

        N_PROV=$(count_records "${STRIPPED}")        # reads that kept a non-empty core
        SVALID=0; SMETRIC="n/a"                      # assume invalid until proven otherwise
        if [ -s "${STRIPPED}" ] && [ "${N_PROV}" -gt 0 ]; then
            # percentage of input bases that were host mask, straight from the coords log
            PCT_N=$(awk -F'\t' 'NR>1 && $2>0 {tot+=$2; mask+=$3+$4} END{if(tot>0) printf "%.1f", 100*mask/tot; else print "0.0"}' "${COORDS}" 2>/dev/null)
            # mark valid and record what the strip actually removed
            SVALID=1; SMETRIC="${N_PROV}/${N_RAW} reads kept a proviral core (${PCT_N:-0}% of bases were host-N flank)"
        fi
        # write the strip step's row to summary.tsv
        append_summary_row "download_qc_pacbio" "strip_hostN" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${SEXIT}" "${SVALID}" "${SMETRIC}"
        # without a proviral core there is nothing left to QC
        if [ "${SVALID}" != "1" ]; then
            echo "WARNING: host-N strip produced no proviral sequence for ${SRR}, skipping its QC." >&2
            continue
        fi
    fi

    # --- Step 0b: put every read on the same strand. ---
    # PacBio CCS names carry the strand as the final path element (.../ccs/0 =
    # forward, .../ccs/1 = reverse), so reads arrive in a mix of orientations.
    # Verified against HXB2 with minimap2 across all four local samples: 77/77
    # reads ending in /0 map to '+' and 87/87 ending in /1 map to '-', with no
    # exceptions, and after this step every mapped read is on '+'. Note the
    # fwd/rev WORD some names also carry contradicts the alignment and is
    # deliberately ignored -- see the orient_forward function above.
    #
    # Doing this before the comparison steps means dedup, MSA and motif mapping
    # never have to reason about strand. It also removes the need for the
    # reverse-complement-aware dedup rule: once everything is forward, a read
    # and its former reverse complement are byte-identical, so the *_exact and
    # *_bothstrands variants below should now agree -- which makes their
    # agreement a useful check that this step did what it claims.
    ORIENTED="${PROV_DIR}/${SRR}.provirus.oriented.${EXT}.gz"  # all-forward proviral cores
    ORREPORT="${PROV_DIR}/${SRR}.orient_report.tsv"      # per-read strand tag + action taken
    OTIME="${RESULTS_DIR}/orient_${SRR}.time"            # timing file for the orientation step
    OLOG="${RESULTS_DIR}/orient_${SRR}.log"              # log for the orientation step
    echo "=== orient reads to forward strand on ${SRR} ==="  # progress marker
    measure_and_run "${OTIME}" -- \
        bash -c 'orient_forward "$@"' _ \
        "${STRIPPED}" "${ORIENTED}" "${ORREPORT}" > "${OLOG}" 2>&1
    OEXIT=$?                                         # capture the orienter's exit status
    # set WALLCLOCK_SEC / PEAK_RSS_MB from the .time file
    parse_time_metrics "${OTIME}"
    N_OR=$(count_records "${ORIENTED}")              # records that came through the orientation
    OVALID=0; OMETRIC="n/a"                          # assume invalid until proven otherwise
    if [ -s "${ORIENTED}" ] && [ "${N_OR}" -gt 0 ]; then
        # how many were flipped vs left alone, straight from the report
        N_FLIP=$(awk -F'\t' 'NR>1 && $3=="revcomp"' "${ORREPORT}" 2>/dev/null | wc -l)
        N_UNTAG=$(awk -F'\t' 'NR>1 && $3=="kept_untagged"' "${ORREPORT}" 2>/dev/null | wc -l)
        # mark valid and record what the orientation actually did
        OVALID=1; OMETRIC="${N_OR} reads all forward (${N_FLIP} reverse-complemented, ${N_UNTAG} untagged)"
    fi
    # write the orientation step's row to summary.tsv
    append_summary_row "download_qc_pacbio" "orient_fwd" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${OEXIT}" "${OVALID}" "${OMETRIC}"
    if [ "${OVALID}" = "1" ]; then                   # orientation succeeded...
        STRIPPED="${ORIENTED}"                       # ...so everything downstream uses the oriented reads
    else
        # fall back to the un-oriented cores rather than losing the sample entirely
        echo "WARNING: orientation failed for ${SRR}, continuing with un-oriented reads." >&2
    fi

    # --- Choose this arm's length bar now that the format is known. ---
    if [ -n "${LEN_MIN}" ]; then                     # an explicit override always wins
        LEN_CUT="${LEN_MIN}"                         # use exactly what the caller asked for
    elif [ "${FORMAT}" = "fasta" ]; then
        # SMRTcap cores: 500 keeps all 165 reads (shortest core is 763bp) while
        # still excluding anything too short to carry a complete U3
        LEN_CUT=500
    else
        LEN_CUT=1000                                 # raw unmasked HiFi reads
    fi
    echo "len_min: ${LEN_CUT}"                       # make the effective bar visible in the log

    # --- Step 1: NanoPlot before and after stripping. ---
    # The pair is the point: the pre-strip length distribution includes the host
    # mask, the post-strip one is the provirus, and the difference reads out how
    # much of each read the SMRTcap mask covered.
    echo "=== NanoPlot (pre-strip) on ${SRR} ==="    # progress marker in the log
    run_nanoplot "${READS}" "${RESULTS_DIR}/nanoplot_pre_out/${SRR}" "${SRR}_pre_" "${FORMAT}"
    echo "=== NanoPlot (post-strip) on ${SRR} ==="   # progress marker in the log
    run_nanoplot "${STRIPPED}" "${RESULTS_DIR}/nanoplot_post_out/${SRR}" "${SRR}_post_" "${FORMAT}"

    # --- Step 2: filter comparison, tool set chosen by format. ---
    if [ "${FORMAT}" = "fastq" ]; then
        FILTER_TOOLS="nanofilt fastp chopper"        # quality+length filters, FASTQ only
    else
        FILTER_TOOLS="seqkit awk"                    # length-only filters that accept FASTA
        # state plainly why the three FASTQ filters are absent, rather than leaving gaps
        append_na_row "nanofilt" "${SRR}" "FASTA input: NanoFilt requires FASTQ (rejects records not starting with '@')"
        append_na_row "chopper"  "${SRR}" "FASTA input: chopper requires FASTQ (fails to parse the record)"
        append_na_row "fastp"    "${SRR}" "FASTA input: fastp exits 0 but writes an empty file"
    fi

    # run each applicable filter on the same stripped input for a head-to-head comparison
    for TOOL in ${FILTER_TOOLS}; do
        # chopper is optional...
        if [ "${TOOL}" = "chopper" ] && ! command -v chopper >/dev/null 2>&1; then
            # ...note its absence...
            echo "NOTE: chopper not installed, skipping (see HIV_U3analysis_env.yml)." >&2
            append_na_row "chopper" "${SRR}" "binary not installed"
            continue                                 # ...and skip it if the binary isn't on PATH
        fi
        OUTDIR="${RESULTS_DIR}/${TOOL}_out"          # per-tool output dir (e.g. fastp_out/)
        mkdir -p "${OUTDIR}"                         # create it if needed
        # this tool's filtered-reads output for this sample
        FILT="${OUTDIR}/${SRR}.filtered.${EXT}.gz"
        # file where measure_and_run records wallclock/RSS
        TIMELOG="${RESULTS_DIR}/${TOOL}_${SRR}.time"
        LOG="${RESULTS_DIR}/${TOOL}_${SRR}.log"      # captured stdout+stderr of the tool

        echo "=== ${TOOL} filter on ${SRR} ==="      # progress marker in the log
        # dispatch to the right command per tool (same thresholds, tool-specific syntax)
        case "${TOOL}" in
            nanofilt)
                # NanoFilt reads stdin, so decompress in, filter, gzip out
                measure_and_run "${TIMELOG}" -- \
                    bash -c 'zcat -f "$1" | NanoFilt -q "$2" -l "$3" | gzip > "$4"' _ \
                    "${STRIPPED}" "${Q_MIN}" "${LEN_CUT}" "${FILT}" > "${LOG}" 2>&1 ;;
            fastp)
                # fastp works on the .gz directly; adapter trimming off (HiFi is adapter-clean);
                # emit JSON/HTML reports
                measure_and_run "${TIMELOG}" -- \
                    fastp -i "${STRIPPED}" -o "${FILT}" \
                        --disable_adapter_trimming \
                        --average_qual "${Q_MIN}" --length_required "${LEN_CUT}" \
                        --json "${OUTDIR}/${SRR}_fastp.json" --html "${OUTDIR}/${SRR}_fastp.html" \
                        --thread "${THREADS}" > "${LOG}" 2>&1 ;;
            chopper)
                # chopper is also stdin/stdout, so same decompress|filter|gzip pattern as NanoFilt
                measure_and_run "${TIMELOG}" -- \
                    bash -c 'zcat -f "$1" | chopper -q "$2" -l "$3" --threads "$4" | gzip > "$5"' _ \
                    "${STRIPPED}" "${Q_MIN}" "${LEN_CUT}" "${THREADS}" "${FILT}" > "${LOG}" 2>&1 ;;
            seqkit)
                # seqkit seq -m is the length filter that does accept FASTA
                measure_and_run "${TIMELOG}" -- \
                    bash -c 'seqkit seq -m "$2" "$1" | gzip > "$3"' _ \
                    "${STRIPPED}" "${LEN_CUT}" "${FILT}" > "${LOG}" 2>&1 ;;
            awk)
                # pure-awk length filter, the independent second implementation
                measure_and_run "${TIMELOG}" -- \
                    bash -c 'filter_len_awk "$@"' _ \
                    "${STRIPPED}" "${FILT}" "${LEN_CUT}" > "${LOG}" 2>&1 ;;
        esac
        # capture the tool's exit status before $? is overwritten
        EXIT_CODE=$?
        # set WALLCLOCK_SEC / PEAK_RSS_MB globals from the .time file
        parse_time_metrics "${TIMELOG}"

        N_FILT=$(count_records "${FILT}")            # how many reads survived this filter
        VALID=0; METRIC="n/a"                        # assume invalid until proven otherwise
        # valid = non-empty, not-corrupt gzip, and at least one read
        if [ -s "${FILT}" ] && gzip -t "${FILT}" 2>/dev/null && [ "${N_FILT}" -gt 0 ]; then
            # the threshold description depends on whether quality was available
            if [ "${FORMAT}" = "fastq" ]; then
                # mark valid and record the length+quality thresholds applied
                VALID=1; METRIC="${N_FILT} reads passed (>=${LEN_CUT}bp, >=Q${Q_MIN})"
            else
                # FASTA: length only, and say so explicitly
                VALID=1; METRIC="${N_FILT} reads passed (>=${LEN_CUT}bp; no quality filter, FASTA)"
            fi
        fi
        # write this filter's row to summary.tsv
        append_summary_row "download_qc_pacbio" "${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${EXIT_CODE}" "${VALID}" "${METRIC}"

        # --- Step 3: Kraken2 host/bacterial removal on this filter's output. ---
        if [ ! -d "${KRAKEN2_DB}" ]; then            # if the Kraken2 DB isn't present...
            # ...note it and skip host removal
            echo "NOTE: Kraken2 DB ${KRAKEN2_DB} absent, skipping host removal for ${TOOL}/${SRR}." >&2
            append_na_row "kraken2_after_${TOOL}" "${SRR}" "Kraken2 DB not present at ${KRAKEN2_DB}"
        # only run Kraken2 when the filter actually produced usable reads
        elif [ "${VALID}" = "1" ]; then
            KOUT="${RESULTS_DIR}/kraken2_${TOOL}_out"          # per-filter Kraken2 output dir
            KTIME="${RESULTS_DIR}/kraken2_${TOOL}_${SRR}.time" # timing file for this Kraken2 run
            KLOG="${RESULTS_DIR}/kraken2_${TOOL}_${SRR}.log"   # log for this Kraken2 run
            echo "=== Kraken2 host removal after ${TOOL} on ${SRR} ==="  # progress marker
            # single-end Kraken2 step: classify reads and drop host/bacterial ones
            measure_and_run "${KTIME}" -- \
                bash -c 'kraken2_filter_se "$@"' _ \
                "${FILT}" "${KOUT}" "${SRR}" "${KRAKEN2_DB}" > "${KLOG}" 2>&1
            KEXIT=$?                                  # capture Kraken2's exit status
            # refresh WALLCLOCK_SEC / PEAK_RSS_MB from the Kraken2 timing file
            parse_time_metrics "${KTIME}"
            # the host-removed reads the wrapper produces, named after the input format
            KFILT="${KOUT}/${SRR}.kraken_filtered.${FORMAT}.gz"
            N_KRAK=$(count_records "${KFILT}")       # reads remaining after host removal
            KVALID=0; KMETRIC="n/a"                  # default to invalid until checked
            # valid if the output exists and has reads
            [ -s "${KFILT}" ] && [ "${N_KRAK}" -gt 0 ] && { KVALID=1; KMETRIC="${N_KRAK}/${N_FILT} reads retained after host removal"; }
            # record the Kraken2 row
            append_summary_row "download_qc_pacbio" "kraken2_after_${TOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${KEXIT}" "${KVALID}" "${KMETRIC}"
        fi
    done

    # --- Step 4: dedup comparison on the cleanest available reads. ---
    # Preferred input is the Kraken2-cleaned output of the format's primary
    # filter, falling back to that filter's own output when Kraken2 was skipped.
    # seqkit and awk are each run under BOTH matching rules, because they are
    # genuinely different questions and the tools do not default to the same one:
    #   *_exact        byte-identical sequences only        (seqkit -s -P)
    #   *_bothstrands  a sequence and its reverse complement (seqkit -s, the default)
    # The distinction is real for SMRTcap data -- the same provirus read in the
    # opposite orientation appears as a reverse complement, and the read names
    # say so (.../ccs/fwd/1, .../ccs/rev/0). Comparing seqkit's default against
    # a plain exact-match awk would have compared two different rules and made
    # awk look like it under-collapsed; pairing them by rule keeps it honest.
    DEDUP_TOOLS="seqkit_exact seqkit_bothstrands awk_exact awk_bothstrands"
    if [ "${FORMAT}" = "fastq" ]; then
        PRIMARY="fastp"                              # fastp is the FASTQ arm's reference filter
        DEDUP_TOOLS="fastp ${DEDUP_TOOLS}"           # fastp --dedup only joins in on FASTQ
    else
        PRIMARY="seqkit"                             # seqkit is the FASTA arm's reference filter
        # record why fastp is missing from the FASTA dedup comparison
        append_na_row "dedup_fastp" "${SRR}" "FASTA input: fastp --dedup exits 0 but writes an empty file"
    fi
    # preferred dedup input: primary-filtered + host-removed reads
    CLEAN="${RESULTS_DIR}/kraken2_${PRIMARY}_out/${SRR}.kraken_filtered.${FORMAT}.gz"
    # fall back to just the primary-filtered reads if Kraken2 was skipped
    [ -s "${CLEAN}" ] || CLEAN="${RESULTS_DIR}/${PRIMARY}_out/${SRR}.filtered.${EXT}.gz"
    if [ -s "${CLEAN}" ]; then                       # only dedup if we actually have a clean input
        # read count before dedup, for the "N duplicates removed" metric
        N_BEFORE=$(count_records "${CLEAN}")
        # compare the two dedup tools on identical input
        for DTOOL in ${DEDUP_TOOLS}; do
            # per-tool dedup output dir (created if needed)
            DOUT="${RESULTS_DIR}/dedup_${DTOOL}_out"; mkdir -p "${DOUT}"
            DFILE="${DOUT}/${SRR}.dedup.${EXT}.gz"   # deduplicated reads output
            DTIME="${RESULTS_DIR}/dedup_${DTOOL}_${SRR}.time"  # timing file for this dedup run
            DLOG="${RESULTS_DIR}/dedup_${DTOOL}_${SRR}.log"    # log for this dedup run
            echo "=== dedup (${DTOOL}) on ${SRR} ==="  # progress marker
            case "${DTOOL}" in                        # tool-specific dedup command
                fastp)
                    # fastp in dedup-only mode: all other filtering disabled so it only removes
                    # duplicates
                    measure_and_run "${DTIME}" -- \
                        fastp -i "${CLEAN}" -o "${DFILE}" --dedup \
                            --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering \
                            --json "${DOUT}/${SRR}_dedup.json" --html "${DOUT}/${SRR}_dedup.html" \
                            --thread "${THREADS}" > "${DLOG}" 2>&1 ;;
                seqkit_exact)
                    # -P restricts the comparison to the positive strand, i.e. exact matches only
                    measure_and_run "${DTIME}" -- \
                        bash -c 'seqkit rmdup -s -P "$1" -o "$2"' _ "${CLEAN}" "${DFILE}" > "${DLOG}" 2>&1 ;;
                seqkit_bothstrands)
                    # seqkit's default: a sequence and its reverse complement count as duplicates
                    measure_and_run "${DTIME}" -- \
                        bash -c 'seqkit rmdup -s "$1" -o "$2"' _ "${CLEAN}" "${DFILE}" > "${DLOG}" 2>&1 ;;
                awk_exact)
                    # pure-awk equivalent of seqkit rmdup -s -P
                    measure_and_run "${DTIME}" -- \
                        bash -c 'rmdup_awk "$@"' _ "${CLEAN}" "${DFILE}" exact > "${DLOG}" 2>&1 ;;
                awk_bothstrands)
                    # pure-awk equivalent of seqkit rmdup -s
                    measure_and_run "${DTIME}" -- \
                        bash -c 'rmdup_awk "$@"' _ "${CLEAN}" "${DFILE}" bothstrands > "${DLOG}" 2>&1 ;;
            esac
            DEXIT=$?                                  # capture the dedup tool's exit status
            parse_time_metrics "${DTIME}"            # refresh timing/RSS from this dedup run
            N_AFTER=$(count_records "${DFILE}")      # read count after dedup
            DVALID=0; DMETRIC="n/a"                  # default to invalid until checked
            # valid if output has reads; report kept vs removed
            [ -s "${DFILE}" ] && [ "${N_AFTER}" -gt 0 ] && { DVALID=1; DMETRIC="${N_AFTER}/${N_BEFORE} reads kept ($((N_BEFORE-N_AFTER)) duplicates removed)"; }
            # record the dedup row
            append_summary_row "download_qc_pacbio" "dedup_${DTOOL}" "${SRR}" "${WALLCLOCK_SEC}" "${PEAK_RSS_MB}" "${DEXIT}" "${DVALID}" "${DMETRIC}"
        done
    fi
done

# final confirmation pointing the user at the results table
echo ""
echo "Done. See ${SUMMARY_TSV}"
