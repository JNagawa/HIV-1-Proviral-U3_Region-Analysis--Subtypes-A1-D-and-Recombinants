#!/bin/bash
# SHIVER assembly (the review's top pick, to avoid reference bias on
# divergent A1/D/recombinant genomes). SHIVER is a 3-step chain, not a single
# tool call: assemble contigs (SPAdes) -> shiver_align_contigs.sh ->
# shiver_map_reads.sh. SPAdes's own cost is broken out into a separate log
# so it isn't double-counted against SHIVER's own runtime/RSS when SPAdes is
# ALSO being benchmarked as its own competing assembler in this same stage.
#
# shiver_align_contigs.sh / shiver_map_reads.sh write their outputs relative
# to the current working directory (no output-dir argument), so this script
# cd's into OUTDIR before calling them.
#
# Usage: run_shiver.sh <SRR> <R1.fastq.gz> <R2.fastq.gz> <OUTDIR>
# -u errors on unset vars, pipefail fails a pipe if any stage fails
set -uo pipefail
# positional args: sample id, paired reads, output dir
SRR="$1" R1="$2" R2="$3" OUTDIR="$4"

# absolute path of this script's own dir, so paths work regardless of launch location
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
# dir holding SHIVER's config.sh and reference alignment
SETUP_DIR="${STAGE_DIR}/shiver_setup"
# repo top-level (three dirs up), the base for tool paths below
REPO_ROOT="$(cd "${STAGE_DIR}/../../.." && pwd)"
SHIVER_BIN="${REPO_ROOT}/scripts/tools/shiver/bin"   # dir with the SHIVER executable scripts

# LANL multi-subtype reference alignment panel
REF_ALIGNMENT="${SETUP_DIR}/HIV1_COM_ref_alignment.fasta"
# the panel is a manual one-time download, so verify it exists
if [ ! -s "${REF_ALIGNMENT}" ]; then
    # tell the user how to obtain it
    echo "ERROR: ${REF_ALIGNMENT} not found. See shiver_setup/SOURCE.md for" \
         "the one-time manual LANL download step this requires." >&2
    exit 1                                           # cannot proceed without the reference panel
fi

# adapter sequences SHIVER trims from reads
ADAPTERS="${REPO_ROOT}/scripts/tools/shiver/data/example_input/adapters_Illumina.fasta"
# PCR primer sequences SHIVER trims from reads
PRIMERS="${REPO_ROOT}/scripts/tools/shiver/data/example_input/primers_GallEtAl2012.fasta"
INIT_DIR="${OUTDIR}/shiver_init"                     # dir for SHIVER's one-time init output
mkdir -p "${OUTDIR}"                                 # ensure the output dir exists

# Resolve R1/R2 to absolute paths before we cd, since callers may pass
# relative paths.
# canonicalize R1 to an absolute path (survives the later cd)
R1="$(cd "$(dirname "${R1}")" && pwd)/$(basename "${R1}")"
# canonicalize R2 to an absolute path (survives the later cd)
R2="$(cd "$(dirname "${R2}")" && pwd)/$(basename "${R2}")"

# shiver_init.sh only needs to run once (its OutDir must not pre-exist)
# only initialise if it hasn't been done yet (dir must not pre-exist)
if [ ! -d "${INIT_DIR}" ]; then
    # build SHIVER's initialisation dir from the reference/adapters/primers
    "${SHIVER_BIN}/shiver_init.sh" "${INIT_DIR}" "${SETUP_DIR}/config.sh" \
        "${REF_ALIGNMENT}" "${ADAPTERS}" "${PRIMERS}" \
        > "${OUTDIR}/shiver_init.log" 2>&1
    if [ $? -ne 0 ]; then                            # if init failed...
        # ...point at the log...
        echo "ERROR: shiver_init.sh failed, see ${OUTDIR}/shiver_init.log" >&2
        exit 1                                       # ...and abort
    fi
fi

# SHIVER writes to cwd, so enter this sample's output dir
cd "${OUTDIR}" || exit 1

# Step 1 of the chain: assemble contigs with SPAdes (cost logged separately).
# separate SPAdes dir so its cost is logged apart from SHIVER's
CONTIGS_DIR="${SRR}_spades_for_shiver"
spades.py --careful -1 "${R1}" -2 "${R2}" -o "${CONTIGS_DIR}" -t "${THREADS:-4}" \
    > "${SRR}_spades_for_shiver.log" 2>&1            # assemble contigs to feed into SHIVER
CONTIGS="${CONTIGS_DIR}/contigs.fasta"               # the contigs SPAdes produced
# abort if no contigs
[ -s "${CONTIGS}" ] || { echo "ERROR: SPAdes (for SHIVER) produced no contigs for ${SRR}" >&2; exit 1; }

# Step 2: align contigs to the reference set, excluding contamination.
# clean/align contigs to the reference panel
"${SHIVER_BIN}/shiver_align_contigs.sh" "${INIT_DIR}" "${SETUP_DIR}/config.sh" \
    "${CONTIGS}" "${SRR}" > "${SRR}_align_contigs.log" 2>&1
if [ $? -ne 0 ]; then                                # if alignment failed...
    # ...point at the log...
    echo "ERROR: shiver_align_contigs.sh failed, see ${OUTDIR}/${SRR}_align_contigs.log" >&2
    exit 1                                           # ...and abort
fi
BLAST_FILE="${SRR}.blast"                            # blast-hits output from the align step
# cleaned contigs + references output from the align step
CUT_WREFS="${SRR}_cut_wRefs.fasta"
# no hits means contigs aren't HIV
[ -s "${BLAST_FILE}" ] || { echo "ERROR: no blast hits for ${SRR} -- contigs may not be HIV" >&2; exit 1; }

# Step 3: map reads using the aligned contigs to build the final assembly.
# map reads to the sample-specific reference to build the consensus
"${SHIVER_BIN}/shiver_map_reads.sh" "${INIT_DIR}" "${SETUP_DIR}/config.sh" \
    "${CONTIGS}" "${SRR}" "${BLAST_FILE}" "${CUT_WREFS}" "${R1}" "${R2}" \
    > "${SRR}_map_reads.log" 2>&1
if [ $? -ne 0 ]; then                                # if mapping failed...
    # ...point at the log...
    echo "ERROR: shiver_map_reads.sh failed, see ${OUTDIR}/${SRR}_map_reads.log" >&2
    exit 1                                           # ...and abort
fi

# SHIVER names its consensus output "<SID>_remap_consensus_MinCov_*.fasta";
# take the first (default coverage threshold) match.
# glob for SHIVER's consensus output, take the first match
CONSENSUS_SRC=$(compgen -G "${SRR}_remap_consensus_MinCov_*.fasta" | head -1)
# if no usable consensus was produced...
if [ -z "${CONSENSUS_SRC}" ] || [ ! -s "${CONSENSUS_SRC}" ]; then
    echo "ERROR: SHIVER did not produce a consensus FASTA for ${SRR}" >&2  # ...report it...
    exit 1                                           # ...and abort
fi
# copy to the standard consensus name the harness expects
cp "${CONSENSUS_SRC}" "${SRR}_consensus.fasta"
sed -i "1s/.*/>${SRR}/" "${SRR}_consensus.fasta"     # rename the FASTA header to the sample id
