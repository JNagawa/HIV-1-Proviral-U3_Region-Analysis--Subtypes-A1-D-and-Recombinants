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
set -uo pipefail
SRR="$1" R1="$2" R2="$3" OUTDIR="$4"

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
SETUP_DIR="${STAGE_DIR}/shiver_setup"
REPO_ROOT="$(cd "${STAGE_DIR}/../.." && pwd)"
SHIVER_BIN="${REPO_ROOT}/scripts/tools/shiver/bin"

REF_ALIGNMENT="${SETUP_DIR}/HIV1_COM_ref_alignment.fasta"
if [ ! -s "${REF_ALIGNMENT}" ]; then
    echo "ERROR: ${REF_ALIGNMENT} not found. See shiver_setup/SOURCE.md for" \
         "the one-time manual LANL download step this requires." >&2
    exit 1
fi

ADAPTERS="${REPO_ROOT}/scripts/tools/shiver/data/example_input/adapters_Illumina.fasta"
PRIMERS="${REPO_ROOT}/scripts/tools/shiver/data/example_input/primers_GallEtAl2012.fasta"
INIT_DIR="${OUTDIR}/shiver_init"
mkdir -p "${OUTDIR}"

# Resolve R1/R2 to absolute paths before we cd, since callers may pass
# relative paths.
R1="$(cd "$(dirname "${R1}")" && pwd)/$(basename "${R1}")"
R2="$(cd "$(dirname "${R2}")" && pwd)/$(basename "${R2}")"

# shiver_init.sh only needs to run once (its OutDir must not pre-exist)
if [ ! -d "${INIT_DIR}" ]; then
    "${SHIVER_BIN}/shiver_init.sh" "${INIT_DIR}" "${SETUP_DIR}/config.sh" \
        "${REF_ALIGNMENT}" "${ADAPTERS}" "${PRIMERS}" \
        > "${OUTDIR}/shiver_init.log" 2>&1
    if [ $? -ne 0 ]; then
        echo "ERROR: shiver_init.sh failed, see ${OUTDIR}/shiver_init.log" >&2
        exit 1
    fi
fi

cd "${OUTDIR}" || exit 1

# Step 1 of the chain: assemble contigs with SPAdes (cost logged separately).
CONTIGS_DIR="${SRR}_spades_for_shiver"
spades.py --careful -1 "${R1}" -2 "${R2}" -o "${CONTIGS_DIR}" -t "${THREADS:-4}" \
    > "${SRR}_spades_for_shiver.log" 2>&1
CONTIGS="${CONTIGS_DIR}/contigs.fasta"
[ -s "${CONTIGS}" ] || { echo "ERROR: SPAdes (for SHIVER) produced no contigs for ${SRR}" >&2; exit 1; }

# Step 2: align contigs to the reference set, excluding contamination.
"${SHIVER_BIN}/shiver_align_contigs.sh" "${INIT_DIR}" "${SETUP_DIR}/config.sh" \
    "${CONTIGS}" "${SRR}" > "${SRR}_align_contigs.log" 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: shiver_align_contigs.sh failed, see ${OUTDIR}/${SRR}_align_contigs.log" >&2
    exit 1
fi
BLAST_FILE="${SRR}.blast"
CUT_WREFS="${SRR}_cut_wRefs.fasta"
[ -s "${BLAST_FILE}" ] || { echo "ERROR: no blast hits for ${SRR} -- contigs may not be HIV" >&2; exit 1; }

# Step 3: map reads using the aligned contigs to build the final assembly.
"${SHIVER_BIN}/shiver_map_reads.sh" "${INIT_DIR}" "${SETUP_DIR}/config.sh" \
    "${CONTIGS}" "${SRR}" "${BLAST_FILE}" "${CUT_WREFS}" "${R1}" "${R2}" \
    > "${SRR}_map_reads.log" 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: shiver_map_reads.sh failed, see ${OUTDIR}/${SRR}_map_reads.log" >&2
    exit 1
fi

# SHIVER names its consensus output "<SID>_remap_consensus_MinCov_*.fasta";
# take the first (default coverage threshold) match.
CONSENSUS_SRC=$(compgen -G "${SRR}_remap_consensus_MinCov_*.fasta" | head -1)
if [ -z "${CONSENSUS_SRC}" ] || [ ! -s "${CONSENSUS_SRC}" ]; then
    echo "ERROR: SHIVER did not produce a consensus FASTA for ${SRR}" >&2
    exit 1
fi
cp "${CONSENSUS_SRC}" "${SRR}_consensus.fasta"
sed -i "1s/.*/>${SRR}/" "${SRR}_consensus.fasta"
