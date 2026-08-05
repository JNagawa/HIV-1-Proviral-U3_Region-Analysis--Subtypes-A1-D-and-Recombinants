#!/bin/bash
# One-shot launcher for the whole PacBio (HIV-SMRTcap) comparison pipeline.
# Submits every step as an sbatch job wired together with afterok dependencies
# (same pattern used for the Illumina re-run), so each step only starts if its
# prerequisite succeeded. Safe to run once (a) the data download job has
# finished (data/raw/pacbio/*.fastq.gz present) and (b) the Illumina
# comparison chain has finished (so installing hifiasm/chopper into the shared
# env can't disrupt live jobs).
#
# Dependency graph (the bracketed part runs once per assembly arm):
#   env-setup ─> download_qc ─> extraction ─> assembly ─┬─> msa ─┬─> subtyping
#                                                       │        └─> motif_mapping
#                                                       └─> biological_filtering
# Usage: bash scripts/utils/run_pacbio_chain.sh
# -e abort on error, -u on unset vars, pipefail on any failed pipe stage
set -euo pipefail

# run from the repo root so relative script paths resolve
cd /etc/ace-data/home/jnagawa/Internship
# generic sbatch wrapper used to launch every comparison step
W=scripts/utils/run_comparison_step.slurm.sh
# the env-setup job that installs hifiasm/chopper/bc
ENVJOB=scripts/utils/setup_pacbio_env.slurm.sh

# Sanity: HiFi data present?
# refuse to submit if the raw HiFi data isn't downloaded yet
if ! ls data/raw/pacbio/*.fastq.gz >/dev/null 2>&1; then
    # tell the user what to run first
    echo "ERROR: no data/raw/pacbio/*.fastq.gz -- run/await scripts/download_qc/pacbio/download_pacbio.slurm.sh first." >&2
    exit 1
fi

# submit env-setup first; --parsable captures just the job id
ENV=$(sbatch --parsable "${ENVJOB}")
echo "env-setup           = ${ENV}"                  # print its job id

# QC step, runs only after env-setup succeeds
QC=$(sbatch --parsable --job-name=pb_download_qc --dependency=afterok:${ENV} \
     "${W}" scripts/download_qc/pacbio/download_qc_pacbio.sh)
echo "download_qc (afterok:${ENV}) = ${QC}"          # print its job id and dependency

# extraction, after QC succeeds
EXT=$(sbatch --parsable --job-name=pb_extraction --dependency=afterok:${QC} \
      "${W}" scripts/proviral_extraction/pacbio/extract_provirus_pacbio.sh)
echo "extraction (afterok:${QC}) = ${EXT}"           # print its job id and dependency

# assembly (longer walltime), after extraction succeeds
ASM=$(sbatch --parsable --job-name=pb_assembly --time=24:00:00 --dependency=afterok:${EXT} \
      "${W}" scripts/assembly/pacbio/compare_assembly_pacbio.sh)
echo "assembly (afterok:${EXT}) = ${ASM}"            # print its job id and dependency

# The assembly step produces two reference-guided consensus sets -- an HXB2
# baseline and a subtype-matched one -- and every downstream stage is run once
# per arm so the effect of the reference choice stays measurable all the way to
# the motif hits. Each arm writes into results/<stage>/pacbio/<arm>/, so the two
# passes never collide. Override with e.g. ASSEMBLY_ARMS="minimap2_consensus".
ARMS="${ASSEMBLY_ARMS:-minimap2_consensus minimap2_bestref}"
# accumulate every submitted job id for the summary line at the end
ALL_IDS="${ENV} ${QC} ${EXT} ${ASM}"

for ARM in ${ARMS}; do
    # multiple-sequence alignment, after assembly succeeds
    MSA=$(sbatch --parsable --job-name="pb_msa_${ARM}" --dependency=afterok:${ASM} \
          --export=ALL,ASSEMBLY_ARM="${ARM}" "${W}" scripts/msa/pacbio/compare_msa_pacbio.sh)
    echo "msa[${ARM}] (afterok:${ASM}) = ${MSA}"      # print its job id and dependency

    # biological filtering, also branches off assembly
    BIO=$(sbatch --parsable --job-name="pb_biofilt_${ARM}" --dependency=afterok:${ASM} \
          --export=ALL,ASSEMBLY_ARM="${ARM}" "${W}" scripts/biological_filtering/pacbio/compare_biological_filtering_pacbio.sh)
    echo "biofilt[${ARM}] (afterok:${ASM}) = ${BIO}"  # print its job id and dependency

    # subtyping, after this arm's MSA succeeds
    SUB=$(sbatch --parsable --job-name="pb_subtyping_${ARM}" --dependency=afterok:${MSA} \
          --export=ALL,ASSEMBLY_ARM="${ARM}" "${W}" scripts/subtyping/pacbio/compare_subtyping_pacbio.sh)
    echo "subtyping[${ARM}] (afterok:${MSA}) = ${SUB}"  # print its job id and dependency

    # motif mapping, also branches off this arm's MSA
    MOT=$(sbatch --parsable --job-name="pb_motif_${ARM}" --dependency=afterok:${MSA} \
          --export=ALL,ASSEMBLY_ARM="${ARM}" "${W}" scripts/motif_mapping/pacbio/compare_motif_mapping_pacbio.sh)
    echo "motif[${ARM}] (afterok:${MSA}) = ${MOT}"    # print its job id and dependency

    ALL_IDS="${ALL_IDS} ${MSA} ${BIO} ${SUB} ${MOT}"  # remember this arm's ids
done

# one line with all submitted job ids for easy scancel/squeue
echo "${ALL_IDS}"
echo "--- PacBio chain submitted ---"                # done marker
