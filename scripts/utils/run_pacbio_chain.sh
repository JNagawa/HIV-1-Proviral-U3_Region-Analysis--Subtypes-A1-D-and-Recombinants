#!/bin/bash
# One-shot launcher for the whole PacBio (HIV-SMRTcap) comparison pipeline.
# Submits every step as an sbatch job wired together with afterok dependencies
# (same pattern used for the Illumina re-run), so each step only starts if its
# prerequisite succeeded. Safe to run once (a) the data download job has
# finished (data/raw/pacbio/*.fastq.gz present) and (b) the Illumina
# comparison chain has finished (so installing hifiasm/chopper into the shared
# env can't disrupt live jobs).
#
# Dependency graph:
#   env-setup ─┬─> download_qc ─> extraction ─> assembly ─┬─> msa ─┬─> subtyping
#              │                                          │        └─> motif_mapping
#              └────────────────────────────────────────>└─> biological_filtering
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

# multiple-sequence alignment, after assembly succeeds
MSA=$(sbatch --parsable --job-name=pb_msa --dependency=afterok:${ASM} \
      "${W}" scripts/msa/pacbio/compare_msa_pacbio.sh)
echo "msa (afterok:${ASM}) = ${MSA}"                 # print its job id and dependency

# biological filtering, also branches off assembly
BIO=$(sbatch --parsable --job-name=pb_biofilt --dependency=afterok:${ASM} \
      "${W}" scripts/biological_filtering/pacbio/compare_biological_filtering_pacbio.sh)
echo "biofilt (afterok:${ASM}) = ${BIO}"             # print its job id and dependency

SUB=$(sbatch --parsable --job-name=pb_subtyping --dependency=afterok:${MSA} \
      "${W}" scripts/subtyping/pacbio/compare_subtyping_pacbio.sh)  # subtyping, after MSA succeeds
echo "subtyping (afterok:${MSA}) = ${SUB}"           # print its job id and dependency

# motif mapping, also branches off MSA
MOT=$(sbatch --parsable --job-name=pb_motif --dependency=afterok:${MSA} \
      "${W}" scripts/motif_mapping/pacbio/compare_motif_mapping_pacbio.sh)
echo "motif (afterok:${MSA}) = ${MOT}"               # print its job id and dependency

# one line with all submitted job ids for easy scancel/squeue
echo "${ENV} ${QC} ${EXT} ${ASM} ${MSA} ${BIO} ${SUB} ${MOT}"
echo "--- PacBio chain submitted ---"                # done marker
