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
set -euo pipefail

cd /etc/ace-data/home/jnagawa/Internship
W=scripts/utils/run_comparison_step.slurm.sh
ENVJOB=scripts/utils/setup_pacbio_env.slurm.sh

# Sanity: HiFi data present?
if ! ls data/raw/pacbio/*.fastq.gz >/dev/null 2>&1; then
    echo "ERROR: no data/raw/pacbio/*.fastq.gz -- run/await scripts/download_qc/pacbio/download_pacbio.slurm.sh first." >&2
    exit 1
fi

ENV=$(sbatch --parsable "${ENVJOB}")
echo "env-setup           = ${ENV}"

QC=$(sbatch --parsable --job-name=pb_download_qc --dependency=afterok:${ENV} \
     "${W}" scripts/download_qc/pacbio/download_qc_pacbio.sh)
echo "download_qc (afterok:${ENV}) = ${QC}"

EXT=$(sbatch --parsable --job-name=pb_extraction --dependency=afterok:${QC} \
      "${W}" scripts/proviral_extraction/pacbio/extract_provirus_pacbio.sh)
echo "extraction (afterok:${QC}) = ${EXT}"

ASM=$(sbatch --parsable --job-name=pb_assembly --time=24:00:00 --dependency=afterok:${EXT} \
      "${W}" scripts/assembly/pacbio/compare_assembly_pacbio.sh)
echo "assembly (afterok:${EXT}) = ${ASM}"

MSA=$(sbatch --parsable --job-name=pb_msa --dependency=afterok:${ASM} \
      "${W}" scripts/msa/pacbio/compare_msa_pacbio.sh)
echo "msa (afterok:${ASM}) = ${MSA}"

BIO=$(sbatch --parsable --job-name=pb_biofilt --dependency=afterok:${ASM} \
      "${W}" scripts/biological_filtering/pacbio/compare_biological_filtering_pacbio.sh)
echo "biofilt (afterok:${ASM}) = ${BIO}"

SUB=$(sbatch --parsable --job-name=pb_subtyping --dependency=afterok:${MSA} \
      "${W}" scripts/subtyping/pacbio/compare_subtyping_pacbio.sh)
echo "subtyping (afterok:${MSA}) = ${SUB}"

MOT=$(sbatch --parsable --job-name=pb_motif --dependency=afterok:${MSA} \
      "${W}" scripts/motif_mapping/pacbio/compare_motif_mapping_pacbio.sh)
echo "motif (afterok:${MSA}) = ${MOT}"

echo "${ENV} ${QC} ${EXT} ${ASM} ${MSA} ${BIO} ${SUB} ${MOT}"
echo "--- PacBio chain submitted ---"
