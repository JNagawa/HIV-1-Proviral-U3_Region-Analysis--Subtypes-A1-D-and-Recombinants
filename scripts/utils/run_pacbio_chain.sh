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
set -euo pipefail                                    # -e abort on error, -u on unset vars, pipefail on any failed pipe stage

cd /etc/ace-data/home/jnagawa/Internship             # run from the repo root so relative script paths resolve
W=scripts/utils/run_comparison_step.slurm.sh         # generic sbatch wrapper used to launch every comparison step
ENVJOB=scripts/utils/setup_pacbio_env.slurm.sh       # the env-setup job that installs hifiasm/chopper/bc

# Sanity: HiFi data present?
if ! ls data/raw/pacbio/*.fastq.gz >/dev/null 2>&1; then  # refuse to submit if the raw HiFi data isn't downloaded yet
    echo "ERROR: no data/raw/pacbio/*.fastq.gz -- run/await scripts/download_qc/pacbio/download_pacbio.slurm.sh first." >&2  # tell the user what to run first
    exit 1
fi

ENV=$(sbatch --parsable "${ENVJOB}")                 # submit env-setup first; --parsable captures just the job id
echo "env-setup           = ${ENV}"                  # print its job id

QC=$(sbatch --parsable --job-name=pb_download_qc --dependency=afterok:${ENV} \
     "${W}" scripts/download_qc/pacbio/download_qc_pacbio.sh)  # QC step, runs only after env-setup succeeds
echo "download_qc (afterok:${ENV}) = ${QC}"          # print its job id and dependency

EXT=$(sbatch --parsable --job-name=pb_extraction --dependency=afterok:${QC} \
      "${W}" scripts/proviral_extraction/pacbio/extract_provirus_pacbio.sh)  # extraction, after QC succeeds
echo "extraction (afterok:${QC}) = ${EXT}"           # print its job id and dependency

ASM=$(sbatch --parsable --job-name=pb_assembly --time=24:00:00 --dependency=afterok:${EXT} \
      "${W}" scripts/assembly/pacbio/compare_assembly_pacbio.sh)  # assembly (longer walltime), after extraction succeeds
echo "assembly (afterok:${EXT}) = ${ASM}"            # print its job id and dependency

MSA=$(sbatch --parsable --job-name=pb_msa --dependency=afterok:${ASM} \
      "${W}" scripts/msa/pacbio/compare_msa_pacbio.sh)  # multiple-sequence alignment, after assembly succeeds
echo "msa (afterok:${ASM}) = ${MSA}"                 # print its job id and dependency

BIO=$(sbatch --parsable --job-name=pb_biofilt --dependency=afterok:${ASM} \
      "${W}" scripts/biological_filtering/pacbio/compare_biological_filtering_pacbio.sh)  # biological filtering, also branches off assembly
echo "biofilt (afterok:${ASM}) = ${BIO}"             # print its job id and dependency

SUB=$(sbatch --parsable --job-name=pb_subtyping --dependency=afterok:${MSA} \
      "${W}" scripts/subtyping/pacbio/compare_subtyping_pacbio.sh)  # subtyping, after MSA succeeds
echo "subtyping (afterok:${MSA}) = ${SUB}"           # print its job id and dependency

MOT=$(sbatch --parsable --job-name=pb_motif --dependency=afterok:${MSA} \
      "${W}" scripts/motif_mapping/pacbio/compare_motif_mapping_pacbio.sh)  # motif mapping, also branches off MSA
echo "motif (afterok:${MSA}) = ${MOT}"               # print its job id and dependency

echo "${ENV} ${QC} ${EXT} ${ASM} ${MSA} ${BIO} ${SUB} ${MOT}"  # one line with all submitted job ids for easy scancel/squeue
echo "--- PacBio chain submitted ---"                # done marker
