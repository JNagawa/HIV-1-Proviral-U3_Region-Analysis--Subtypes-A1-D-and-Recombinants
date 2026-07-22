#!/bin/bash
#SBATCH --job-name=compare_step
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --time=14:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# Generic sbatch wrapper for any scripts/<step>/illumina/compare_*.sh
# tool-comparison harness script, so none of them need their own #SBATCH
# headers (they're also meant to be runnable directly for quick local
# testing, per their own "Usage: ./compare_X.sh" docstrings).
# Usage: sbatch scripts/utils/run_comparison_step.slurm.sh <path/to/compare_script.sh>
set -euo pipefail
TARGET_SCRIPT="$1"

CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"
if [ -f "$CONDA_SH" ]; then
    source "$CONDA_SH"
else
    source "$(conda info --base)/etc/profile.d/conda.sh"
fi
conda activate HIV_U3analysis

export THREADS="${SLURM_CPUS_PER_TASK:-8}"
cd /etc/ace-data/home/jnagawa/Internship
bash "${TARGET_SCRIPT}"
