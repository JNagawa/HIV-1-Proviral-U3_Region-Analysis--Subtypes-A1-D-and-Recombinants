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
# -e abort on error, -u on unset vars, pipefail on any failed pipe stage
set -euo pipefail
# the comparison script to run, passed as the sbatch argument
TARGET_SCRIPT="$1"

CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"   # expected conda init script location
if [ -f "$CONDA_SH" ]; then                          # prefer the known miniconda path if present
    # load conda so `conda activate` works in this non-interactive shell
    source "$CONDA_SH"
else
    # fall back to wherever conda reports its base
    source "$(conda info --base)/etc/profile.d/conda.sh"
fi
# activate the shared analysis env with all the tools
conda activate HIV_U3analysis

# expose the allocated core count so tools thread correctly (default 8)
export THREADS="${SLURM_CPUS_PER_TASK:-8}"
# run from the repo root so the script's relative paths resolve
cd /etc/ace-data/home/jnagawa/Internship
bash "${TARGET_SCRIPT}"                              # execute the requested comparison step
