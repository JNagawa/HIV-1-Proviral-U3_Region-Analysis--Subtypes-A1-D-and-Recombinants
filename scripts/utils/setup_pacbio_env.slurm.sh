#!/bin/bash
#SBATCH --job-name=pacbio_env_setup
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# Install the long-read tools the PacBio arm needs (hifiasm, chopper) into the
# HIV_U3analysis env, plus `bc` (missing on the compute nodes, which blocks
# SHIVER's shiver_init on both the Illumina and PacBio arms). All pinned in
# HIV_U3analysis_env.yml. Run as its own sbatch job (not on the login node)
# AND only once the Illumina comparison chain has finished -- installing into
# a shared env while those jobs are live could swap libraries mid-run.
# run_pacbio_chain.sh submits this first and makes every PacBio step depend on
# it (afterok).
# Usage: sbatch scripts/utils/setup_pacbio_env.slurm.sh
set -euo pipefail                                    # -e abort on error, -u on unset vars, pipefail on any failed pipe stage

CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"   # expected conda init script location
if [ -f "$CONDA_SH" ]; then source "$CONDA_SH"; else source "$(conda info --base)/etc/profile.d/conda.sh"; fi  # load conda from the known path or its reported base
conda activate HIV_U3analysis                        # activate the shared env we're installing into

INSTALLER="conda"                                    # default package installer
command -v mamba >/dev/null 2>&1 && INSTALLER="mamba"  # prefer mamba if available (much faster solver)

echo "=== installing hifiasm + chopper + bc via ${INSTALLER} ==="  # progress marker
"${INSTALLER}" install -y -n HIV_U3analysis -c bioconda -c conda-forge hifiasm chopper bc  # install the long-read tools + bc into the env

echo "=== verifying ==="                             # progress marker
for t in hifiasm chopper bc; do                      # confirm each tool actually landed on PATH
    if command -v "$t" >/dev/null 2>&1; then          # is it resolvable now?
        echo "OK   $t -> $(command -v $t)"           # report success and its location
    else
        echo "WARN $t still not on PATH after install" >&2  # warn (don't fail) if it's still missing
    fi
done
echo "Done."                                         # final marker
