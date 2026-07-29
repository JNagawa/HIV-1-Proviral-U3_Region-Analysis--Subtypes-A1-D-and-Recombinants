#!/bin/bash
CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"   # expected conda init script location
if [ -f "$CONDA_SH" ]; then                          # prefer the known miniconda path if present
    source "$CONDA_SH"                               # load conda so `conda activate` works here
elif command -v conda >/dev/null 2>&1; then          # else, if conda is on PATH...
    source "$(conda info --base)/etc/profile.d/conda.sh"  # ...load it from its reported base
else
    echo "ERROR: Conda not found." >&2               # no conda at all
    exit 1                                           # nothing to check against, so abort
fi
# activate the env whose tools we want to verify
conda activate HIV_U3analysis

# the tools the pipeline expects to be installed
tools=("NanoFilt" "NanoStat" "porechop" "porechop_abi" "nanoQC" "prefetch" "fasterq-dump" "trimmomatic" "cutadapt" "kraken2" "ccs" "NanoPlot" "shiver_init.sh" "shiver_align.sh")

echo "Checking tools in HIV_U3analysis environment:"  # header for the report
for tool in "${tools[@]}"; do                        # check each expected tool in turn
    if command -v "$tool" >/dev/null 2>&1; then       # is it resolvable on PATH?
        echo "FOUND: $tool ($(command -v "$tool"))"  # report it found, with its path
    else
        echo "MISSING: $tool"                        # report it missing
    fi
done
