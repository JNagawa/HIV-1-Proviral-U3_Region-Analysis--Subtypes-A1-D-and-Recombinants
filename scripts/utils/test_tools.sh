#!/bin/bash
CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"
if [ -f "$CONDA_SH" ]; then
    source "$CONDA_SH"
elif command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
else
    echo "ERROR: Conda not found." >&2
    exit 1
fi
conda activate HIV_U3analysis

tools=("NanoFilt" "NanoStat" "porechop" "porechop_abi" "nanoQC" "prefetch" "fasterq-dump" "trimmomatic" "cutadapt" "kraken2" "ccs" "NanoPlot" "shiver_init.sh" "shiver_align.sh")

echo "Checking tools in HIV_U3analysis environment:"
for tool in "${tools[@]}"; do
    if command -v "$tool" >/dev/null 2>&1; then
        echo "FOUND: $tool ($(command -v "$tool"))"
    else
        echo "MISSING: $tool"
    fi
done
