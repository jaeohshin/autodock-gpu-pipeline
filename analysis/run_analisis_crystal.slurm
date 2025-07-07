#!/bin/bash
#SBATCH --job-name=debug_analyze
#SBATCH --output=logs/debug.out
#SBATCH --error=logs/debug.err
#SBATCH --time=2-23:00:00
#SBATCH --partition=normal
#SBATCH --nodelist=ada-004

# === Load conda and activate the dock environment ===
source ~/miniconda3/etc/profile.d/conda.sh
conda activate dock

# Define the list of kinase directories
KINASES=(akt2 mapk1 mapk10 mapkapk2 jak2 prkcb tgfbr1 abl1 akt1 met kit lck)

# Get kinase name for this array task
KINASE=${KINASES[$SLURM_ARRAY_TASK_ID]}

# Loop over all kinases
for KINASE in "${KINASES[@]}"; do
    echo "[$(date)] Running analysis for kinase: $KINASE"
    python analyze_crystal_docking.py --kinase "$KINASE"
done
