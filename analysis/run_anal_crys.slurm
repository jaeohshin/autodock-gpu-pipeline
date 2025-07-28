#!/bin/bash
#SBATCH --job-name=debug_analyze
#SBATCH --output=logs/debug.out
#SBATCH --error=logs/debug.err
#SBATCH --time=24:00:00
#SBATCH --partition=normal

# === Load conda and activate the dock environment ===
source ~/miniconda3/etc/profile.d/conda.sh
conda activate dock

# Define the list of kinase directories
KINASES=(abl1 akt1 akt2 braf cdk2 csf1r egfr fak1 fgfr1 igf1r jak2 kpcb kit lck mapk2 met mk01 mk10 mk14 mp2k1 plk1 rock1 src tgfr1 vgfr2 wee1)
#KINASES=(egfr fak1 fgfr1 igf1r jak2 kpcb kit lck mapk2 met mk01 mk10 mk14 mp2k1 plk1 rock1 src tgfr1 vgfr2 wee1)
# Get kinase name for this array task
KINASE=${KINASES[$SLURM_ARRAY_TASK_ID]}

# Loop over all kinases
for KINASE in "${KINASES[@]}"; do
    echo "[$(date)] Running analysis for kinase: $KINASE"
    python analyze_crystal_docking.py --kinase "$KINASE"
done
