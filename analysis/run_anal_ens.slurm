#!/bin/bash
#SBATCH --job-name=ens_analyze
#SBATCH --output=ens_logs/kinase_%A_%a.out
#SBATCH --error=ens_logs/kinase_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --partition=normal
#SBATCH --array=0-18

# === Load conda and activate the dock environment ===
source ~/miniconda3/etc/profile.d/conda.sh
conda activate dock

# Create logs directory if it doesn't exist
mkdir -p ens_logs

# Define successful kinases
KINASES=(abl1 braf fak1 fgfr1 igf1r jak2 kit mapk2 met mk01 wee1 mp2k1 csf1r mk10 rock1 tgfr1 plk1 akt1 kpcb)
# Get kinase name for this array task
KINASE=${KINASES[$SLURM_ARRAY_TASK_ID]}

echo "[$(date)] Running analysis for kinase: $KINASE"
#python anal_ens_dock.py --kinase "$KINASE" --max_receptor 50
python anal_mean.py --kinase "$KINASE" --max_receptor 50 --method mean

