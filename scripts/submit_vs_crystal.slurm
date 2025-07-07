#!/bin/sh
#SBATCH --job-name=vs_array
#SBATCH --output=vs_%A_%a.out
#SBATCH --error=vs_%A_%a.err
#SBATCH --partition=normal
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --time=120:00:00
#SBATCH --array=0-17%9  # 18 kinases total, max 9 concurrent jobs

KINASES=(abl1 akt1 akt2 braf egfr fgfr1 igf1r jak2 kit lck mapk1 mapk10 mapkapk2 met prkcb tgfbr1 wee1)

source ~/miniconda3/etc/profile.d/conda.sh
conda activate dock

cd /store/jaeohshin/work/dock/scripts

KINASE=${KINASES[$SLURM_ARRAY_TASK_ID]}

echo "Running docking for kinase: $KINASE"

nvidia-smi
python --version
which python

python run_vs_crystal.py --kinase $KINASE
