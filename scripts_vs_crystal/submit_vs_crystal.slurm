#!/bin/sh
#SBATCH --job-name=vs_array
#SBATCH --output=vs_%A_%a.out
#SBATCH --error=vs_%A_%a.err
#SBATCH --partition=normal
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --time=120:00:00
#SBATCH --array=0-25%12  # 18 kinases total, max 9 concurrent jobs

KINASES=(abl1 akt1 akt2 braf cdk2 csf1r egfr fak1 fgfr1 igf1r jak2 kit kpcb lck mapk2 met mk01 mk10 mk14 mp2k1 plk1 rock1 src tgfr1 vgfr2 wee1)

source ~/miniconda3/etc/profile.d/conda.sh
conda activate dock

cd /store/jaeohshin/work/dock/scripts

KINASE=${KINASES[$SLURM_ARRAY_TASK_ID]}

echo "Running docking for kinase: $KINASE"

nvidia-smi
python --version
which python

python run_vs_crystal.py --kinase $KINASE
