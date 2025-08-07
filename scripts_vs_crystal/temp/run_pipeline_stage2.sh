#!/bin/bash
#SBATCH --job-name=nvt_$1
#SBATCH --output=logs/nvt_%x.out
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00

set -e
source ~/miniconda3/etc/profile.d/conda.sh
conda activate gro

bash scripts/gromacs_scripts.sh "$1"
