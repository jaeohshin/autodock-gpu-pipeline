#!/bin/bash
#SBATCH --job-name=process_all_kinases
#SBATCH --output=process_kinases_%j.out
#SBATCH --error=process_kinases_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=normal

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "=============================="

source ~/.bashrc
conda activate dock

python -u dlg_processor.py \
    --base-dir /store/jaeohshin/work/dock/virtual_screening/docking_output \
    --all \
    --processes $SLURM_CPUS_PER_TASK

echo "=============================="
echo "Job finished: $(date)"
