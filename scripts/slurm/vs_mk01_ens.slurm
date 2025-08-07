#!/bin/bash
#SBATCH --job-name=vs_mk01_ens
#SBATCH --output=logs/slurm_vs_mk01_ens.out
#SBATCH --error=logs/slurm_vs_mk01_ens.err
#SBATCH --partition=normal
#SBATCH --gres-flags=enforce-binding
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-task=1
#SBATCH --time=48:00:00

source ~/.bashrc
conda activate dock

echo "[INFO] Starting ensemble docking for mk01 on $(hostname)"
start_time=$(date +%s)

python run_vs.py \
  --project virtual_screening \
  --mode ensemble \
  --kinase mk01 \
  --nprocs 8

status=$?
end_time=$(date +%s)
duration=$((end_time - start_time))

if [ $status -eq 0 ]; then
    echo -e "mk01\tensemble\tSUCCESS\t${duration} sec" >> logs/vs_status.tsv
else
    echo -e "mk01\tensemble\tFAIL\t${duration} sec" >> logs/vs_status.tsv
fi
