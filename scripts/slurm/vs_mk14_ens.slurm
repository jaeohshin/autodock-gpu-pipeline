#!/bin/bash
#SBATCH --job-name=vs_mk14_ens
#SBATCH --output=logs/slurm_vs_mk14_ens.out
#SBATCH --error=logs/slurm_vs_mk14_ens.err
#SBATCH --partition=normal
#SBATCH --gres-flags=enforce-binding
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-task=1
#SBATCH --time=48:00:00

source ~/.bashrc
conda activate dock

echo "[INFO] Starting ensemble docking for mk14 on $(hostname)"
start_time=$(date +%s)

python run_vs.py \
  --project virtual_screening \
  --mode ensemble \
  --kinase mk14 \
  --nprocs 8

status=$?
end_time=$(date +%s)
duration=$((end_time - start_time))

if [ $status -eq 0 ]; then
    echo -e "mk14\tensemble\tSUCCESS\t${duration} sec" >> logs/vs_status.tsv
else
    echo -e "mk14\tensemble\tFAIL\t${duration} sec" >> logs/vs_status.tsv
fi
