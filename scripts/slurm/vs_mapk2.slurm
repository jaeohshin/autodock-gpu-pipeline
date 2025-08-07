#!/bin/bash
#SBATCH --job-name=vs_mapk2
#SBATCH --output=logs/slurm_vs_mapk2.out
#SBATCH --error=logs/slurm_vs_mapk2.err
#SBATCH --partition=normal
#SBATCH --gres-flags=enforce-binding
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-task=1
#SBATCH --time=48:00:00

source ~/.bashrc
conda activate dock

echo "[INFO] Starting virtual screening for mapk2 on $(hostname)"
start_time=$(date +%s)

python run_vs.py \
  --project vs_crystal \
  --mode crystal \
  --kinase mapk2 \
  --nprocs 8

status=$?
end_time=$(date +%s)
duration=$((end_time - start_time))

if [ $status -eq 0 ]; then
    echo -e "mapk2\tcrystal\tSUCCESS\t${duration} sec" >> logs/vs_status.tsv
else
    echo -e "mapk2\tcrystal\tFAIL\t${duration} sec" >> logs/vs_status.tsv
fi
