#!/bin/bash
#SBATCH --job-name=vs_braf
#SBATCH --output=logs/slurm_vs_braf.out
#SBATCH --error=logs/slurm_vs_braf.err
#SBATCH --partition=normal
#SBATCH --gres-flags=enforce-binding
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-task=1
#SBATCH --time=48:00:00

source ~/.bashrc
conda activate dock

echo "[INFO] Starting virtual screening for braf on $(hostname)"
start_time=$(date +%s)

python run_vs.py \
  --project vs_crystal \
  --mode crystal \
  --kinase braf \
  --nprocs 8

status=$?
end_time=$(date +%s)
duration=$((end_time - start_time))

if [ $status -eq 0 ]; then
    echo -e "braf\tcrystal\tSUCCESS\t${duration} sec" >> logs/vs_status.tsv
else
    echo -e "braf\tcrystal\tFAIL\t${duration} sec" >> logs/vs_status.tsv
fi
