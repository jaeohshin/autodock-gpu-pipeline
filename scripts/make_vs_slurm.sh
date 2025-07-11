#!/bin/bash

KINASES=(
  "abl1" "akt1" "akt2" "braf" "cdk2" "csf1r" "egfr" "fak1" "fgfr1"
  "igf1r" "jak2" "kpcb" "kit" "lck" "mapk2" "met" "mk01" "mk10"
  "mk14" "mp2k1" "plk1" "rock1" "src" "tgfr1" "vgfr2" "wee1"
)

for kinase in "${KINASES[@]}"; do
cat <<EOF > vs_${kinase}.slurm
#!/bin/bash
#SBATCH --job-name=vs_${kinase}
#SBATCH --output=logs/slurm_vs_${kinase}.out
#SBATCH --error=logs/slurm_vs_${kinase}.err
#SBATCH --partition=normal
#SBATCH --gres-flags=enforce-binding
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-task=1
#SBATCH --time=48:00:00

source ~/.bashrc
conda activate dock

echo "[INFO] Starting virtual screening for ${kinase} on \$(hostname)"
start_time=\$(date +%s)

python run_vs.py \\
  --project vs_crystal \\
  --mode crystal \\
  --kinase ${kinase} \\
  --nprocs 8

status=\$?
end_time=\$(date +%s)
duration=\$((end_time - start_time))

if [ \$status -eq 0 ]; then
    echo -e "${kinase}\tcrystal\tSUCCESS\t\${duration} sec" >> logs/vs_status.tsv
else
    echo -e "${kinase}\tcrystal\tFAIL\t\${duration} sec" >> logs/vs_status.tsv
fi
EOF
done

echo "[✓] All SLURM scripts generated: vs_<kinase>.slurm"

