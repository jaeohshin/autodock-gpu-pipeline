#!/bin/bash

KINASES=(
  "abl1" "akt1" "akt2" "braf" "cdk2" "csf1r" "egfr" "fak1" "fgfr1"
  "igf1r" "jak2" "kpcb" "kit" "lck" "mapk2" "met" "mk01" "mk10"
  "mk14" "mp2k1" "plk1" "rock1" "src" "tgfr1" "vgfr2" "wee1"
)

for kinase in "${KINASES[@]}"; do
cat <<EOF > vs_${kinase}_ens.slurm
#!/bin/bash
#SBATCH --job-name=vs_${kinase}_ens
#SBATCH --output=logs/slurm_vs_${kinase}_ens.out
#SBATCH --error=logs/slurm_vs_${kinase}_ens.err
#SBATCH --partition=normal
#SBATCH --gres-flags=enforce-binding
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-task=1
#SBATCH --time=48:00:00

source ~/.bashrc
conda activate dock

echo "[INFO] Starting ensemble docking for ${kinase} on \$(hostname)"
start_time=\$(date +%s)

python run_vs.py \\
  --project virtual_screening \\
  --mode ensemble \\
  --kinase ${kinase} \\
  --nprocs 8

status=\$?
end_time=\$(date +%s)
duration=\$((end_time - start_time))

if [ \$status -eq 0 ]; then
    echo -e "${kinase}\tensemble\tSUCCESS\t\${duration} sec" >> logs/vs_status.tsv
else
    echo -e "${kinase}\tensemble\tFAIL\t\${duration} sec" >> logs/vs_status.tsv
fi
EOF
done

echo "[✓] All ensemble SLURM scripts generated: vs_<kinase>_ens.slurm"

