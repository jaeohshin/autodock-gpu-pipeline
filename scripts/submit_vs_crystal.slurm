#!/bin/bash
#SBATCH --job-name=vs_array
#SBATCH --output=logs/slurm_vs_%A_%a.out
#SBATCH --error=logs/slurm_vs_%A_%a.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --array=0-25

source ~/.bashrc
conda activate dock

KINASES=(
  "abl1" "akt1" "akt2" "braf" "cdk2" "csf1r" "egfr" "fak1" "fgfr1"
  "igf1r" "jak2" "kpcb" "kit" "lck" "mapk2" "met" "mk01" "mk10"
  "mk14" "mp2k1" "plk1" "rock1" "src" "tgfr1" "vgfr2" "wee1"
)

KINASE=${KINASES[$SLURM_ARRAY_TASK_ID]}
LOGFILE="logs/vs_status.tsv"
MODE="crystal"
PROJECT="vs_crystal"

echo "[INFO] Starting $KINASE at $(date)"

START_TIME=$(date +%s)

python scripts/run_vs.py \
    --project "$PROJECT" \
    --mode "$MODE" \
    --kinase "$KINASE" \
    --nprocs 8

STATUS=$?
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

if [ $STATUS -eq 0 ]; then
    echo -e "$KINASE\t$MODE\tSUCCESS\t$DURATION sec" >> "$LOGFILE"
    echo "[✓] Finished $KINASE in $DURATION sec"
else
    echo -e "$KINASE\t$MODE\tFAIL\t$DURATION sec" >> "$LOGFILE"
    echo "[✗] $KINASE failed after $DURATION sec"
fi

