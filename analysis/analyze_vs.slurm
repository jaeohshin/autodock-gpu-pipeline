#!/bin/bash
#SBATCH --job-name=debug_analyze
#SBATCH --output=logs/debug.out
#SBATCH --error=logs/debug.err
#SBATCH --time=2-23:00:00
#SBATCH --partition=normal
#SBATCH --nodelist=ada-004
#SBATCH --cpus-per-task=1

# Debug prints
echo "SLURM job started on $(hostname)"
echo "Running as user: $(whoami)"
echo "Current directory: $(pwd)"
echo "Loading conda env..."

# Most robust way to activate conda env
eval "$(conda shell.bash hook)"
conda activate dock

echo "Python path: $(which python)"
python --version

# Test simple Python run
echo "Testing simple python call..."
python -c "print('✓ Python is working!')"

# Now test your script
echo "Now running analyze script..."
python analyze_ensemble_docking.py --kinase braf --max_receptor 15

echo "Job completed"

