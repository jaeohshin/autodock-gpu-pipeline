#!/bin/bash
#SBATCH --job-name=dlg_process
#SBATCH --output=dlg_process_%j.out
#SBATCH --error=dlg_process_%j.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=normal

# Job info
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"

# Load required modules (adjust based on Ada's modules)
module load python/3.8
# or if using conda:
# module load anaconda3

# Activate conda environment if using one
# conda activate base

# Install required packages if not already installed
pip install --user pandas numpy

# Create results directory
mkdir -p results

# Run the processing
echo "Starting DLG processing..."
python ada_dlg_processor.py --processes 16

# Create a package for easy download
echo "Creating download package..."
mkdir -p download_package
cp docking_results_clean_*.csv download_package/
cp docking_results_summary_*.json download_package/

# Create README
cat > download_package/README.txt << EOF
DUD-E Kinase Docking Results
Processed on: $(date)
Job ID: $SLURM_JOB_ID

Files:
- docking_results_clean_*.csv: Main results file for dashboard
- docking_results_summary_*.json: Summary statistics

Total kinases processed: 26
Total structures: ~1274 (49 receptors × 26 kinases)

To use:
1. Download this entire folder to your local computer
2. Run the dashboard with the CSV file
EOF

# Compress for download
tar -czf docking_results_${SLURM_JOB_ID}.tar.gz download_package/

echo "Job completed at: $(date)"
echo "Download: docking_results_${SLURM_JOB_ID}.tar.gz"