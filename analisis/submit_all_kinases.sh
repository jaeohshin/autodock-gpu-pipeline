#!/bin/bash
#SBATCH --job-name=process_all_kinases
#SBATCH --output=process_kinases_%j.out
#SBATCH --error=process_kinases_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16

# Job information
echo "=========================================="
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Running on host: $(hostname)"
echo "Job started at: $(date)"
echo "Working directory: $(pwd)"
echo "=========================================="

# Load required modules (adjust based on your cluster)
# module load python/3.8
# or if using conda:
# module load anaconda3
# conda activate your_env

# Ensure required packages are installed
echo "Checking Python packages..."
pip install --user pandas numpy tqdm

# Set the base directory
BASE_DIR="/store/jaeohshin/work/dock/virtual_screening/docking_output"

# Create output directory for results
OUTPUT_DIR="docking_results_${SLURM_JOB_ID}"
mkdir -p $OUTPUT_DIR

# Run the processing
echo "Starting DLG processing for all kinases..."
echo "Base directory: $BASE_DIR"
echo "Using $SLURM_CPUS_PER_TASK processes"

python ada-cluster-processor.py \
    --base-dir $BASE_DIR \
    --all \
    --processes $SLURM_CPUS_PER_TASK

# Check if processing was successful
if [ $? -eq 0 ]; then
    echo "Processing completed successfully!"
    
    # Move results to output directory
    mv docking_results_clean_*.csv $OUTPUT_DIR/
    mv docking_results_summary_*.json $OUTPUT_DIR/
    
    # Create a summary file
    echo "Creating job summary..."
    cat > $OUTPUT_DIR/job_summary.txt <<EOF
SLURM Job Summary
=================
Job ID: $SLURM_JOB_ID
Start time: $(date)
Base directory: $BASE_DIR
Processes used: $SLURM_CPUS_PER_TASK
Output directory: $OUTPUT_DIR

Files created:
$(ls -lh $OUTPUT_DIR/)
EOF
    
    # Compress results for easy download
    echo "Compressing results..."
    tar -czf docking_results_${SLURM_JOB_ID}.tar.gz $OUTPUT_DIR/
    
    echo "=========================================="
    echo "RESULTS READY FOR DOWNLOAD:"
    echo "docking_results_${SLURM_JOB_ID}.tar.gz"
    echo "=========================================="
    
else
    echo "ERROR: Processing failed!"
    exit 1
fi

echo "Job completed at: $(date)"
echo "=========================================="
