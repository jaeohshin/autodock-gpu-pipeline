#!/bin/bash

# Array of all kinases
kinases=("abl1" "akt1" "akt2" "braf" "csf1r" "fak1" "fgfr1" "igf1r" "jak2" "kit" "kpcb" "mapk2" "met" "mk01" "mk10" "mp2k1" "plk1" "rock1" "tgfr1")

echo "Running anal_mean.py for ${#kinases[@]} kinases..."
echo "=========================================="

# Loop through each kinase
for kinase in "${kinases[@]}"; do
    echo "Processing $kinase..."
    python anal_mean.py --kinase "$kinase" --tmscore_file ranked_report.txt --min_tmscore 0.8 --method mean
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "✓ $kinase completed successfully"
    else
        echo "✗ $kinase failed"
    fi
    echo "----------"
done

echo "All kinases processed!"
