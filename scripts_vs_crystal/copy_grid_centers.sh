#!/bin/bash
set -e

COGNATE_DIR="../cognate"
TARGET_BASE="../virtual_screening/preprocessed/grid_crystal_centers"

for kinase_dir in "$COGNATE_DIR"/*/ ; do
    kinase=$(basename "$kinase_dir")
    src="$kinase_dir/grid_center.txt"
    tgt_dir="$TARGET_BASE/${kinase,,}"   # lowercase kinase name
    tgt="$tgt_dir/receptor_crystal.txt"
    
    if [ -f "$src" ]; then
        mkdir -p "$tgt_dir"
        cp "$src" "$tgt"
        echo "Copied $src to $tgt"
    else
        echo "Warning: $src not found, skipping."
    fi
done

