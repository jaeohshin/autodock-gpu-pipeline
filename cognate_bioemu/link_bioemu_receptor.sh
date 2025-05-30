#!/bin/bash

REFERENCE_LIST="kinase.list"
base_src="/data/work/flowpacker/samples"
base_dst="/data/work/dock/cognate_bioemu"

# Default number of structures
NUM_STRUCTURES=10

# Parse optional --num_structures argument
while [[ $# -gt 0 ]]; do
    case "$1" in
        --num_structures)
            NUM_STRUCTURES="$2"
            shift 2
            ;;
        *)
            echo "[ERROR] Unknown argument: $1"
            exit 1
            ;;
    esac
done

echo "[INFO] Linking up to $NUM_STRUCTURES structures per kinase."

while read -r kinase pdbid; do
    echo "[INFO] Processing $kinase..."

    kinase_lower=$(echo "$kinase" | tr '[:upper:]' '[:lower:]')
    src_dir="$base_src/$kinase_lower/run_1"
    dst_dir="$base_dst/$kinase/receptor"

    mkdir -p "$dst_dir"

    for i in $(seq 1 "$NUM_STRUCTURES"); do
        src_file="$src_dir/frame_re_${i}.pdb"
        dst_file="$dst_dir/receptor_$(printf "%04d" $i).pdb"

        if [[ -f "$src_file" ]]; then
            ln -sf "$src_file" "$dst_file"
            echo "  [✓] Linked: $dst_file"
        else
            echo "  [!] Missing: $src_file"
        fi
    done

done < "$REFERENCE_LIST"

