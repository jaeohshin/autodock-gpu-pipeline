#!/bin/bash

REFERENCE_LIST="/data/work/dock/virtual_screening/input/abl1.txt"
SRC_ROOT="/data/work/flowpacker/samples"
DST_ROOT="/data/work/dock/virtual_screening/input/receptors"

while read -r kinase _; do
    kinase_lower="${kinase,,}"
    dst_dir="${DST_ROOT}/${kinase_lower}"
    receptor_list="${dst_dir}/receptorlist.txt"

    echo "[INFO] Processing kinase: $kinase_lower"
    mkdir -p "$dst_dir"
    > "$receptor_list"

    for i in $(seq 1 20); do
        src="${SRC_ROOT}/${kinase_lower}/new_relaxation/frame_re_${i}.pdb"
        dst_filename="receptor_$(printf "%04d" $i).pdb"
        dst="${dst_dir}/${dst_filename}"

        if [ -f "$src" ]; then
            cp "$src" "$dst"
            echo "$dst_filename" >> "$receptor_list"
        else
            echo "[WARN] Missing: $src"
        fi
    done
done < "$REFERENCE_LIST"

