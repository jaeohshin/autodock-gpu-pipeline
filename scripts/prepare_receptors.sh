#!/bin/bash
# Make structures using BioEmu etc. then copy them to dock directory

REFERENCE_LIST="/store/jaeohshin/work/dock/kinase.txt"
SRC_ROOT="/store/jaeohshin/work/kine/output"
DST_ROOT="/store/jaeohshin/work/dock/virtual_screening/input/receptors"

while read -r kinase _; do
    kinase_lower="${kinase,,}"
    dst_dir="${DST_ROOT}/${kinase_lower}"
    receptor_list="${dst_dir}/receptorlist.txt"

    exec > >(tee "${dst_dir}/${kinase_lower}_copy.log") 2>&1

    echo "[INFO] Processing kinase: $kinase_lower"
    mkdir -p "$dst_dir"
    > "$receptor_list"

    for i in $(seq 1 50); do
	padded_i=$(printf "%03d" $i)
        src="${SRC_ROOT}/${kinase_lower}/final_str/receptor_${padded_i}.pdb"
        dst_filename="receptor_${padded_i}.pdb"
        dst="${dst_dir}/${dst_filename}"

        if [ -f "$src" ]; then
            cp -u "$src" "$dst"
            echo "$dst_filename" >> "$receptor_list"
        else
            echo "[WARN] Missing: $src"
        fi
    done
done < "$REFERENCE_LIST"

