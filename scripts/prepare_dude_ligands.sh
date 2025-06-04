#!/bin/bash

ROOT="/data/work/dock/dude_raw"
OUT_ROOT="/data/work/dock/virtual_screening/input/ligands"


echo "[INFO] Starting SDF-based .pdbqt conversion..."

for target_path in "$ROOT"/*/; do
    target=$(basename "$target_path")
    echo "[TARGET] $target"

    for type in actives decoys; do
        sdf_gz="$target_path/${type}_final.sdf.gz"
        sdf="${sdf_gz%.gz}"
        out_dir="$OUT_ROOT/$target/$type"
        mkdir -p "$out_dir"

        # Step 1: decompress if needed
        if [[ -f "$sdf_gz" ]]; then
            echo "  → Decompressing $sdf_gz"
            gunzip -kf "$sdf_gz"
        fi

        # Step 2: convert to .pdbqt
        if [[ -f "$sdf" ]]; then
            echo "  → Converting $sdf to .pdbqt"
            obabel "$sdf" -O "$out_dir/ligand.pdbqt" -m --partialcharge gasteiger
        else
            echo "  ⚠️  Missing SDF file: $sdf"
        fi
    done
done

echo "[DONE] All ligands converted using SDF input."

