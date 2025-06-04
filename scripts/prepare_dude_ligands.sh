#!/bin/bash

ROOT="/data/work/dock/dude_raw"
OUT_ROOT="/data/work/dock/virtual_screening/input/ligands"

echo "[INFO] Starting SDF-based .pdbqt conversion with CHEMBL ID naming..."

for target_path in "$ROOT"/*/; do
    target=$(basename "${target_path%/}")
    echo "[TARGET] $target"

    for type in actives decoys; do
        sdf_gz="${target_path%/}/${type}_final.sdf.gz"
        sdf="${sdf_gz%.gz}"
        out_dir="${OUT_ROOT%/}/$target/$type"
        list_out="${OUT_ROOT%/}/$target/ligands_${type}.list"
        mkdir -p "$out_dir"

        if [[ -f "$sdf_gz" ]]; then
            echo "  → Decompressing $sdf_gz"
            gunzip -kf "$sdf_gz"
        fi

        if [[ -f "$sdf" ]]; then
            echo "  → Converting $sdf to temporary .pdbqt files"
            obabel "$sdf" -O "$out_dir/tmp.pdbqt" -m --partialcharge gasteiger

            echo "  → Renaming files and writing list $list_out"
            cd "$out_dir" || { echo "[ERROR] Cannot enter $out_dir"; exit 1; }
            > "$list_out"
            i=1
            for f in tmp*.pdbqt; do
                chembl_id=$(head -n 1 "$f" | tr -d ' \t\r\n')
                if [ -z "$chembl_id" ]; then
                    chembl_id="unknown"
                fi
                newname=$(printf "%s_%05d_%s.pdbqt" "$type" "$i" "$chembl_id")
                mv "$f" "$newname"
                echo "$newname" >> "$list_out"
                ((i++))
            done
        else
            echo "  ⚠️ Missing SDF file: $sdf"
        fi
    done
done

echo "[DONE] Conversion and renaming finished."
