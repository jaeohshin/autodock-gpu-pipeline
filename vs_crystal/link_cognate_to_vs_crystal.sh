#!/bin/bash

KINASES=(
  "abl1" "akt1" "akt2" "braf" "cdk2" "csf1r" "egfr" "fak1" "fgfr1"
  "igf1r" "jak2" "kpcb" "kit" "lck" "mapk2" "met" "mk01" "mk10"
  "mk14" "mp2k1" "plk1" "rock1" "src" "tgfr1" "vgfr2" "wee1"
)

COGNATE_BASE="/store/jaeohshin/work/dock/cognate"
VS_BASE="/store/jaeohshin/work/dock/vs_crystal"
VS_PIPE_BASE="/store/jaeohshin/work/dock/virtual_screening"

echo "[INFO] Linking from COGNATE → VS_CRYSTAL..."

for kinase in "${KINASES[@]}"; do
    echo "→ [$kinase]"

    ## 1. receptor.pdb (raw PDB)
    src_pdb="${COGNATE_BASE}/${kinase}/receptor/receptor.pdb"
    dst_pdb_dir="${VS_BASE}/input/receptors/${kinase}"
    mkdir -p "$dst_pdb_dir"
    if [ -f "$src_pdb" ]; then
        ln -sfn "$src_pdb" "${dst_pdb_dir}/receptor.pdb"
        echo "   [✓] Linked receptor.pdb"
    else
        echo "   [⚠] Missing receptor.pdb at $src_pdb"
    fi

    ## 2. receptor.pdbqt
    src_pdbqt="${COGNATE_BASE}/${kinase}/grid/receptor.pdbqt"
    dst_pdbqt_dir="${VS_BASE}/preprocessed/receptors_pdbqt/${kinase}"
    mkdir -p "$dst_pdbqt_dir"
    if [ -f "$src_pdbqt" ]; then
        ln -sfn "$src_pdbqt" "${dst_pdbqt_dir}/receptor.pdbqt"
        echo "   [✓] Linked receptor.pdbqt"
    else
        echo "   [⚠] Missing receptor.pdbqt at $src_pdbqt"
    fi

    ## 3. grid_center.txt
    src_center="${COGNATE_BASE}/${kinase}/grid/grid_center.txt"
    dst_center_dir="${VS_BASE}/preprocessed/grid_centers/${kinase}"
    mkdir -p "$dst_center_dir"
    if [ -f "$src_center" ]; then
        ln -sfn "$src_center" "${dst_center_dir}/grid_center.txt"
        ln -sfn "$src_center" "${dst_center_dir}/receptor_crystal.txt"
        echo "   [✓] Linked grid_center.txt"
    else
        echo "   [⚠] Missing grid_center.txt at $src_center"
    fi

    ## 4. grid files (.map, .fld, and alias receptor_crystal.maps.fld)
    src_grid_dir="${COGNATE_BASE}/${kinase}/grid"
    dst_grid_dir="${VS_BASE}/grids/${kinase}"
    mkdir -p "$dst_grid_dir"
    map_found=false
    for mapfile in "$src_grid_dir"/*.map; do
        if [ -f "$mapfile" ]; then
            ln -sfn "$mapfile" "$dst_grid_dir/"
            map_found=true
        fi
    done

    fld_file="$src_grid_dir/receptor.maps.fld"
    if [ -f "$fld_file" ]; then
        ln -sfn "$fld_file" "$dst_grid_dir/receptor.maps.fld"
        ln -sfn "$fld_file" "$dst_grid_dir/receptor_crystal.maps.fld"
        echo "   [✓] Linked receptor.maps.fld and receptor_crystal.maps.fld"
    else
        echo "   [⚠] receptor.maps.fld missing in $src_grid_dir"
    fi

    if [ "$map_found" = true ]; then
        echo "   [✓] Linked .map files"
    else
        echo "   [⚠] No .map files found in $src_grid_dir"
    fi

done

## 5. Link ligands_sdf2meeko as ligands_pdbqt
src_ligands="${VS_PIPE_BASE}/input/ligands_sdf2meeko"
dst_ligands="${VS_BASE}/preprocessed/ligands_pdbqt"
ln -sfn "$src_ligands" "$dst_ligands"
echo "[✓] Linked ligands_pdbqt → $dst_ligands"


## 6. Link ligand list files (.list)
    for kinase in "${KINASES[@]}"; do
        src_list_dir="${VS_PIPE_BASE}/input/ligands_sdf2meeko/${kinase}"
        dst_list_dir="${VS_BASE}/preprocessed/ligands_pdbqt/${kinase}"
        mkdir -p "$dst_list_dir"

        for list_type in ligands_actives.list ligands_decoys.list; do
            src_list="${src_list_dir}/${list_type}"
            dst_list="${dst_list_dir}/${list_type}"
            if [ -f "$src_list" ]; then
                ln -sfn "$src_list" "$dst_list"
                echo "→ [✓] Linked $list_type for $kinase"
            else
                echo "→ [⚠] Missing $list_type for $kinase"
            fi
        done
    done


echo "[✅ DONE] All links created successfully."

