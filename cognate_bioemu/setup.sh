#!/bin/bash

REFERENCE_LIST="kinase.list"

# Residues to exclude from ligand extraction
EXCLUDE_RESNS="HOH|NA|CL|SO4|PO4|MG|CA|K|ZN|FE|MN|GOL|DMS|ACT|EDO|PEG|TRS|MSE"

while read -r kinase pdbid; do
    echo "[INFO] Processing $kinase (PDB: $pdbid)..."

    mkdir -p "$kinase/receptor"
    mkdir -p "$kinase/ligands"

    ### 1. Skip receptor copy — assume you manually put receptor files in ./receptor/
    echo "  [OK] Assuming receptor PDBs are already placed manually in $kinase/receptor/"

    ### 2. Download original PDB from RCSB for ligand extraction
    pdb_file="$kinase/receptor/${pdbid}.pdb"
    if [[ ! -f "$pdb_file" ]]; then
        wget -q "https://files.rcsb.org/download/${pdbid}.pdb" -O "$pdb_file"
        echo "  [OK] Downloaded $pdbid.pdb from RCSB"
    fi

    # Link to original for reference
    ln -sf "$pdbid.pdb" "$kinase/receptor/original.pdb"

    ### 3. Extract largest ligand from chain A
    ligand_out="$kinase/ligands/ligand.pdb"

    awk -v exclude="$EXCLUDE_RESNS" '
    $1 == "HETATM" {
        resn = substr($0, 18, 3)
        chain = substr($0, 22, 1)
        resi = substr($0, 23, 4)
        key = resn "_" resi "_" chain
        if (chain == "A" && match(resn, exclude) == 0) {
            count[key]++
            line[key] = line[key] $0 "\n"
        }
    }
    END {
        max = 0
        for (k in count) {
            if (count[k] > max) {
                max = count[k]
                best = k
            }
        }
        if (best in line) {
            printf "%s", line[best]
        }
    }' "$pdb_file" > "$ligand_out"

    if [[ ! -s "$ligand_out" ]]; then
        echo "  [WARNING] No ligand found in chain A of $pdb_file"
        rm -f "$ligand_out"
    else
        echo "  [OK] Ligand extracted to $ligand_out"
    fi

done < "$REFERENCE_LIST"
