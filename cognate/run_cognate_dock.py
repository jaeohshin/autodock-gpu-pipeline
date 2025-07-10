"""
2025/07/10
Cognate docking script using Meeko + docking_utils.py
"""

import os
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

# Make sure we can import docking_utils
from docking_utils import (
    generate_gpf,
    run_autogrid,
    run_docking,
    read_grid_center,
    get_atom_types_from_ligand,
    GRID_SIZE,
    prepare_receptor
)

KINASE_LIST = "kinase.txt"

def prepare_ligand_meeko(pdb_path, pdbqt_out):
    mol = Chem.MolFromPDBFile(pdb_path, removeHs=True)
    if mol is None:
        raise ValueError(f"Failed to load ligand: {pdb_path}")

    mol = Chem.AddHs(mol)

    if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
        raise RuntimeError(f"RDKit failed to embed 3D coordinates for: {pdb_path}")

    AllChem.UFFOptimizeMolecule(mol)

    prep = MoleculePreparation()
    mols_prepped = prep.prepare(mol)
    writer = PDBQTWriterLegacy()
    pdbqt_str = writer.write_string(mols_prepped[0])[0]

    with open(pdbqt_out, "w") as f:
        f.write(pdbqt_str)

    print(f"[INFO] Ligand PDBQT saved to {pdbqt_out}")

def fix_pdbqt_atom_names(pdbqt_file):
    fixed_lines = []
    with open(pdbqt_file, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                raw_name = line[12:16].strip()
                # Truncate long names (like HG A1 → GA1 or A1)
                cleaned = raw_name[-4:].rjust(4)
                fixed_line = line[:12] + cleaned + line[16:]
                fixed_lines.append(fixed_line)
            else:
                fixed_lines.append(line)
    with open(pdbqt_file, 'w') as f:
        f.writelines(fixed_lines)
    print(f"[INFO] Atom names cleaned in {pdbqt_file}")

def run_cognate_docking(kinase):
    print(f"\n=== Running {kinase} ===")

    base = os.path.join(kinase)
    rec_pdb = os.path.join(base, "receptor", "receptor.pdb")
    rec_pdbqt = os.path.join(base, "receptor", "receptor.pdbqt")
    lig_pdb = os.path.join(base, "ligands", "ligand.pdb")
    lig_pdbqt = os.path.join(base, "ligands", "ligand.pdbqt")
    grid_txt = os.path.join(base, "grid", "grid_center.txt")
    grid_dir = os.path.join(base, "grid")
    docking_dir = os.path.join(base, "docking")

    os.makedirs(docking_dir, exist_ok=True)

    # Step 1 – Prepare ligand
    prepare_ligand_meeko(lig_pdb, lig_pdbqt)

    # Step 2 – Prepare receptor.pdbqt
#    prepare_receptor(rec_pdb, rec_pdbqt)

    # Step 3 – Read grid center
    center = read_grid_center(grid_txt)

    # Step 4 – Generate GPF
    gpf_file = os.path.join(grid_dir, "grid.gpf")
    atom_types = get_atom_types_from_ligand(lig_pdbqt)
    generate_gpf(lig_pdbqt, rec_pdbqt, gpf_file, center, GRID_SIZE, atom_types)

    # Step 5 – Run AutoGrid
    run_autogrid(gpf_file)
    fld_file = os.path.join(grid_dir, "receptor.maps.fld")

    # Step 6 – Run docking
    output_prefix = os.path.join(docking_dir, "ligand_docked")
    run_docking(lig_pdbqt, fld_file, output_prefix)

    print(f"[DONE] {kinase}")

if __name__ == "__main__":
    with open(KINASE_LIST) as f:
        for line in f:
            if not line.strip():
                continue
            kinase = line.strip().split()[0]
            run_cognate_docking(kinase)
