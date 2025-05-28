#!/usr/bin/env python3
"""
Fix ligand PDB extracted from an experimental structure:
- Removes CONECT records
- Converts to SMILES via Open Babel
- Uses RDKit to rebuild 3D structure with hydrogens
- Outputs a clean, dockable PDB file
"""

import sys
import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

def run_cmd_get_output(cmd):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] Command failed: {cmd}")
        print(result.stderr)
        sys.exit(1)
    return result.stdout.strip()

def remove_conect(input_file, temp_file):
    with open(input_file, "r") as fin, open(temp_file, "w") as fout:
        for line in fin:
            if not line.startswith("CONECT"):
                fout.write(line)
    print(f"[INFO] CONECT records removed: {temp_file}")

def generate_smiles_from_pdb(pdb_file):
    cmd = f"obabel {pdb_file} -osmi -xh"
    print(f"[RUN] {cmd}")
    smiles = run_cmd_get_output(cmd)
    print(f"[INFO] SMILES: {smiles}")
    return smiles

def smiles_to_3d_pdb(smiles, output_pdb):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("[ERROR] RDKit failed to parse SMILES")
        sys.exit(1)

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)
    Chem.MolToPDBFile(mol, output_pdb)

    total = mol.GetNumAtoms()
    heavy = mol.GetNumHeavyAtoms()
    print(f"[INFO] Heavy atoms: {heavy}, Hydrogens: {total - heavy}")
    print(f"[INFO] Final PDB written to: {output_pdb}")

def main(input_pdb, output_pdb):
    temp_clean = "ligand_no_conect.pdb"
    remove_conect(input_pdb, temp_clean)
    smiles = generate_smiles_from_pdb(temp_clean)
    smiles_to_3d_pdb(smiles, output_pdb)
    os.remove(temp_clean)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fix_ligand_pdb_smiles.py ligand.pdb ligand_fixed.pdb")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
