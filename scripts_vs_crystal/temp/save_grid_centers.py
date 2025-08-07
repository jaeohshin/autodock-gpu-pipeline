#!/usr/bin/env python3
"""
Compute grid centers from aligned ligand PDB files.

Expected input files:
    ligand_aligned_0001.pdb ~ ligand_aligned_0010.pdb

Output:
    grid_center_0001.txt ~ grid_center_0010.txt

Usage:
    cd /data/work/dock/cognate_bioemu/ABL1/ligands
    python ../../scripts/save_grid_centers.py
"""

import os
import numpy as np
from Bio.PDB import PDBParser


def calc_ligand_center(pdb_file):
    """Compute center of mass from PDB atoms."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("lig", pdb_file)
    coords = [atom.coord for atom in structure.get_atoms()]
    if not coords:
        raise ValueError(f"No atoms found in: {pdb_file}")
    return np.mean(coords, axis=0)


def main():
    aligned_dir = "./"
    num_structures = 50

    for i in range(1, num_structures + 1):
        aligned_pdb = os.path.join(aligned_dir, f"ligand_aligned_{i:04d}.pdb")
        center_txt = os.path.join(aligned_dir, f"grid_center_{i:04d}.txt")

        if os.path.exists(aligned_pdb):
            try:
                center = calc_ligand_center(aligned_pdb)
                with open(center_txt, "w") as f:
                    f.write(f"{center[0]:.3f} {center[1]:.3f} {center[2]:.3f}\n")
                print(f"[INFO] Saved: {center_txt}")
            except Exception as e:
                print(f"[ERROR] Failed to process {aligned_pdb}: {e}")
        else:
            print(f"[WARNING] File not found: {aligned_pdb}")


if __name__ == "__main__":
    main()
