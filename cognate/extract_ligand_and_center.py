import os
from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np
from collections import defaultdict

KINASE_LIST = "kinase.list"

EXCLUDED = {"HOH", "WAT", "SO4", "CL", "NA", "K", "MG", "GOL", "CA", "ZN", "PO4", "PEG"}

class SpecificLigandSelect(Select):
    def __init__(self, target_residue):
        self.target_residue = target_residue

    def accept_residue(self, residue):
        return residue is self.target_residue

def find_largest_ligand_in_chainA(structure):
    ligands = defaultdict(list)
    for model in structure:
        for chain in model:
            if chain.id != "A":
                continue
            for residue in chain:
                if residue.id[0] != " " and residue.resname not in EXCLUDED:
                    ligands[residue].extend([atom.coord for atom in residue])
    if not ligands:
        return None, None
    # Select residue with the most atoms
    best_residue = max(ligands.items(), key=lambda x: len(x[1]))[0]
    coords = np.array(ligands[best_residue])
    com = coords.mean(axis=0)
    return best_residue, com

def process_kinase(kinase):
    folder = kinase
    original_pdb = os.path.join(folder, "original.pdb")
    ligand_pdb = os.path.join(folder, "ligand.pdb")
    center_txt = os.path.join(folder, "grid_center.txt")

    if not os.path.isfile(original_pdb):
        print(f"[WARNING] {original_pdb} not found. Skipping.")
        return

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", original_pdb)

    residue, center = find_largest_ligand_in_chainA(structure)
    if residue is None:
        print(f"[WARNING] No valid ligand found in chain A for {kinase}")
        return

    # Save ligand PDB
    io = PDBIO()
    io.set_structure(structure)
    io.save(ligand_pdb, SpecificLigandSelect(residue))

    # Save grid center
    with open(center_txt, "w") as f:
        f.write(f"{center[0]:.3f} {center[1]:.3f} {center[2]:.3f}\n")

    print(f"[INFO] {kinase}: saved ligand.pdb ({residue.resname}) and grid_center.txt")

if __name__ == "__main__":
    with open(KINASE_LIST) as f:
        for line in f:
            kinase = line.strip().split()[0]
            process_kinase(kinase)
