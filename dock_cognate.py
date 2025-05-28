"""
2025/05/20
Jaeoh Shin, Korea Institute for Advanced Study.
"""

import os
from glob import glob
from docking_utils import (
    prepare_receptor,
    generate_gpf,
    run_autogrid,
    run_docking,
    calc_ligand_center,
    run_cmd,
    PREPARE_LIGAND
)

ROOT_DIR = "./cognate"
GRID_SIZE = (60, 60, 60)

def run_cognate_docking(kinase_path):
    receptor_pdb = os.path.join(kinase_path, "receptor.pdb")
    ligand_pdb   = os.path.join(kinase_path, "ligand.pdb")
    output_dir   = os.path.join(kinase_path, "docking_output")
    os.makedirs(output_dir, exist_ok=True)

    receptor_pdbqt = os.path.join(output_dir, "receptor.pdbqt")
    ligand_pdbqt   = os.path.join(output_dir, "ligand.pdbqt")
    gpf_file       = os.path.join(output_dir, "grid.gpf")
    fld_file       = os.path.join(output_dir, "receptor.maps.fld")
    result_base    = os.path.join(output_dir, "ligand_docked")

    print(f"[INFO] Processing: {kinase_path}")
    prepare_receptor(receptor_pdb, receptor_pdbqt)

    # Always use MGLTools for ligand preparation (skip Meeko)
    print(f"[INFO] Preparing ligand using MGLTools only.")
    run_cmd(f"{PREPARE_LIGAND} -l {ligand_pdb} -o {ligand_pdbqt} -A hydrogens")

    center = calc_ligand_center(ligand_pdb)
    generate_gpf(ligand_pdbqt, receptor_pdbqt, gpf_file, center, GRID_SIZE)
    run_autogrid(gpf_file)
    run_docking(ligand_pdbqt, fld_file, result_base)

def main():
    kinase_dirs = sorted(glob(os.path.join(ROOT_DIR, "*")))
    for kinase_dir in kinase_dirs:
        if os.path.isdir(kinase_dir):
            run_cognate_docking(kinase_dir)

if __name__ == "__main__":
    main()
