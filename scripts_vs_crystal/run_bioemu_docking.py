"""
2025/05/29
Jaeoh Shin, Korea Institute for Advanced Study.
"""
# Run: python run_cognate_docking.py

import os
from docking_utils import (
    GRID_SIZE,
    prepare_receptor, prepare_ligand,
    calc_ligand_center, read_grid_center,
    generate_gpf, run_autogrid, run_docking
)

top20_mode = True
ROOT_DIR = "../cognate_bioemu"  # Change to your actual target folder

def run_ensemble_docking(kinase_dir):
    print(f"\n[INFO] Processing ensemble docking: {kinase_dir}")

    receptor_dir = os.path.join(kinase_dir, "receptor")
    ligand_dir = os.path.join(kinase_dir, "ligands")
    output_base = os.path.join(kinase_dir, "output")
    if top20_mode:
        top20_path = os.path.join(receptor_dir, "top20.txt")
        with open(top20_path, "r") as f:
            indices = [
                os.path.splitext(line.strip())[0].split("_")[1]
                for line in f if line.strip().startswith("receptor_")
            ]
    else:
        indices = [f"{i:04d}" for i in range(1, 101)]

    for idx in indices:
        
        

        receptor_pdb = os.path.join(receptor_dir, f"receptor_{idx}.pdb")
        ligand_pdb   = os.path.join(ligand_dir, f"ligand_aligned_{idx}.pdb")
        center_file  = os.path.join(ligand_dir, f"grid_center_{idx}.txt")
        
        suffix = f"top20_{idx}" if top20_mode else f"ensemble_{idx}"
        output_dir   = os.path.join(output_base, suffix)
        
        os.makedirs(output_dir, exist_ok=True)

        if not os.path.exists(receptor_pdb):
            print(f"[WARN] Missing receptor file: {receptor_pdb}")
            continue
        if not os.path.exists(ligand_pdb):
            print(f"[WARN] Missing ligand file: {ligand_pdb}")
            continue
        if not os.path.exists(center_file):
            print(f"[WARN] Missing center file: {center_file}")
            continue

        grid_center = read_grid_center(center_file)

        receptor_pdbqt = os.path.join(output_dir, "receptor.pdbqt")
        ligand_pdbqt   = os.path.join(output_dir, "ligand.pdbqt")

        print(f"[INFO] Running docking for pair {idx}")
        prepare_receptor(receptor_pdb, receptor_pdbqt)
        prepare_ligand(ligand_pdb, ligand_pdbqt)

        gpf_file = os.path.join(output_dir, "grid.gpf")
        generate_gpf(ligand_pdbqt, receptor_pdbqt, gpf_file, grid_center, GRID_SIZE)

        run_autogrid(gpf_file)
        fld_file = os.path.join(output_dir, "receptor.maps.fld")

        output_basename = os.path.join(output_dir, "ligand_docked")
        run_docking(ligand_pdbqt, fld_file, output_basename)


def main():
    if not os.path.isdir(ROOT_DIR):
        raise FileNotFoundError(f"[ERROR] Cannot find docking root directory: {ROOT_DIR}")

    kinase_dirs = sorted([f.path for f in os.scandir(ROOT_DIR) if f.is_dir()])
    print(f"[INFO] Found {len(kinase_dirs)} kinase folders.")
    for kinase_dir in kinase_dirs:
        try:
            run_ensemble_docking(kinase_dir)
        except Exception as e:
            print(f"[ERROR] Docking failed for {kinase_dir}: {e}")


if __name__ == "__main__":
    main()
