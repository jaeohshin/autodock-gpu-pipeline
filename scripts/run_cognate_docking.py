"""
2025/05/20
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

ROOT_DIR = "../cognate"

def run_cognate_docking(kinase_dir):
    print(f"\n[INFO] Processing: {kinase_dir}")

    receptor_pdb = os.path.join(kinase_dir, "receptor", "receptor.pdb")
    ligand_file  = os.path.join(kinase_dir, "ligands", "ligand.pdb")
    center_file  = os.path.join(kinase_dir, "grid_center.txt")
    output_dir   = os.path.join(kinase_dir, "output")

    if not os.path.exists(receptor_pdb) or not os.path.exists(ligand_file):
        print(f"[WARN] Missing receptor.pdb or ligand.pdb in {kinase_dir}")
        return

    if not os.path.exists(center_file):
        try:
            grid_center = calc_ligand_center(ligand_file)
            with open(center_file, "w") as f:
                f.write(f"{grid_center[0]:.3f} {grid_center[1]:.3f} {grid_center[2]:.3f}\n")
            print(f"[INFO] grid_center.txt generated for {kinase_dir}")
        except Exception as e:
            print(f"[ERROR] Failed to compute grid center for {kinase_dir}: {e}")
            return
    else:
        grid_center = read_grid_center(center_file)

    os.makedirs(output_dir, exist_ok=True)

    receptor_pdbqt = os.path.join(output_dir, "receptor.pdbqt")
    prepare_receptor(receptor_pdb, receptor_pdbqt)

    ligand_pdbqt = os.path.join(output_dir, "ligand.pdbqt")
    prepare_ligand(ligand_file, ligand_pdbqt)

    if not os.path.exists(ligand_pdbqt):
        raise RuntimeError(f"[ERROR] Ligand PDBQT not generated: {ligand_pdbqt}")

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
            run_cognate_docking(kinase_dir)
        except Exception as e:
            print(f"[ERROR] Docking failed for {kinase_dir}: {e}")

if __name__ == "__main__":
    main()
