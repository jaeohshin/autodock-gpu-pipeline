"""
2025/06/04
Jaeoh Shin, Korea Institute for Advanced Study.

python run_vs.py --kinase abl1 
or
python run_vs.py --all
    
"""
import os
import sys
import argparse

from docking_utils import (
    prepare_receptor, generate_gpf, run_autogrid, run_docking, read_grid_center, GRID_SIZE
)

INPUT_DIR = "../virtual_screening/input"
PREPROCESS_DIR = "../virtual_screening/preprocessed"
GRID_DIR = "../virtual_screening/grids"
OUT_DIR = "../virtual_screening/docking_output"


def run_vs_for_kinase(kinase):
    kinase_lig_dir = os.path.join(INPUT_DIR, "ligands", kinase)
    receptor_pdb_dir = os.path.join(INPUT_DIR, "receptors", kinase)
    receptor_pdbqt_dir = os.path.join(PREPROCESS_DIR, "receptors_pdbqt", kinase)
    grid_center_dir = os.path.join(PREPROCESS_DIR, "grid_centers", kinase)
    fld_dir = os.path.join(GRID_DIR, kinase)
    out_dir = os.path.join(OUT_DIR, kinase)

    os.makedirs(receptor_pdbqt_dir, exist_ok=True)
    os.makedirs(fld_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    for mode in ["actives", "decoys"]:
        list_path = os.path.join(kinase_lig_dir, f"ligands_{mode}.list")
        lig_pdbqt_dir = os.path.join(kinase_lig_dir, mode)

        with open(list_path) as f:
            ligand_names = [x.strip() for x in f if x.strip()]

        for receptor_file in sorted(os.listdir(receptor_pdb_dir)):
            if not receptor_file.endswith(".pdb"):
                continue

            receptor_idx = os.path.splitext(receptor_file)[0].split("_")[-1]

            receptor_pdb = os.path.join(receptor_pdb_dir, receptor_file)
            receptor_pdbqt = os.path.join(receptor_pdbqt_dir, f"receptor_{receptor_idx}.pdbqt")
            center_file = os.path.join(grid_center_dir, f"receptor_{receptor_idx}.txt")
            fld_file = os.path.join(fld_dir, f"receptor_{receptor_idx}.maps.fld")

            # Convert receptor if needed
            if not os.path.exists(receptor_pdbqt):
                prepare_receptor(receptor_pdb, receptor_pdbqt)

            # Read grid center
            if not os.path.exists(center_file):
                print(f"[WARN] Missing center: {center_file}")
                continue
            center = read_grid_center(center_file)

            # Generate GPF and Grid maps
            gpf_file = fld_file.replace(".maps.fld", ".gpf")
            if not os.path.exists(fld_file):
                ligand_path = os.path.join(lig_pdbqt_dir, ligand_names[0])
                generate_gpf(ligand_path, receptor_pdbqt, gpf_file, center, GRID_SIZE)
                run_autogrid(gpf_file)

            for ligand_name in ligand_names:
                ligand_pdbqt = os.path.join(lig_pdbqt_dir, ligand_name)
                if not os.path.exists(ligand_pdbqt):
                    print(f"[WARN] Missing ligand: {ligand_pdbqt}")
                    continue

                out_subdir = os.path.join(out_dir, f"receptor_{receptor_idx}")
                os.makedirs(out_subdir, exist_ok=True)

                out_prefix = os.path.join(out_subdir, os.path.splitext(ligand_name)[0])
                dlg_file = out_prefix + ".dlg"

                if os.path.exists(dlg_file):
                    print(f"[SKIP] Already docked: {dlg_file}")
                    continue

                print(f"[DOCK] {kinase} | receptor_{receptor_idx} | {ligand_name}")
                run_docking(ligand_pdbqt, fld_file, out_prefix)


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--kinase", type=str, help="Run virtual screening for a single kinase (e.g., abl1)")
    group.add_argument("--all", action="store_true", help="Run for all kinases in kinase_list.txt")
    args = parser.parse_args()

    if args.kinase:
        run_vs_for_kinase(args.kinase.lower())
    elif args.all:
        kinase_list = os.path.join(INPUT_DIR, "kinase_list.txt")
        with open(kinase_list) as f:
            kinases = [x.strip().lower() for x in f if x.strip()]
        for kinase in kinases:
            run_vs_for_kinase(kinase)


if __name__ == "__main__":
    main()
