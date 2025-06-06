"""
2025/06/04
Jaeoh Shin, Korea Institute for Advanced Study.

python run_vs.py --kinase abl1 
or
python run_vs.py --all
    
"""
import os
import sys
import re
import argparse
import multiprocessing

from docking_utils import (
    prepare_receptor,
    generate_gpf,
    run_autogrid,
    run_docking,
    read_grid_center,
    GRID_SIZE,
    extract_atom_types_from_dir
)

INPUT_DIR = "../virtual_screening/input"
PREPROCESS_DIR = "../virtual_screening/preprocessed"
GRID_DIR = "../virtual_screening/grids"
OUT_DIR = "../virtual_screening/docking_output"

def dock_ligand_wrapper(args):
    ligand_name, lig_pdbqt_dir, fld_file, out_subdir, receptor_idx, kinase = args

    ligand_pdbqt = os.path.join(lig_pdbqt_dir, ligand_name)
    if not os.path.exists(ligand_pdbqt):
        print(f"[WARN] Missing ligand: {ligand_pdbqt}")
        return

    out_prefix = os.path.join(out_subdir, os.path.splitext(ligand_name)[0])
    dlg_file = out_prefix + ".dlg"

    if os.path.exists(dlg_file):
        print(f"[SKIP] Already docked: {dlg_file}")
        return

    print(f"[DOCK] {kinase} | receptor_{receptor_idx} | {ligand_name}")
    run_docking(ligand_pdbqt, fld_file, out_prefix)


def is_valid_map(map_path, nx, ny, nz):
    expected_lines = (nx + 1) * (ny + 1) * (nz + 1) + 6
    try:
        with open(map_path) as f:
            lines = f.readlines()
            if len(lines) != expected_lines:
                print(f"[ERROR] Invalid line count in {map_path}: {len(lines)} (expected {expected_lines})")
                return False
            for i, line in enumerate(lines[6:], 7):  # skip header
                if not re.match(r'^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$', line.strip()):
                    print(f"[ERROR] Malformed float in {map_path} at line {i}: {line.strip()}")
                    return False
    except Exception as e:
        print(f"[ERROR] Failed reading {map_path}: {e}")
        return False
    return True


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

    # --- Extract ligand atom types once ---
    atom_types = set()
    for mode in ["actives", "decoys"]:
        mode_dir = os.path.join(kinase_lig_dir, mode)
        if os.path.isdir(mode_dir):
            atom_types.update(extract_atom_types_from_dir(mode_dir))
    atom_types = sorted(atom_types)
    
    
    print(f"[INFO] Atom types for {kinase}: {atom_types}")    
    
    # --- PRE-GENERATE fld/gpf files ---

    # Use first ligand found across actives/decoys
    example_ligand_path = None
    for mode in ["actives", "decoys"]:
        mode_dir = os.path.join(kinase_lig_dir, mode)
        if not os.path.isdir(mode_dir):
            continue
        ligands = [f for f in os.listdir(mode_dir) if f.endswith(".pdbqt")]
        if ligands:
            example_ligand_path = os.path.join(mode_dir, ligands[0])
            break

    if example_ligand_path is None:
        raise RuntimeError(f"No ligand files found in actives/decoys for {kinase}")

    for receptor_file in sorted(os.listdir(receptor_pdb_dir)):

        print(f"[CHECK] Scanning: {receptor_file}")
        if not receptor_file.endswith(".pdb"):
            continue

        receptor_idx = os.path.splitext(receptor_file)[0].split("_")[-1]
        receptor_pdb = os.path.join(receptor_pdb_dir, receptor_file)
        receptor_pdbqt = os.path.join(receptor_pdbqt_dir, f"receptor_{receptor_idx}.pdbqt")
        center_file = os.path.join(grid_center_dir, f"receptor_{receptor_idx}.txt")
        fld_file = os.path.join(fld_dir, f"receptor_{receptor_idx}.maps.fld")
        gpf_file = fld_file.replace(".maps.fld", ".gpf")

        if not os.path.exists(center_file):
            print(f"[WARN] Missing center: {center_file}")
            continue

        nx, ny, nz = GRID_SIZE

        missing_or_invalid_maps = False
        for atom_type in atom_types:
            map_path = os.path.join(fld_dir, f"receptor_{receptor_idx}.{atom_type}.map")
            if not is_valid_map(map_path, nx, ny, nz):
                missing_or_invalid_maps = True
                break

        if missing_or_invalid_maps:
            center = read_grid_center(center_file)
            
            shifts = [(0.0, 0.0, 0.0), (0.1, 0.1, 0.1), (-0.1, -0.1, -0.1), (0.05, -0.05, 0.05)]
            for attempt_idx, shift in enumerate(shifts):
                current_center = [c + s for c, s in zip(center, shift)]
                if not os.path.exists(receptor_pdbqt):
                    print(f"[INFO] Preparing receptor: {receptor_pdb} → {receptor_pdbqt}")
                    prepare_receptor(receptor_pdb, receptor_pdbqt)
                    
                generate_gpf(example_ligand_path, receptor_pdbqt, gpf_file, current_center, GRID_SIZE, atom_types)
                run_autogrid(gpf_file)

                maps_ok = all(
                    is_valid_map(os.path.join(fld_dir, f"receptor_{receptor_idx}.{t}.map"), nx, ny, nz)
                    for t in atom_types
                )
                if maps_ok:
                    print(f"[GEN] Created fld for receptor_{receptor_idx} (attempt {attempt_idx}, shift={shift})")
                    break
            else:
                print(f"[FAIL] receptor_{receptor_idx}: All {len(shifts)} center shifts failed to yield valid maps.")
                    
                    

    # --- DOCKING ---
    for mode in ["actives", "decoys"]:
        list_path = os.path.join(kinase_lig_dir, f"ligands_{mode}.list")
        lig_pdbqt_dir = os.path.join(kinase_lig_dir, mode)

        with open(list_path) as f:
            ligand_names = [x.strip() for x in f if x.strip()]

        for receptor_file in sorted(os.listdir(receptor_pdb_dir)):
            if not receptor_file.endswith(".pdb"):
                continue

            receptor_idx = os.path.splitext(receptor_file)[0].split("_")[-1]
            fld_file = os.path.join(fld_dir, f"receptor_{receptor_idx}.maps.fld")
            out_subdir = os.path.join(out_dir, f"receptor_{receptor_idx}")
            os.makedirs(out_subdir, exist_ok=True)

            tasks = [
                (ligand_name, lig_pdbqt_dir, fld_file, out_subdir, receptor_idx, kinase)
                for ligand_name in ligand_names
            ]

            with multiprocessing.Pool(processes=1) as pool:
                pool.map(dock_ligand_wrapper, tasks)


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