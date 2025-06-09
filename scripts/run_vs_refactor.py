"""
2025/06/04
Jaeoh Shin, Korea Institute for Advanced Study.

python run_vs_refactor.py --kinase abl1 
or
python run_vs_refactor.py --all
    
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


def get_atom_types(kinase_lig_dir):
    print(f"[DEBUG] Scanning ligand directory: {kinase_lig_dir}")
    atom_types = set()
    for mode in ["actives", "decoys"]:
        mode_dir = os.path.join(kinase_lig_dir, mode)
        if os.path.isdir(mode_dir):
            atom_types.update(extract_atom_types_from_dir(mode_dir))
    result = sorted(atom_types)
    print(f"[DEBUG] Extracted atom types: {result}")
    return result

def find_example_ligand(kinase_lig_dir):
    print(f"[DEBUG] Searching for example ligand in: {kinase_lig_dir}")
    for mode in ["actives", "decoys"]:
        mode_dir = os.path.join(kinase_lig_dir, mode)
        if os.path.isdir(mode_dir):
            ligands = [f for f in os.listdir(mode_dir) if f.endswith(".pdbqt")]
            if ligands:
                path = os.path.join(mode_dir, ligands[0])
                print(f"[DEBUG] Found example ligand: {path}")
                return path
    raise RuntimeError("No ligand files found in actives/decoys")


def generate_grid_if_needed(receptor_idx, receptor_pdb, receptor_pdbqt, center_file, fld_dir, gpf_file, example_ligand_path, atom_types):
    if not os.path.exists(center_file):
        print(f"[WARN] Missing center: {center_file}")
        return False

    center = read_grid_center(center_file)
    nx, ny, nz = GRID_SIZE
    shifts = [(0.0, 0.0, 0.0), (0.1, 0.1, 0.1), (-0.1, -0.1, -0.1), (0.05, -0.05, 0.05)]
    for attempt_idx, shift in enumerate(shifts):
        current_center = [c + s for c, s in zip(center, shift)]
        if not os.path.exists(receptor_pdbqt):
            print(f"[DEBUG] Preparing receptor: {receptor_pdb} → {receptor_pdbqt}")
            prepare_receptor(receptor_pdb, receptor_pdbqt)
        print(f"[DEBUG] Attempting center shift {attempt_idx}: {shift}")
        generate_gpf(example_ligand_path, receptor_pdbqt, gpf_file, current_center, GRID_SIZE, atom_types)
        run_autogrid(gpf_file)
        maps_ok = all(
            is_valid_map(os.path.join(fld_dir, f"receptor_{receptor_idx}.{t}.map"), nx, ny, nz)
            for t in atom_types
        )
        if maps_ok:
            patch_gpf_ligand_types(gpf_path=gpf_file, atom_types=atom_types)
            print(f"[SUCCESS] Valid maps created for receptor_{receptor_idx}")
            return True
    print(f"[ERROR] receptor_{receptor_idx}: All {len(shifts)} center shifts failed")
    return False


def patch_gpf_ligand_types(gpf_path, atom_types):
    """
    Overwrite the 'ligand_types' line in a .gpf file with a full list of atom types.

    Parameters:
    - gpf_path (str): Path to the .gpf file.
    - atom_types (list of str): Atom types to include, e.g., ['C', 'N', 'O', 'F', ...]
    """
    patched_lines = []
    with open(gpf_path, 'r') as f:
        for line in f:
            if line.startswith("ligand_types"):
                atom_str = ' '.join(atom_types)
                patched_lines.append(f"ligand_types  {atom_str}\n")
            else:
                patched_lines.append(line)
    with open(gpf_path, 'w') as f:
        f.writelines(patched_lines)
    print(f"[INFO] Patched ligand_types in {gpf_path} with: {atom_str}")


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

    atom_types = get_atom_types(kinase_lig_dir)
    print(f"[INFO] Atom types for {kinase}: {atom_types}")
    example_ligand_path = find_example_ligand(kinase_lig_dir)

    for receptor_file in sorted(os.listdir(receptor_pdb_dir)):
        if not receptor_file.endswith(".pdb"):
            continue
        receptor_idx = os.path.splitext(receptor_file)[0].split("_")[-1]
        receptor_pdb = os.path.join(receptor_pdb_dir, receptor_file)
        receptor_pdbqt = os.path.join(receptor_pdbqt_dir, f"receptor_{receptor_idx}.pdbqt")
        center_file = os.path.join(grid_center_dir, f"receptor_{receptor_idx}.txt")
        fld_file = os.path.join(fld_dir, f"receptor_{receptor_idx}.maps.fld")
        gpf_file = fld_file.replace(".maps.fld", ".gpf")
                # Example use
        generate_grid_if_needed(receptor_idx, receptor_pdb, receptor_pdbqt, center_file, fld_dir, gpf_file, example_ligand_path, atom_types)
        
        
        

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
            print(f"[INFO] Starting docking for receptor_{receptor_idx} with {len(ligand_names)} ligands")
            tasks = [
                (ligand_name, lig_pdbqt_dir, fld_file, out_subdir, receptor_idx, kinase)
                for ligand_name in ligand_names
            ]
            with multiprocessing.Pool(processes=4) as pool:
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
