"""
2025/07/10
Unified Virtual Screening Script
Author: Jaeoh Shin, Korea Institute for Advanced Study

Supports two use cases with one script:
  1. Crystal structure docking (vs_crystal/) with one receptor per kinase
  2. Ensemble docking (virtual_screening/) using multiple BioEmu-generated structures

Usage:
  python run_vs.py --project vs_crystal --mode crystal --kinase abl1
  python run_vs.py --project virtual_screening --mode ensemble --all --nprocs 16
  python run_vs.py --project vs_crystal --mode crystal --kinase abl1 --force-grid  # Force grid regeneration

Required Directory Layout:
  ../{project}/input/
    ├── receptors/<kinase>/receptor.pdb (crystal) or receptor_*.pdb (ensemble)
    ├── ligands/ or ligands_sdf2meeko/<kinase>/{actives,decoys}/
    └── kinase_list.txt

  ../{project}/preprocessed/
    ├── receptors_pdbqt/<kinase>/receptor_*.pdbqt
    ├── ligands_pdbqt/<kinase>/{actives,decoys}/
    └── grid_centers/<kinase>/receptor_*.txt

  ../{project}/grids/<kinase>/ → contains .map files and .fld
  ../{project}/docking_output/<kinase>/ → docking results
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

def is_valid_map(map_path, *args, **kwargs):
    return os.path.exists(map_path)

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
    for mode in ["actives", "decoys"]:
        mode_dir = os.path.join(kinase_lig_dir, mode)
        if os.path.isdir(mode_dir):
            ligands = [f for f in os.listdir(mode_dir) if f.endswith(".pdbqt")]
            if ligands:
                return os.path.join(mode_dir, ligands[0])
    raise RuntimeError("No ligand files found in actives/decoys")

def patch_gpf_ligand_types(gpf_path, atom_types):
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

def check_grid_files_exist(fld_dir, receptor_idx, atom_types):
    """Check if all required grid files exist for a receptor"""
    # Check for the main .fld file
    fld_file = os.path.join(fld_dir, f"receptor_{receptor_idx}.maps.fld")
    if not os.path.exists(fld_file):
        return False
    
    # Check for the .gpf file
    gpf_file = os.path.join(fld_dir, f"receptor_{receptor_idx}.gpf")
    if not os.path.exists(gpf_file):
        return False
    
    # Check for all required .map files
    for atom_type in atom_types:
        map_file = os.path.join(fld_dir, f"receptor_{receptor_idx}.{atom_type}.map")
        if not os.path.exists(map_file):
            return False
    
    # Check for electrostatic and desolvation maps
    e_map = os.path.join(fld_dir, f"receptor_{receptor_idx}.e.map")
    d_map = os.path.join(fld_dir, f"receptor_{receptor_idx}.d.map")
    if not os.path.exists(e_map) or not os.path.exists(d_map):
        return False
    
    return True

def generate_grid_if_needed(receptor_idx, receptor_pdb, receptor_pdbqt, center_file, fld_dir, gpf_file, example_ligand_path, atom_types, force_grid=False):
    """Generate grid files only if they don't exist or if force_grid is True"""
    
    if not os.path.exists(center_file):
        print(f"[WARN] Missing center: {center_file}")
        return False
    
    # Check if grid files already exist
    if not force_grid and check_grid_files_exist(fld_dir, receptor_idx, atom_types):
        print(f"[SKIP] Grid files already exist for receptor_{receptor_idx}. Use --force-grid to regenerate.")
        return True
    
    # If we reach here, either force_grid is True or files are missing
    if force_grid:
        print(f"[INFO] Force regenerating grid for receptor_{receptor_idx}")
    else:
        print(f"[INFO] Generating missing grid files for receptor_{receptor_idx}")
    
    center = read_grid_center(center_file)
    nx, ny, nz = GRID_SIZE
    shifts = [(0.0, 0.0, 0.0), (0.1, 0.1, 0.1), (-0.1, -0.1, -0.1), (0.05, -0.05, 0.05)]
    
    for attempt_idx, shift in enumerate(shifts):
        current_center = [c + s for c, s in zip(center, shift)]
        
        # Prepare receptor PDBQT if needed
        if not os.path.exists(receptor_pdbqt):
            print(f"[INFO] Preparing receptor PDBQT: {receptor_pdbqt}")
            prepare_receptor(receptor_pdb, receptor_pdbqt)
        
        # Generate GPF and run autogrid
        generate_gpf(example_ligand_path, receptor_pdbqt, gpf_file, current_center, GRID_SIZE, atom_types)
        run_autogrid(gpf_file)
        
        # Check if all maps were created successfully
        maps_ok = all(is_valid_map(os.path.join(fld_dir, f"receptor_{receptor_idx}.{t}.map")) for t in atom_types)
        
        if maps_ok:
            patch_gpf_ligand_types(gpf_path=gpf_file, atom_types=atom_types)
            print(f"[SUCCESS] Valid maps created for receptor_{receptor_idx} (attempt {attempt_idx + 1})")
            return True
    
    print(f"[ERROR] receptor_{receptor_idx}: All {len(shifts)} center shifts failed")
    return False

def run_vs_for_kinase(kinase, args, INPUT_DIR, PREPROCESS_DIR, GRID_DIR, OUT_DIR):
    ROOT = f"../{args.project}"
    input_dir = os.path.join(ROOT, "input")
    pre_dir = os.path.join(ROOT, "preprocessed")
    grid_dir = os.path.join(ROOT, "grids")
    out_dir = os.path.join(ROOT, "docking_output")

    # Define receptor directories
    receptor_dir = os.path.join(input_dir, "receptors", kinase)
    receptor_pdbqt_dir = os.path.join(pre_dir, "receptors_pdbqt", kinase)
    center_dir = os.path.join(pre_dir, "grid_centers", kinase)
    fld_dir = os.path.join(grid_dir, kinase)
    out_dir = os.path.join(out_dir, kinase)

    os.makedirs(receptor_pdbqt_dir, exist_ok=True)
    os.makedirs(fld_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    # Ligand directory logic (fallback for ensemble mode)
    if args.mode == "crystal":
        lig_dir = os.path.join(pre_dir, "ligands_pdbqt", kinase)
    else:  # ensemble mode
        lig_dir = os.path.join(pre_dir, "ligands_pdbqt", kinase)
        if not os.path.isdir(lig_dir):
            fallback_dir = os.path.join(ROOT, "input", "ligands_sdf2meeko", kinase)
            if os.path.isdir(fallback_dir):
                lig_dir = fallback_dir
                print(f"[INFO] [ensemble] Using fallback ligand path: {lig_dir}")
            else:
                raise FileNotFoundError(f"No ligand directory found for {kinase}")

    # Automatically create ligand list files if missing
    if args.mode == "ensemble":
        for mode in ["actives", "decoys"]:
            mode_dir = os.path.join(lig_dir, mode)
            list_path = os.path.join(lig_dir, f"ligands_{mode}.list")
            if not os.path.exists(list_path) and os.path.isdir(mode_dir):
                ligand_names = sorted(f for f in os.listdir(mode_dir) if f.endswith(".pdbqt"))
                with open(list_path, "w") as f:
                    f.write("\n".join(ligand_names) + "\n")
                print(f"[INFO] Created ligand list: {list_path}")

    # Atom types for grid map generation
    atom_types = get_atom_types(lig_dir)

    print(f"[DEBUG] Looking for ligands in: {lig_dir}")
    for mode in ["actives", "decoys"]:
        mode_dir = os.path.join(lig_dir, mode)
        print(f"→ Checking: {mode_dir}")
        if not os.path.isdir(mode_dir):
            print(f"[WARN] Missing dir: {mode_dir}")
        else:
            print(f"[OK] Found {len(os.listdir(mode_dir))} files in {mode_dir}")

    example_ligand = find_example_ligand(lig_dir)

    receptor_files = sorted(f for f in os.listdir(receptor_dir) if f.endswith(".pdb"))
    for receptor_file in receptor_files:
        receptor_pdb = os.path.join(receptor_dir, receptor_file)
        
        if args.mode == "crystal":
            receptor_idx = "crystal"
            receptor_pdbqt = os.path.join(receptor_pdbqt_dir, "receptor.pdbqt")
        else:
            receptor_idx = os.path.splitext(receptor_file)[0].split("_")[-1]
            receptor_pdbqt = os.path.join(receptor_pdbqt_dir, f"receptor_{receptor_idx}.pdbqt")

        center_file = os.path.join(center_dir, f"receptor_{receptor_idx}.txt")
        fld_file = os.path.join(fld_dir, f"receptor_{receptor_idx}.maps.fld")
        gpf_file = fld_file.replace(".maps.fld", ".gpf")

        # Pass the force_grid flag from args
        generate_grid_if_needed(
            receptor_idx, receptor_pdb, receptor_pdbqt, center_file, 
            fld_dir, gpf_file, example_ligand, atom_types, 
            force_grid=args.force_grid
        )

    # Run docking
    for mode in ["actives", "decoys"]:
        list_path = os.path.join(lig_dir, f"ligands_{mode}.list")
        lig_pdbqt_dir = os.path.join(lig_dir, mode)
        with open(list_path) as f:
            ligand_names = [x.strip() for x in f if x.strip()]
        for receptor_file in receptor_files:
            receptor_idx = "crystal" if args.mode == "crystal" else os.path.splitext(receptor_file)[0].split("_")[-1]
            fld_file = os.path.join(fld_dir, f"receptor_{receptor_idx}.maps.fld")
            out_subdir = out_dir if args.mode == "crystal" else os.path.join(out_dir, f"receptor_{receptor_idx}")
            os.makedirs(out_subdir, exist_ok=True)
            tasks = [
                (ligand_name, lig_pdbqt_dir, fld_file, out_subdir, receptor_idx, kinase)
                for ligand_name in ligand_names
            ]
            with multiprocessing.Pool(processes=args.nprocs) as pool:
                pool.map(dock_ligand_wrapper, tasks)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=str, default="virtual_screening", choices=["virtual_screening", "vs_crystal"])
    parser.add_argument("--mode", type=str, default="ensemble", choices=["ensemble", "crystal"])
    parser.add_argument("--nprocs", type=int, default=4)
    parser.add_argument("--force-grid", action="store_true", 
                        help="Force regeneration of grid files even if they exist")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--kinase", type=str)
    group.add_argument("--all", action="store_true")
    args = parser.parse_args()
    
    BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "vs_crystal"))
    INPUT_DIR = os.path.join(BASE_DIR, "input")
    PREPROCESS_DIR = os.path.join(BASE_DIR, "preprocessed")
    GRID_DIR = os.path.join(BASE_DIR, "grids")
    OUT_DIR = os.path.join(BASE_DIR, "docking_output")

    if args.kinase:
        run_vs_for_kinase(args.kinase.lower(), args, INPUT_DIR, PREPROCESS_DIR, GRID_DIR, OUT_DIR)
    elif args.all:
        kinase_list_path = os.path.join(f"../{args.project}", "input", "kinase.txt")
        with open(kinase_list_path) as f:
            kinases = [x.strip().lower() for x in f if x.strip()]
        for kinase in kinases:
            run_vs_for_kinase(kinase, args, INPUT_DIR, PREPROCESS_DIR, GRID_DIR, OUT_DIR)

if __name__ == "__main__":
    main()