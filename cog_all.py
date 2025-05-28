"""
Docking automation script for multiple kinases.

2025/05/28
Jaeoh Shin, Korea Institute for Advanced Study.
"""

import os
import subprocess
import numpy as np
from Bio.PDB import PDBParser

# === CONFIGURATION ===
ROOT_DIR = "./cognate"
GRID_SIZE = (60, 60, 60)
AUTODOCK_GPU_BIN = "/home/jaeohshin/programs/AutoDock-GPU/bin/autodock_gpu_128wi"

MGL_HOME = os.path.expanduser("~/programs/mgltools_x86_64Linux2_1.5.7")
PREPARE_RECEPTOR = f"{MGL_HOME}/bin/pythonsh {MGL_HOME}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
PREPARE_LIGAND   = f"{MGL_HOME}/bin/pythonsh {MGL_HOME}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
PREPARE_GPF      = f"{MGL_HOME}/bin/pythonsh {MGL_HOME}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py"
SPLIT_ALTLOCS    = f"{MGL_HOME}/bin/pythonsh {MGL_HOME}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_pdb_split_alt_confs.py"
AUTOGRID_BIN     = "autogrid4"

# === Utilities ===
def run_cmd(cmd):
    print(f"[RUN] {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def has_altlocs(pdb_path):
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and line[16] not in (" ", "A"):
                return True
    return False

def remove_waters(pdb_in, pdb_out):
    with open(pdb_in, 'r') as fin, open(pdb_out, 'w') as fout:
        for line in fin:
            if not (line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == "HOH"):
                fout.write(line)

def read_grid_center(path):
    with open(path, 'r') as f:
        return tuple(map(float, f.readline().strip().split()))

def calc_ligand_center(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("lig", pdb_file)
    coords = [atom.coord for atom in structure.get_atoms()]
    if not coords:
        raise ValueError(f"[ERROR] No atoms found in ligand PDB: {pdb_file}")
    center = np.mean(coords, axis=0)
    return tuple(center)

# === Preparation Steps ===
def prepare_receptor(pdb_file, output_pdbqt):
    pdb_file_abs = os.path.abspath(pdb_file)
    receptor_dir = os.path.dirname(pdb_file_abs)
    receptor_stem = os.path.splitext(os.path.basename(pdb_file_abs))[0]
    nowat_path = os.path.join(receptor_dir, f"{receptor_stem}_nowat.pdb")
    remove_waters(pdb_file_abs, nowat_path)
    pdb_file_abs = nowat_path

    if has_altlocs(pdb_file_abs):
        print(f"[INFO] Alternate locations detected in {pdb_file}. Running split...")
        split_stem = f"{receptor_stem}_alt"
        run_cmd(f"cd {receptor_dir} && {SPLIT_ALTLOCS} -r {os.path.basename(pdb_file_abs)} -o {split_stem}")
        alt_file = os.path.join(receptor_dir, f"{split_stem}_A.pdb")
        if not os.path.exists(alt_file):
            raise FileNotFoundError(f"[ERROR] Altloc file not found: {alt_file}")
        pdb_file_abs = alt_file
        print(f"[INFO] Using alternate: {pdb_file_abs}")

    run_cmd(f"{PREPARE_RECEPTOR} -r {pdb_file_abs} -o {output_pdbqt} -A checkhydrogens -U nphs,lps,waters")

def prepare_ligand(pdb_file, output_pdbqt):
    lig_dir = os.path.dirname(os.path.abspath(pdb_file))
    lig_name = os.path.basename(pdb_file)
    output_pdbqt_abs = os.path.abspath(output_pdbqt)
    run_cmd(f"cd {lig_dir} && {PREPARE_LIGAND} -l {lig_name} -o {output_pdbqt_abs} -A hydrogens")


def generate_gpf(ligand_pdbqt, receptor_pdbqt, output_gpf, center, size):
    output_dir = os.path.dirname(output_gpf)
    ligand_name = os.path.basename(ligand_pdbqt)
    receptor_name = os.path.basename(receptor_pdbqt)
    gpf_name = os.path.basename(output_gpf)
    center_str = ",".join(map(str, center))
    size_str = ",".join(map(str, size))
    run_cmd(
        f"cd {output_dir} && {PREPARE_GPF} "
        f"-l {ligand_name} -r {receptor_name} -y -o {gpf_name} "
        f"-p npts={size_str} -p gridcenter={center_str}"
    )

def run_autogrid(gpf_file):
    fld_dir = os.path.dirname(gpf_file)
    gpf_name = os.path.basename(gpf_file)
    run_cmd(f"cd {fld_dir} && {AUTOGRID_BIN} -p {gpf_name} -l autogrid.log")

def run_docking(lig_pdbqt, fld_file, output_basename):
    run_cmd(
        f"{AUTODOCK_GPU_BIN} "
        f"--ffile {fld_file} "
        f"--lfile {lig_pdbqt} "
        f"--nrun 50 --nev 2500000 --ngen 42000 "
        f"--heuristics 1 --autostop 1 --lsrat 100 "
        f"--resnam {output_basename}"
    )

# === Main Docking Workflow for One Kinase ===
def run_cognate_docking(kinase_dir):
    print(f"\n[INFO] Processing: {kinase_dir}")

    receptor_pdb = os.path.join(kinase_dir, "receptor", "receptor.pdb")
    ligand_file  = os.path.join(kinase_dir, "ligands", "ligand.pdb")
    center_file  = os.path.join(kinase_dir, "grid_center.txt")
    output_dir   = os.path.join(kinase_dir, "output")

    if not os.path.exists(receptor_pdb) or not os.path.exists(ligand_file):
        print(f"[WARN] Missing receptor.pdb or ligand.pdb in {kinase_dir}")
        return

    # Auto-generate grid center if missing
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

# === Batch Execution ===
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
