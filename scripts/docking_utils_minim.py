# docking_utils_minim.py
import os
import subprocess
import shutil
from Bio.PDB import PDBParser
import numpy as np

# === CONFIGURATION ===
GRID_SIZE = (60, 60, 60)
AUTODOCK_GPU_BIN = "/home/jaeohshin/programs/AutoDock-GPU/bin/autodock_gpu_128wi"

MGL_HOME = os.path.expanduser("~/programs/mgltools_x86_64Linux2_1.5.7")
PREPARE_RECEPTOR = f"{MGL_HOME}/bin/pythonsh {MGL_HOME}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
PREPARE_LIGAND   = f"{MGL_HOME}/bin/pythonsh {MGL_HOME}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
PREPARE_GPF      = f"{MGL_HOME}/bin/pythonsh {MGL_HOME}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py"
AUTOGRID_BIN     = "autogrid4"

def run_cmd(cmd):
    print(f"[RUN] {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def has_altlocs(pdb_path):
    with open(pdb_path, 'r') as f:
        return any(line.startswith(("ATOM", "HETATM")) and line[16] not in (" ", "A") for line in f)

def extract_atom_types_from_dir(ligand_dir):
    atom_types = set()
    for fname in os.listdir(ligand_dir):
        if fname.endswith(".pdbqt"):
            with open(os.path.join(ligand_dir, fname)) as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")) and len(line) >= 78:
                        atom_type = line[77:79].strip()
                        if atom_type:
                            atom_types.add(atom_type)
    return sorted(atom_types)

def read_grid_center(path):
    with open(path, 'r') as f:
        return tuple(map(float, f.readline().strip().split()))

def prepare_receptor(pdb_file, output_pdbqt):
    pdb_file_abs = os.path.abspath(pdb_file)
    receptor_dir = os.path.dirname(pdb_file_abs)
    receptor_stem = os.path.splitext(os.path.basename(pdb_file_abs))[0]

    if has_altlocs(pdb_file_abs):
        split_script = f"{MGL_HOME}/bin/pythonsh {MGL_HOME}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_pdb_split_alt_confs.py"
        run_cmd(f"cd {receptor_dir} && {split_script} -r {os.path.basename(pdb_file_abs)} -o {receptor_stem}_alt")
        alt_file = os.path.join(receptor_dir, f"{receptor_stem}_alt_A.pdb")
        if not os.path.exists(alt_file):
            raise FileNotFoundError(f"[ERROR] Altloc file not found: {alt_file}")
        pdb_file_abs = alt_file
        print(f"[INFO] Using alternate: {pdb_file_abs}")

    run_cmd(f"{PREPARE_RECEPTOR} -r {pdb_file_abs} -o {output_pdbqt} -A checkhydrogens -U nphs,lps,waters")

def generate_gpf(ligand_pdbqt, receptor_pdbqt, output_gpf, center, size, atom_types):
    output_dir = os.path.dirname(output_gpf)
    os.makedirs(output_dir, exist_ok=True)

    ligand_name = os.path.basename(ligand_pdbqt)
    receptor_name = os.path.basename(receptor_pdbqt)
    gpf_name = os.path.basename(output_gpf)

    # Ensure ligand/receptor files are present in same folder
    for src in [ligand_pdbqt, receptor_pdbqt]:
        dst = os.path.join(output_dir, os.path.basename(src))
        if not os.path.exists(dst):
            shutil.copy(os.path.abspath(src), dst)

    center_str = ",".join(map(str, center))
    size_str = ",".join(map(str, size))
    types_str = " ".join(atom_types)

    # Step 1: Run prepare_gpf4.py
    run_cmd(
        f"cd {output_dir} && {PREPARE_GPF} "
        f"-l {ligand_name} -r {receptor_name} -y -o {gpf_name} "
        f"-p npts={size_str} -p gridcenter={center_str} -p ligand_types='{types_str}'"
    )

    # Step 2: Explicitly add map lines
    receptor_idx = os.path.splitext(receptor_name)[0].split('_')[-1]
    map_lines = [f"map receptor_{receptor_idx}.{atom}.map\n" for atom in atom_types]

    with open(output_gpf, 'r') as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        new_lines.append(line)
        if line.strip().startswith("gridcenter"):
            new_lines.extend(map_lines)

    with open(output_gpf, 'w') as f:
        f.writelines(new_lines)


def run_autogrid(gpf_file):
    fld_dir = os.path.dirname(gpf_file)
    gpf_name = os.path.basename(gpf_file)
    run_cmd(f"cd {fld_dir} && {AUTOGRID_BIN} -p {gpf_name} -l autogrid.log")

def run_docking(lig_pdbqt, fld_file, output_basename):
    run_cmd(
        f"{AUTODOCK_GPU_BIN} "
        f"--ffile {fld_file} "
        f"--lfile {lig_pdbqt} "
        f"--nrun 20 --nev 2500000 "
        f"--ngen 42000 --heuristics 1 --autostop 1 --lsrat 25 "
        f"--resnam {output_basename}"
    )
