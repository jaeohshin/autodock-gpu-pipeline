"""
2025/05/20
Jaeoh Shin, Korea Institute for Advanced Study.
"""
# docking_utils.py
import os
import subprocess
import numpy as np
import shutil 
from Bio.PDB import PDBParser

# === CONFIGURATION ===
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
        return any(line.startswith(("ATOM", "HETATM")) and line[16] not in (" ", "A") for line in f)

def remove_waters(pdb_in, pdb_out):
    with open(pdb_in, 'r') as fin, open(pdb_out, 'w') as fout:
        for line in fin:
            if not (line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == "HOH"):
                fout.write(line)

def extract_atom_types_from_dir(ligand_dir):
    atom_types = set()
    for fname in os.listdir(ligand_dir):
        if not fname.endswith(".pdbqt"):
            continue
        path = os.path.join(ligand_dir, fname)
        with open(path) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    if len(line) >= 78:
                        atom_type = line[77:79].strip()
                        if atom_type:
                            atom_types.add(atom_type)
                    
    return sorted(atom_types)

def get_atom_types_from_ligand(pdbqt_path):
    atom_types = set()
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_type = line[77:79].strip()  # correct fixed-width column for type
                if atom_type:
                    atom_types.add(atom_type)
    return sorted(atom_types)


def read_grid_center(path):
    with open(path, 'r') as f:
        return tuple(map(float, f.readline().strip().split()))

def calc_ligand_center(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("lig", pdb_file)
    coords = [atom.coord for atom in structure.get_atoms()]
    if not coords:
        raise ValueError(f"[ERROR] No atoms found in ligand PDB: {pdb_file}")
    return tuple(np.mean(coords, axis=0))

# === Docking Steps ===
def prepare_receptor(pdb_file, output_pdbqt):
    pdb_file_abs = os.path.abspath(pdb_file)
    receptor_dir = os.path.dirname(pdb_file_abs)
    receptor_stem = os.path.splitext(os.path.basename(pdb_file_abs))[0]

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


def patch_gpf_center(gpf_path, center): ## This code is need to fix problem with prepare_gpf4.py function's problem with setting grid center.
    with open(gpf_path, 'r') as f:
        lines = f.readlines()

    center_line = f"gridcenter {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}\n"

    with open(gpf_path, 'w') as f:
        for line in lines:
            if line.strip().startswith("gridcenter"):
                f.write(center_line)
            else:
                f.write(line)
                

def get_receptor_idx(gpf_path):
    # Extract receptor index like 0001 from receptor_0001.gpf
    base = os.path.basename(gpf_path)
    return base.split("_")[1].split(".")[0]

def patch_gpf_maps(gpf_path, atom_types):
    with open(gpf_path, 'r') as f:
        lines = f.readlines()

    patched_lines = []
    receptor_idx = get_receptor_idx(gpf_path)
    map_lines = [f"map receptor_{receptor_idx}.{atom}.map\n" for atom in atom_types if atom not in ['e', 'd']]

    for line in lines:
        patched_lines.append(line)
        if line.strip().startswith("ligand_types"):
            patched_lines.append(" ".join(["ligand_types"] + atom_types) + "\n")
        if line.strip().startswith("gridcenter"):
            patched_lines.extend(map_lines)

    with open(gpf_path, 'w') as f:
        f.writelines(patched_lines)


                
def generate_gpf(ligand_pdbqt, receptor_pdbqt, output_gpf, center, size, atom_types):
    output_dir = os.path.dirname(output_gpf)
    os.makedirs(output_dir, exist_ok=True)
    types_str = ",".join(atom_types)

    # === Prepare Ligand ===
    ligand_src = os.path.abspath(ligand_pdbqt)
    ligand_dst = os.path.join(output_dir, os.path.basename(ligand_pdbqt))
    if not os.path.exists(ligand_dst):
        shutil.copy(ligand_src, ligand_dst)

    # === Prepare Receptor ===
    receptor_src = os.path.abspath(receptor_pdbqt)
    receptor_dst = os.path.join(output_dir, os.path.basename(receptor_pdbqt))
    if not os.path.exists(receptor_dst):
        shutil.copy(receptor_src, receptor_dst)

    # === Run GPF Generation ===
    ligand_name = os.path.basename(ligand_dst)
    receptor_name = os.path.basename(receptor_dst)
    gpf_name = os.path.basename(output_gpf)
    center_str = ",".join(map(str, center))
    size_str = ",".join(map(str, size))
    
    types_str = " ".join(atom_types)
    run_cmd(
        f"cd {output_dir} && {PREPARE_GPF} "
        f"-l {ligand_name} -r {receptor_name} -o {gpf_name} "
        f"-p npts={size_str} "
        f"-p gridcenter={center_str} "
        f" -p ligand_types=Br,A,C,Cl,F,HD,I,N,NA,OA,P,S,SA "
        )
    #patch_gpf_center(output_gpf, center)
    #patch_gpf_maps(output_gpf, atom_types)

#f"-p ligand_types={types_str}"


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
