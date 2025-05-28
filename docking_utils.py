"""
2025/05/20
Jaeoh Shin, Korea Institute for Advanced Study.
"""

import os
import subprocess
from Bio.PDB import PDBParser
from glob import glob

# === Paths ===
MGL_HOME = os.path.expanduser("~/programs/mgltools_x86_64Linux2_1.5.7")
AUTODOCK_GPU_BIN = "/home/jaeohshin/programs/AutoDock-GPU/bin/autodock_gpu_128wi"
AUTOGRID_BIN = "autogrid4"

PREPARE_RECEPTOR = f"{MGL_HOME}/bin/pythonsh {MGL_HOME}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
PREPARE_LIGAND   = f"{MGL_HOME}/bin/pythonsh {MGL_HOME}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
PREPARE_GPF      = f"{MGL_HOME}/bin/pythonsh {MGL_HOME}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py"
SPLIT_ALTLOCS    = f"{MGL_HOME}/bin/pythonsh {MGL_HOME}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_pdb_split_alt_confs.py"

# === Utility ===
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

# === Step 1: Receptor Prep ===
def prepare_receptor(pdb_file, output_pdbqt):
    pdb_file_abs = os.path.abspath(pdb_file)
    receptor_dir = os.path.dirname(pdb_file_abs)
    receptor_stem = os.path.splitext(os.path.basename(pdb_file_abs))[0]

    # Step 1: Remove waters before anything
    nowat_path = os.path.join(receptor_dir, f"{receptor_stem}_nowat.pdb")
    remove_waters(pdb_file_abs, nowat_path)

    # Step 2: Use water-free file going forward
    pdb_file_abs = nowat_path

    # Step 3: If altlocs exist, split conformers
    if has_altlocs(pdb_file_abs):
        print(f"[INFO] Alternate locations detected in {pdb_file}. Running split...")
        split_stem = f"{receptor_stem}_alt"
        run_cmd(f"cd {receptor_dir} && {SPLIT_ALTLOCS} -r {os.path.basename(pdb_file_abs)} -o {split_stem}")
        alt_file = os.path.join(receptor_dir, f"{split_stem}_A.pdb")
        if not os.path.exists(alt_file):
            raise FileNotFoundError(f"[ERROR] Expected alternate file not found: {alt_file}")
        pdb_file_abs = alt_file
        print(f"[INFO] Using alternate structure: {pdb_file_abs}")

    # Step 4: Prepare receptor
    output_pdbqt_abs = os.path.abspath(output_pdbqt)
    run_cmd(f"{PREPARE_RECEPTOR} -r {pdb_file_abs} -o {output_pdbqt_abs} -A checkhydrogens -U nphs,lps,waters")

# === Step 2: Ligand Prep — Only MGLTools ===
def prepare_ligand(pdb_file, output_pdbqt):
    run_cmd(f"{PREPARE_LIGAND} -l {pdb_file} -o {output_pdbqt} -A hydrogens")

# === Step 3: GPF ===
def generate_gpf(ligand_pdbqt, receptor_pdbqt, output_gpf, center, size):
    center_str = f"{center[0]},{center[1]},{center[2]}"
    size_str = f"{size[0]},{size[1]},{size[2]}"
    cmd = (
        f"{PREPARE_GPF} -l {ligand_pdbqt} -r {receptor_pdbqt} -y -o {output_gpf} "
        f"-p npts={size_str} -p gridcenter={center_str}"
    )
    run_cmd(cmd)

# === Step 4: AutoGrid ===
def run_autogrid(gpf_file):
    run_cmd(f"{AUTOGRID_BIN} -p {gpf_file} -l {os.path.dirname(gpf_file)}/autogrid.log")

# === Step 5: AutoDock-GPU ===
def run_docking(lig_pdbqt, fld_file, output_basename):
    cmd = (
        f"{AUTODOCK_GPU_BIN} --ffile {fld_file} --lfile {lig_pdbqt} "
        f"--nrun 10 --nev 2500000 --ngen 42000 --heuristics 1 --autostop 1 --lsrat 100 "
        f"--resnam {output_basename}"
    )
    run_cmd(cmd)

# === Optional: Ligand center for docking box ===
def calc_ligand_center(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("LIG", pdb_file)
    coords = [atom.coord for atom in structure.get_atoms()]
    center = sum(coords) / len(coords)
    return tuple(center)
