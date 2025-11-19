import os
import subprocess
import numpy as np
import pdbtools
from Bio.PDB import PDBParser, PDBIO, Select

KINASE_LIST = "kinase.txt"

# === MGLTools paths ===
MGL_HOME = os.path.expanduser("~/programs/mgltools_x86_64Linux2_1.5.7")
PREPARE_RECEPTOR = f"{MGL_HOME}/bin/pythonsh {MGL_HOME}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
SPLIT_ALTLOCS = f"{MGL_HOME}/bin/pythonsh {MGL_HOME}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_pdb_split_alt_confs.py"

# === Selectors for cleaning ===
class ProteinOnlySelect(Select):
    def accept_residue(self, residue):
        return residue.id[0] == " "  # Only standard residues (ATOM)

    def accept_atom(self, atom):
        return atom.get_altloc() in (" ", "A")  # Keep first altloc

# === Utility functions ===
def run_cmd(cmd):
    print(f"[RUN] {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def has_altlocs(pdb_path):
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and line[16] not in (" ", "A"):
                return True
    return False  

def compute_ligand_center(pdb_path):
    coords = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("HETATM"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    if not coords:
        raise ValueError(f"No HETATM atoms found in {pdb_path}")
    coords = np.array(coords)
    return coords.mean(axis=0)

def clean_receptor_structure(pdb_in, pdb_out, keep_chain="A"):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", pdb_in)

    # Keep only chain A
    for model in structure:
        chains_to_remove = [chain.id for chain in model if chain.id != keep_chain]
        for cid in chains_to_remove:
            model.detach_child(cid)

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_out, ProteinOnlySelect())
    print(f"[INFO] Cleaned receptor (chain {keep_chain}) saved to {pdb_out}")

def prepare_receptor(pdb_file, output_pdbqt):
    pdb_file_abs = os.path.abspath(pdb_file)  # should be receptor.pdb (already cleaned)
    receptor_dir = os.path.dirname(pdb_file_abs)
    receptor_stem = os.path.splitext(os.path.basename(pdb_file_abs))[0]

    # if has_altlocs(pdb_file_abs):  # checks cleaned file
    #     print(f"[INFO] Alternate locations detected. Splitting...")
    #     split_stem = f"{receptor_stem}_alt"
    #     run_cmd(f"cd {receptor_dir} && {SPLIT_ALTLOCS} -r {os.path.basename(pdb_file_abs)} -o {split_stem}")
    #     alt_file = os.path.join(receptor_dir, f"{split_stem}_A.pdb")
    #     if not os.path.exists(alt_file):
    #         raise FileNotFoundError(f"[ERROR] Expected altloc file not found: {alt_file}")
    #     pdb_file_abs = alt_file  # use cleaned + altloc-filtered version
        
    altloc_filtered = os.path.join(receptor_dir, f"{receptor_stem}_altA.pdb")
    run_cmd(f"pdb_selaltloc -A {pdb_file_abs} > {altloc_filtered}")
    pdb_file_abs = altloc_filtered

    output_pdbqt_abs = os.path.abspath(output_pdbqt)
    run_cmd(f"{PREPARE_RECEPTOR} -r {pdb_file_abs} -o {output_pdbqt_abs} -A checkhydrogens -U nphs,lps,waters")



# === Main per-kinase logic ===
def process_kinase(kinase):
    orig_pdb = os.path.join(kinase, "receptor", "original.pdb")
    clean_pdb = os.path.join(kinase, "receptor", "receptor.pdb")
    receptor_pdbqt = os.path.join(kinase, "receptor", "receptor.pdbqt")
    ligand_pdb = os.path.join(kinase, "ligands", "ligand.pdb")
    center_txt = os.path.join(kinase, "grid", "grid_center.txt")

    if not os.path.exists(orig_pdb):
        print(f"[WARNING] {orig_pdb} not found. Skipping {kinase}.")
        return
    if not os.path.exists(ligand_pdb):
        print(f"[WARNING] {ligand_pdb} not found. Skipping {kinase}.")
        return

    # Step 1: Clean dirty receptor
    os.makedirs(os.path.dirname(clean_pdb), exist_ok=True)
    clean_receptor_structure(orig_pdb, clean_pdb)

    # Step 2: Compute grid center from ligand
    os.makedirs(os.path.dirname(center_txt), exist_ok=True)
    try:
        center = compute_ligand_center(ligand_pdb)
        with open(center_txt, "w") as f:
            f.write(f"{center[0]:.3f} {center[1]:.3f} {center[2]:.3f}\n")
        print(f"[INFO] {kinase}: grid center saved to {center_txt}")
    except Exception as e:
        print(f"[ERROR] Failed to compute center for {kinase}: {e}")
        return

    # Step 3: Prepare receptor.pdbqt
    prepare_receptor(clean_pdb, receptor_pdbqt)
    print(f"[INFO] {kinase}: receptor.pdbqt ready.")

if __name__ == "__main__":
    with open(KINASE_LIST) as f:
        for line in f:
            kinase = line.strip().split()[0]
            process_kinase(kinase)
