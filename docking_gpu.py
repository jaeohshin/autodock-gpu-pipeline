"""
2025/ 05/ 20
Jaeoh Shin, Korea Institute for Advanced Study.

"""
import os
import subprocess
from glob import glob

from rdkit import Chem
from meeko import MoleculePreparation

# === CONFIGURATION ===
RECEPTOR_PDB = "receptor.pdb"
LIGAND_PDB_DIR = "./ligands_pdb"
OUTPUT_DIR = "./docking_output"
AUTODOCK_GPU_BIN = "/home/jaeohshin/programs/AutoDock-GPU/bin/autodock_gpu_128wi"

MGL_HOME = os.path.expanduser("~/programs/mgltools_x86_64Linux2_1.5.7")

PREPARE_RECEPTOR = f"{MGL_HOME}/bin/pythonsh {MGL_HOME}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
PREPARE_LIGAND   = f"{MGL_HOME}/bin/pythonsh {MGL_HOME}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
PREPARE_GPF      = f"{MGL_HOME}/bin/pythonsh {MGL_HOME}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py"

AUTOGRID_BIN     = "autogrid4"

# 👇 Update this to your actual binding site center
GRID_CENTER = (10.0, 15.0, -5.0)
GRID_SIZE   = (60, 60, 60)

# === COMMAND EXECUTION ===
def run_cmd(cmd):
    print(f"[RUN] {cmd}")
    subprocess.run(cmd, shell=True, check=True)

# === STEP 1: Prepare receptor ===
def prepare_receptor(pdb_file, output_pdbqt):
    pdb_file_abs = os.path.abspath(pdb_file)
    output_pdbqt_abs = os.path.abspath(output_pdbqt)
    run_cmd(f"{PREPARE_RECEPTOR} -r {pdb_file_abs} -o {output_pdbqt_abs} -A checkhydrogens")



# === STEP 2: Prepare ligand ===
def prepare_ligand(pdb_file, output_pdbqt, use_meeko=True):
    if use_meeko:
        try:
            mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
            if mol is None:
                raise ValueError
            prep = MoleculePreparation()
            pdbqt_string = prep.prepare(mol)
            with open(output_pdbqt, "w") as f:
                f.write(pdbqt_string)
            print(f"[INFO] Ligand prepared using Meeko: {output_pdbqt}")
            return
        except:
            print(f"[WARN] Meeko failed on {pdb_file}, falling back to MGLTools...")

    # Fallback to MGLTools
    pdb_file_abs = os.path.abspath(pdb_file)
    output_pdbqt_abs = os.path.abspath(output_pdbqt)
    lig_dir = os.path.dirname(pdb_file_abs)
    lig_name = os.path.basename(pdb_file_abs)
    cmd = f"cd {lig_dir} && {PREPARE_LIGAND} -l {lig_name} -o {output_pdbqt_abs}"
    run_cmd(cmd)


# === STEP 3: Generate GPF ===
def generate_gpf(ligand_pdbqt, receptor_pdbqt, output_gpf, center, size):
    output_dir = os.path.dirname(os.path.abspath(output_gpf))
    ligand_name = os.path.basename(ligand_pdbqt)
    receptor_name = os.path.basename(receptor_pdbqt)
    gpf_name = os.path.basename(output_gpf)
    center_str = f"{center[0]},{center[1]},{center[2]}"
    size_str = f"{size[0]},{size[1]},{size[2]}"
    cmd = (
        f"cd {output_dir} && {PREPARE_GPF} "
        f"-l {ligand_name} -r {receptor_name} -y -o {gpf_name} "
        f"-p npts={size_str} -p gridcenter={center_str}"
    )
    run_cmd(cmd)

# === STEP 4: Run AutoGrid ===
def run_autogrid(gpf_file):
    fld_dir = os.path.dirname(os.path.abspath(gpf_file))
    gpf_name = os.path.basename(gpf_file)
    cmd = f"cd {fld_dir} && {AUTOGRID_BIN} -p {gpf_name} -l autogrid.log"
    run_cmd(cmd)

# === STEP 5: Run AutoDock-GPU ===
def run_docking(lig_pdbqt, fld_file, output_basename):
    fld_file_abs = os.path.abspath(fld_file)
    lig_pdbqt_abs = os.path.abspath(lig_pdbqt)
    output_base_abs = os.path.abspath(output_basename)
    cmd = (
        f"{AUTODOCK_GPU_BIN} "
        f"--ffile {fld_file_abs} "
        f"--lfile {lig_pdbqt_abs} "
        f"--nrun 50 "
        f"--resnam {output_base_abs}"
    )
    run_cmd(cmd)

# === MAIN PIPELINE ===
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    receptor_pdbqt = os.path.join(OUTPUT_DIR, "receptor.pdbqt")
    prepare_receptor(RECEPTOR_PDB, receptor_pdbqt)

    ligand_pdbs = glob(os.path.join(LIGAND_PDB_DIR, "*.pdb"))
    ligand_pdbqts = []

    for lig_pdb in ligand_pdbs:
        name = os.path.splitext(os.path.basename(lig_pdb))[0]
        lig_pdbqt = os.path.join(OUTPUT_DIR, f"{name}.pdbqt")
        prepare_ligand(lig_pdb, lig_pdbqt)
        ligand_pdbqts.append((lig_pdbqt, name))

    gpf_file = os.path.join(OUTPUT_DIR, "grid.gpf")
    generate_gpf(
        ligand_pdbqt=ligand_pdbqts[0][0],
        receptor_pdbqt=receptor_pdbqt,
        output_gpf=gpf_file,
        center=GRID_CENTER,
        size=GRID_SIZE,
    )

    run_autogrid(gpf_file)
    fld_file = os.path.join(OUTPUT_DIR, "receptor.maps.fld")

    for lig_pdbqt, name in ligand_pdbqts:
        out_base = os.path.join(OUTPUT_DIR, f"{name}_docked")
        run_docking(lig_pdbqt, fld_file, out_base)

if __name__ == "__main__":
    main()
