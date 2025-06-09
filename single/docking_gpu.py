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
SPLIT_ALTLOCS    = f"{MGL_HOME}/bin/pythonsh {MGL_HOME}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_pdb_split_alt_confs.py"

AUTOGRID_BIN     = "autogrid4"

# 👇 Update this to your actual binding site center
GRID_CENTER = (25, 25.0, 20.0)
GRID_SIZE   = (40, 40, 40)

# === UTILITY ===
def run_cmd(cmd):
    print(f"[RUN] {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def has_altlocs(pdb_path):
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and line[16] not in (" ", "A"):  # altLoc field (column 17)
                return True
    return False  

def remove_waters(pdb_in, pdb_out):
    with open(pdb_in, 'r') as fin, open(pdb_out, 'w') as fout:
        for line in fin:
            if not (line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == "HOH"):
                fout.write(line)


# === STEP 1: Prepare receptor ===
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
    run_cmd(f"{PREPARE_RECEPTOR} -r {pdb_file_abs} -o {output_pdbqt_abs} -A checkhydrogens")

# -U nphs,lps,waters -e 

# === STEP 2: Prepare ligand ===

def prepare_ligand(pdb_file, output_pdbqt):
    pdb_file_abs = os.path.abspath(pdb_file)
    output_pdbqt_abs = os.path.abspath(output_pdbqt)
    lig_dir = os.path.dirname(pdb_file_abs)
    lig_name = os.path.basename(pdb_file_abs)
    
    cmd = f"cd {lig_dir} && {PREPARE_LIGAND} -l {lig_name} -o {output_pdbqt_abs} -A hydrogens"
    run_cmd(cmd)

"""

def prepare_ligand(pdb_file, output_pdbqt, use_meeko=True):
    import os
    import shutil
    import tempfile
    from rdkit import Chem
    from meeko import MoleculePreparation

    def try_meeko(mol, label=""):
        try:
            prep = MoleculePreparation()
            prep.automatic_detection = True
            prep.add_polar_hydrogens = True
            prep.assign_charges = True
            prep.uff_energy_minimize = True
            pdbqt_string = prep.prepare(mol)
            with open(output_pdbqt, "w") as f:
                f.write(pdbqt_string)
            print(f"[INFO] Ligand prepared using Meeko ({label}): {output_pdbqt}")
            return True
        except Exception as e:
            print(f"[WARN] Meeko failed ({label}): {e}")
            return False

    if use_meeko:
        # === Attempt 1: Raw PDB ===
        mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
        if mol is not None and try_meeko(mol, "raw PDB"):
            return

        # === Attempt 2: .mol2 preprocessing ===
        print(f"[INFO] Attempting Open Babel preprocessing to .mol2 for Meeko...")
        tmpdir = tempfile.mkdtemp()
        mol2_path = os.path.join(tmpdir, "ligand_fixed.mol2")
        os.system(f"obabel {pdb_file} -O {mol2_path} --gen3D --addh")

        if os.path.exists(mol2_path):
            mol2_mol = Chem.MolFromMol2File(mol2_path, sanitize=True, removeHs=False)
            if mol2_mol is not None and try_meeko(mol2_mol, ".mol2"):
                shutil.rmtree(tmpdir)
                return

            # === Attempt 3: Full Open Babel → PDBQT fallback ===
            print("[WARN] RDKit failed on MOL2. Using Open Babel to convert directly to PDBQT.")
            os.system(f"obabel {pdb_file} -O {output_pdbqt} --gen3D --addh --partialcharge gasteiger")
            shutil.rmtree(tmpdir)
            return

        shutil.rmtree(tmpdir)

    # === Final fallback: MGLTools ===
    print("[WARN] Meeko failed completely. Falling back to MGLTools...")
    pdb_file_abs = os.path.abspath(pdb_file)
    output_pdbqt_abs = os.path.abspath(output_pdbqt)
    lig_dir = os.path.dirname(pdb_file_abs)
    lig_name = os.path.basename(pdb_file_abs)
    cmd = f"cd {lig_dir} && {PREPARE_LIGAND} -l {lig_name} -o {output_pdbqt_abs}"
    run_cmd(cmd)

"""

# === STEP 3: Generate GPF ===
def generate_gpf(ligand_pdbqt, receptor_pdbqt, output_gpf, center, size):
    output_dir = os.path.dirname(os.path.abspath(output_gpf))
    ligand_name = os.path.basename(ligand_pdbqt)
    receptor_name = os.path.basename(receptor_pdbqt)
    gpf_name = os.path.basename(output_gpf)


    # 🔥 Force remove old GPF if it exists
    if os.path.exists(output_gpf):
        print(f"[INFO] Removing old GPF: {output_gpf}")
        os.remove(output_gpf)
            
    center_str = f"{center[0]},{center[1]},{center[2]}"
    size_str = f"{size[0]},{size[1]},{size[2]}"
    
    cmd = (
        f"cd {output_dir} && {PREPARE_GPF} "
        f"-l {ligand_name} -r {receptor_name} -o {gpf_name} "
        f"-p npts={size_str} -p gridcenter={center_str} "
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
    lig_pdbqt_abs = os.path.abspath(lig_pdbqt)
    fld_file_abs = os.path.abspath(fld_file)
    output_base_abs = os.path.abspath(output_basename)
    
    cmd = (
        f"{AUTODOCK_GPU_BIN} "
        f"--ffile {fld_file_abs} "
        f"--lfile {lig_pdbqt_abs} "
        f"--nrun 50 "
        f"--nev 2500000 " 
        f"--ngen 42000 " 
        f"--heuristics 1 "
        f"--autostop 1 "
        f"--lsrat 50 "
        f"--resnam {output_base_abs} "    
    )
    run_cmd(cmd)

# === MAIN PIPELINE ===
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    receptor_pdbqt = os.path.join(OUTPUT_DIR, "receptor.pdbqt")
    prepare_receptor(RECEPTOR_PDB, receptor_pdbqt)

    """
    ligand_pdbqts = []
    for lig_pdbqt in sorted(glob(os.path.join(LIGAND_PDB_DIR, "*.pdbqt"))):
        name = os.path.splitext(os.path.basename(lig_pdbqt))[0]
        ligand_pdbqts.append((lig_pdbqt, name))

    """
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
