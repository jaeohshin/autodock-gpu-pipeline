"""
/usr/bin/python3 -m pymol -cq prepare_grids.py --kinase abl1

"""
import os
import numpy as np
import requests
import sys
from pymol import cmd

INPUT_DIR = "../vs_crystal/input"
REF_LIST = os.path.join(INPUT_DIR, "kinase.txt")
RECEPTOR_DIR = os.path.join(INPUT_DIR, "receptors")
REFERENCE_PDB_DIR = os.path.join(INPUT_DIR, "references")
GRID_OUT_DIR = "../vs_crystal/preprocessed/grid_centers"

def download_pdb(pdb_id, save_path):
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        with open(save_path, "w") as f:
            f.write(response.text)
        print(f"[INFO] Downloaded {pdb_id} to {save_path}")
    else:
        raise FileNotFoundError(f"[ERROR] Failed to download PDB {pdb_id} from RCSB")

def auto_detect_ligand(obj_name):
    hetatms = cmd.get_model(f"{obj_name} and not polymer and not solvent")
    ligands = {}
    for atom in hetatms.atom:
        key = (atom.resn, atom.chain)
        ligands.setdefault(key, []).append(atom.coord)

    if not ligands:
        raise ValueError("No ligand found in reference structure")

    best_lig = max(ligands.items(), key=lambda x: len(x[1]))
    resn, chain = best_lig[0]
    coords = best_lig[1]
    center = np.mean(coords, axis=0)
    return resn, chain, center

def process_kinase(kinase, pdb_id):
    kinase = kinase.lower()
    receptor_path = os.path.join(RECEPTOR_DIR, kinase)
    reference_pdb = os.path.join(REFERENCE_PDB_DIR, f"{pdb_id}.pdb")
    grid_out_path = os.path.join(GRID_OUT_DIR, kinase)
    os.makedirs(grid_out_path, exist_ok=True)

    if not os.path.exists(reference_pdb):
        os.makedirs(REFERENCE_PDB_DIR, exist_ok=True)
        download_pdb(pdb_id, reference_pdb)

    for fname in sorted(os.listdir(receptor_path)):
        if not fname.endswith(".pdb"):
            continue

        idx = os.path.splitext(fname)[0].split("_")[-1]
        receptor_file = os.path.join(receptor_path, fname)
        center_outfile = os.path.join(grid_out_path, f"receptor_{idx}.txt")
        complex_outfile = os.path.join(grid_out_path, f"receptor_{idx}_complex.pdb")

        if os.path.exists(center_outfile) and os.path.exists(complex_outfile):
            print(f"[SKIP] {center_outfile} exists")
            continue

        cmd.reinitialize()
        cmd.load(reference_pdb, "ref")
        cmd.load(receptor_file, "mob")

        cmd.align("ref and polymer and chain A", "mob and polymer")

        try:
            resn, chain, ligand_center = auto_detect_ligand("ref")
        except ValueError as e:
            print(f"[WARN] {kinase} receptor_{idx}: {e}")
            continue

        with open(center_outfile, "w") as f:
            f.write("{:.3f} {:.3f} {:.3f}\n".format(*ligand_center))

        cmd.create("lig", f"ref and resn {resn} and chain {chain}")
        cmd.create("rec", "mob")
        cmd.save(complex_outfile, "rec or lig")

        print(f"[OK] {kinase} receptor_{idx} → {center_outfile}, {complex_outfile}")

def main():
    target_kinase = None
    if len(sys.argv) > 2 and sys.argv[1] == "--kinase":
        target_kinase = sys.argv[2].lower()

    with open(REF_LIST) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            kinase, pdb_id = parts
            if target_kinase and kinase.lower() != target_kinase:
                continue
            process_kinase(kinase, pdb_id)

main()
cmd.quit()
