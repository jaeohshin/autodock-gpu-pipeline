import os
import pymol2
import numpy as np
import requests

INPUT_DIR = "../input"
REF_LIST = os.path.join(INPUT_DIR, "reference.list")
RECEPTOR_DIR = os.path.join(INPUT_DIR, "receptors")
REFERENCE_PDB_DIR = os.path.join(INPUT_DIR, "references")  # expected to contain downloaded PDBs
GRID_OUT_DIR = "../preprocessed/grid_centers"


def download_pdb(pdb_id, save_path):
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        with open(save_path, "w") as f:
            f.write(response.text)
        print(f"[INFO] Downloaded {pdb_id} to {save_path}")
    else:
        raise FileNotFoundError(f"[ERROR] Failed to download PDB {pdb_id} from RCSB")


def auto_detect_ligand(cmd, obj_name):
    hetatms = cmd.get_model(f"{obj_name} and not polymer and not solvent")
    ligands = {}
    for atom in hetatms.atom:
        key = (atom.resn, atom.chain)
        ligands.setdefault(key, []).append(atom.coord)

    if not ligands:
        raise ValueError("No ligand found in reference structure")

    # Choose the largest ligand group (by number of atoms)
    best_lig = max(ligands.items(), key=lambda x: len(x[1]))
    resn, chain = best_lig[0]
    coords = best_lig[1]
    center = np.mean(coords, axis=0)
    return resn, chain, center


def process_kinase(cmd, kinase, pdb_id):
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

        # Align receptor chains
        cmd.align("ref and polymer and chain A", "mob and polymer")

        # Auto-detect ligand and compute center
        try:
            resn, chain, ligand_center = auto_detect_ligand(cmd, "ref")
        except ValueError as e:
            print(f"[WARN] {kinase} receptor_{idx}: {e}")
            continue

        # Save grid center
        with open(center_outfile, "w") as f:
            f.write("{:.3f} {:.3f} {:.3f}\n".format(*ligand_center))

        # Save complex (receptor + aligned ligand)
        cmd.create("lig", f"ref and resn {resn} and chain {chain}")
        cmd.create("rec", "mob")
        cmd.save(complex_outfile, "rec or lig")

        print(f"[OK] {kinase} receptor_{idx} → {center_outfile}, {complex_outfile}")


def main():
    with open(REF_LIST) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            kinase, pdb_id = parts
            process_kinase(pymol2.PyMOL(), kinase, pdb_id)


if __name__ == "__main__":
    main()
