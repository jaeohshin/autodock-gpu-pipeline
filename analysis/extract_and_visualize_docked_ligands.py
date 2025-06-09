import os
import subprocess

KINASE = "ABL1"
ENSEMBLE_DIR = f"./{KINASE}/output"
NUM_STRUCTURES = 10
POSE_FILENAME = "ligand_docked_pose.pdb"
PML_FILENAME = f"{KINASE}_docked_visualization.pml"

def extract_top_pose(dlg_file, out_pdbqt):
    with open(dlg_file, "r") as f:
        lines = f.readlines()

    pose_lines = []
    in_pose = False
    for line in lines:
        if line.startswith("DOCKED:"):
            content = line[len("DOCKED:"):].lstrip()
            if content.startswith(("ROOT", "ENDROOT", "ATOM", "HETATM", "TORSDOF")):
                pose_lines.append(content)
                in_pose = True
        elif in_pose and line.startswith("DOCKED: ENDMDL"):
            break

    with open(out_pdbqt, "w") as f:
        f.writelines(pose_lines)
    print(f"[INFO] Extracted pose to: {out_pdbqt}")

def convert_pdbqt_to_pdb(pdbqt_file, pdb_file):
    subprocess.run(["obabel", pdbqt_file, "-O", pdb_file], check=True)
    print(f"[INFO] Converted to: {pdb_file}")

def generate_pymol_script(pml_path, receptor_paths, ligand_paths):
    with open(pml_path, "w") as f:
        f.write("bg_color white\n")
        for i, (rec, lig) in enumerate(zip(receptor_paths, ligand_paths), 1):
            f.write(f"load {rec}, receptor_{i:04d}\n")
            f.write(f"load {lig}, ligand_{i:04d}\n")
            f.write(f"show cartoon, receptor_{i:04d}\n")
            f.write(f"show sticks, ligand_{i:04d}\n")
            f.write(f"color cyan, receptor_{i:04d}\n")
            f.write(f"color red, ligand_{i:04d}\n\n")
        f.write("hide everything, hydro\n")
        f.write("set ray_opaque_background, off\n")
    print(f"[INFO] Wrote PyMOL script to: {pml_path}")

def main():
    receptor_paths = []
    ligand_paths = []

    for i in range(1, NUM_STRUCTURES + 1):
        idx = f"{i:04d}"
        pair_dir = os.path.join(ENSEMBLE_DIR, f"ensemble_{idx}")
        dlg_file = os.path.join(pair_dir, "ligand_docked.dlg")
        rec_pdbqt = os.path.join(pair_dir, "receptor.pdbqt")
        lig_pdbqt = os.path.join(pair_dir, "ligand_docked_pose.pdbqt")
        lig_pdb = os.path.join(pair_dir, POSE_FILENAME)

        if not os.path.exists(dlg_file):
            print(f"[WARN] Missing {dlg_file}, skipping {idx}.")
            continue

        extract_top_pose(dlg_file, lig_pdbqt)
        convert_pdbqt_to_pdb(lig_pdbqt, lig_pdb)

        receptor_paths.append(rec_pdbqt)
        ligand_paths.append(lig_pdb)

    generate_pymol_script(PML_FILENAME, receptor_paths, ligand_paths)

if __name__ == "__main__":
    main()

