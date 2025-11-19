import os
import sys
from pymol import cmd

def main():
    # === Parse kinase name ===
    if len(sys.argv) != 2 or not sys.argv[1].startswith("--"):
        print("Usage: pymol -cq align_all_docked_poses.py --KINASE_NAME")
        sys.exit(1)

    kinase = sys.argv[1][2:]
    print(f"[INFO] Processing kinase: {kinase}")

    # === Define paths ===
    ligand_ref_path = os.path.join(kinase, "ligands", "ligand.pdb")
    docked_dir = os.path.join(kinase, "docking")
    receptor_path = os.path.join(kinase, "receptor", "receptor.pdb")
    rmsd_log_path = os.path.join(docked_dir, "rmsd_report.txt")
    output_img = os.path.join(docked_dir, f"{kinase}_aligned.png")

    # === Check reference ligand ===
    if not os.path.exists(ligand_ref_path):
        print(f"[ERROR] ligand.pdb not found: {ligand_ref_path}")
        sys.exit(1)

    print(f"[DEBUG] Loading reference ligand from: {ligand_ref_path}")
    cmd.load(ligand_ref_path, "ref_ligand")
    cmd.hide("everything", "ref_ligand")
    cmd.show("sticks", "ref_ligand")
    cmd.color("green", "ref_ligand")

    if not cmd.count_atoms("ref_ligand"):
        print("[ERROR] Reference ligand failed to load or is empty.")
        sys.exit(1)

    # === Optional: Load receptor ===
    if os.path.exists(receptor_path):
        print(f"[DEBUG] Loading receptor from: {receptor_path}")
        cmd.load(receptor_path, "receptor")
        cmd.hide("everything", "receptor")
        cmd.show("cartoon", "receptor")
        cmd.color("cyan", "receptor")
    else:
        print(f"[WARNING] No receptor.pdb found at: {receptor_path}")

    # === Start logging RMSD ===
    print(f"[INFO] Writing RMSDs to: {rmsd_log_path}")
    with open(rmsd_log_path, "w") as log:
        log.write("Pose\tRMSD (Å)\n")

        for i in range(1, 21):
            pose_name = f"ligand_docked_pose_{i}"
            pose_path = os.path.join(docked_dir, f"{pose_name}.pdbqt")

            if not os.path.exists(pose_path):
                print(f"[WARNING] Missing pose file: {pose_path}")
                continue

            print(f"[DEBUG] Loading: {pose_path}")
            cmd.load(pose_path, pose_name)
            cmd.hide("everything", pose_name)
            cmd.show("sticks", pose_name)
            cmd.color("orange", pose_name)

            result = cmd.align(pose_name, "ref_ligand")
            rmsd = result[0]
            print(f"[RESULT] {pose_name}: RMSD = {rmsd:.3f} Å")
            log.write(f"{pose_name}\t{rmsd:.3f}\n")

    # === Finish up ===
    cmd.bg_color("white")
    cmd.orient("ref_ligand or ligand_docked_pose_*")

    try:
        cmd.png(output_img, dpi=300, ray=1)
        print(f"[INFO] Saved image to: {output_img}")
    except Exception as e:
        print(f"[WARNING] Could not save image: {e}")

    print("[DONE] PyMOL alignment script completed.")

if __name__ == "__main__":
    if len(sys.argv) == 2:
        main()
    else:
        # Allow calling main() directly from inside PyMOL
        main_arg = "--ABL1"  # change this to your default test target
        sys.argv = ["align_all_docked_poses.py", main_arg]
        main()
