from pymol import cmd
import os

def main(kinase="WEE1"):
    base_dir = kinase
    receptor_path = os.path.join(base_dir, "receptor", "receptor.pdbqt")
    ligand_ref_path = os.path.join(base_dir, "ligands", "ligand.pdb")
    docking_dir = os.path.join(base_dir, "docking")

    # Load receptor (optional)
    if os.path.exists(receptor_path):
        cmd.load(receptor_path, "receptor")
        cmd.hide("everything", "receptor")
        cmd.show("cartoon", "receptor")
        cmd.color("cyan", "receptor")
        print(f"[INFO] Loaded receptor: {receptor_path}")
    else:
        print(f"[WARNING] No receptor found: {receptor_path}")

    # Load reference ligand
    cmd.load(ligand_ref_path, "ligand_ref")
    cmd.hide("everything", "ligand_ref")
    cmd.show("sticks", "ligand_ref")
    cmd.color("green", "ligand_ref")
    print(f"[INFO] Loaded reference ligand: {ligand_ref_path}")

    # Align all poses
    for i in range(1, 21):
        pose_name = f"ligand_docked_pose_{i}"
        pose_file = os.path.join(docking_dir, f"{pose_name}.pdbqt")

        if not os.path.exists(pose_file):
            print(f"[WARNING] Missing pose: {pose_file}")
            continue

        cmd.load(pose_file, pose_name)
        cmd.hide("everything", pose_name)
        cmd.show("sticks", pose_name)
        cmd.color("orange", pose_name)

        rmsd = cmd.align(pose_name, "ligand_ref")[0]
        print(f"[RMSD] {kinase} Pose {i:02d}: RMSD = {rmsd:.3f} Å")

    # Visualization
    cmd.bg_color("white")
    cmd.orient("ligand_ref or ligand_docked_pose_*")

    # Save image and session
    out_image = os.path.join(docking_dir, f"{kinase}_aligned.png")
    out_session = os.path.join(docking_dir, f"{kinase}_aligned.pse")
    cmd.png(out_image, dpi=300, ray=1)
    cmd.save(out_session)

    print(f"[DONE] Saved: {out_image}")
    print(f"[DONE] Saved: {out_session}")

# Optional entry point for headless use
if __name__ == "__main__":
    import sys
    kinase = sys.argv[1] if len(sys.argv) > 1 else "WEE1"
    main(kinase)

