import os
import sys
import subprocess

MGLTOOLS = os.path.expanduser("~/programs/mgltools_x86_64Linux2_1.5.7")
PYTHONSH = os.path.join(MGLTOOLS, "bin/pythonsh")
SCRIPT = os.path.join(MGLTOOLS, "MGLToolsPckgs/AutoDockTools/Utilities24/write_conformations_from_dlg.py")
def extract_all_poses(dlg_file):
    if not os.path.exists(dlg_file):
        print(f"[ERROR] DLG file not found: {dlg_file}")
        return

    dlg_dir = os.path.dirname(os.path.abspath(dlg_file))
    dlg_base = os.path.splitext(os.path.basename(dlg_file))[0]

    output_stem = os.path.join(dlg_dir, dlg_base + "_pose")

    cmd = f"{PYTHONSH} {SCRIPT} -d {dlg_file} -o {output_stem}"
    print(f"[RUN] {cmd}")
    subprocess.run(cmd, shell=True, check=True)

    print(f"[INFO] Extracted poses to: {output_stem}_*.pdbqt")
"""
def extract_all_poses(dlg_file, output_stem="ligand_docked_pose"):
    if not os.path.exists(dlg_file):
        print(f"[ERROR] DLG file not found: {dlg_file}")
        return

    dlg_dir = os.path.dirname(os.path.abspath(dlg_file))
    dlg_base = os.path.splitext(os.path.basename(dlg_file))[0]

    # If output_stem not passed, use the same base name as dlg
    if output_stem == "ligand_docked_pose":
        output_stem = os.path.join(dlg_dir, dlg_base + "_pose")
    else:
        output_stem = os.path.join(dlg_dir, output_stem)

    cmd = f"{PYTHONSH} {SCRIPT} -d {dlg_file} -o {output_stem}"
    print(f"[RUN] {cmd}")
    subprocess.run(cmd, shell=True, check=True)

    print(f"[INFO] Docked poses written to: {output_stem}_*.pdbqt")
"""
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python extract_all_docked_poses.py path/to/ligand_docked.dlg")
        sys.exit(1)

    dlg_file = sys.argv[1]
    extract_all_poses(dlg_file)
