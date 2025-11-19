import os
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

def run_all_from_list(kinase_list_file):
    with open(kinase_list_file) as f:
        for line in f:
            if not line.strip():
                continue
            kinase = line.strip().split()[0]
            dlg_path = os.path.join(kinase, "docking", "ligand_docked.dlg")
            print(f"\n=== Processing {kinase} ===")
            extract_all_poses(dlg_path)

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python extract_all_kinase_poses.py kinase.txt")
        sys.exit(1)

    kinase_list_file = sys.argv[1]
    run_all_from_list(kinase_list_file)

