#!/usr/bin/env python3
"""
dlg2pdb_all_poses.py

Extracts all ligand poses from a DLG file into a multi-model PDB file.

Usage:
    python dlg2pdb_all_poses.py docking.dlg output_all_poses.pdb
"""

import sys

if len(sys.argv) != 3:
    print("Usage: python dlg2pdb_all_poses.py docking.dlg output_all_poses.pdb")
    sys.exit(1)

dlg_file = sys.argv[1]
output_pdb = sys.argv[2]

poses = []
current_pose = []
inside_model = False

with open(dlg_file, 'r') as f:
    for line in f:
        if line.startswith("DOCKED: MODEL"):
            inside_model = True
            current_pose = ["MODEL\n"]
        elif line.startswith("DOCKED: ENDMDL") and inside_model:
            current_pose.append("ENDMDL\n")
            poses.append(current_pose)
            inside_model = False
        elif inside_model and (line.startswith("DOCKED: ATOM") or line.startswith("DOCKED: HETATM")):
            # Strip "DOCKED: " and keep the rest
            pdb_line = line.replace("DOCKED: ", "", 1)
            current_pose.append(pdb_line)

if not poses:
    print("[ERROR] No valid ligand poses found in the DLG file.")
    sys.exit(1)

with open(output_pdb, 'w') as f:
    for pose in poses:
        f.writelines(pose)

print(f"[INFO] Extracted {len(poses)} poses to {output_pdb}")
