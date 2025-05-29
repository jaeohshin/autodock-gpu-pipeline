#!/usr/bin/env python3
"""
Extract grid center from receptor/ligand complex or via alignment to a reference.

Usage:
    pymol_sys -cq extract_grid_center.py -- receptor.pdb reference_with_ligand.pdb LIG output.txt

    For cognate docking:
    pymol_sys -cq extract_grid_center.py -- receptor_with_ligand.pdb cognate LIG output.txt
"""

import sys
import os
from pymol import cmd

# === Input arguments ===
if len(sys.argv) < 5:
    print("Usage:\n"
          "  For cognate: pymol_sys -cq extract_grid_center.py -- receptor_with_ligand.pdb cognate LIG output.txt\n"
          "  For non-cognate: pymol_sys -cq extract_grid_center.py -- receptor.pdb reference.pdb LIG output.txt")
    sys.exit(1)

receptor_pdb  = sys.argv[1]
mode_or_ref   = sys.argv[2]  # either 'cognate' or reference pdb
ligand_resn   = sys.argv[3]
output_file   = sys.argv[4]
ligand_chain  = "A"  # default chain

cmd.load(receptor_pdb, "receptor")

# === Cognate docking mode ===
if mode_or_ref.lower() == "cognate":
    ligand_sel = f"receptor and resn {ligand_resn} and chain {ligand_chain}"
    model = cmd.get_model(ligand_sel)
    
    if not model.atom:
        print(f"[ERROR] Ligand '{ligand_resn}' in chain '{ligand_chain}' not found in receptor.")
        sys.exit(1)

    print("[INFO] Running in cognate mode.")

        # Optional output of ligand PDB
    if len(sys.argv) > 5:
        ligand_outfile = sys.argv[5]
        cmd.select("ligand_selected", ligand_sel)
        cmd.save(ligand_outfile, "ligand_selected")
        print(f"[INFO] Ligand saved to: {ligand_outfile}")


else:
    # === Non-cognate: align reference with ligand ===
    reference_pdb = mode_or_ref
    cmd.load(reference_pdb, "reference")
    rmsd = cmd.align("reference and polymer", "receptor and polymer")[0]
    print(f"[INFO] Alignment RMSD: {rmsd:.3f} Å")

    ligand_sel = f"reference and resn {ligand_resn} and chain {ligand_chain}"
    model = cmd.get_model(ligand_sel)

    if not model.atom:
        print(f"[ERROR] Ligand '{ligand_resn}' in chain '{ligand_chain}' not found in reference.")
        sys.exit(1)

# === Compute grid center ===
coords = [atom.coord for atom in model.atom]
center = [sum(c[i] for c in coords) / len(coords) for i in range(3)]
center_str = f"{center[0]:.3f},{center[1]:.3f},{center[2]:.3f}"

# === Write to file ===
output_dir = os.path.dirname(output_file)
if output_dir:
    os.makedirs(output_dir, exist_ok=True)

with open(output_file, "w") as f:
    f.write(center_str + "\n")

print(f"[INFO] Grid center saved to {output_file}: {center_str}")
cmd.quit()
