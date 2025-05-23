#!/usr/bin/env python3
"""
Extracts grid box centers for each aligned kinase conformation
based on the cognate ligand in the reference structure.

Usage:
    pymol -cq extract_grid_centers.py -- receptors/kinase_001.pdb references/ref_with_ligand.pdb LIG centers/kinase_001.txt
"""

import sys
from pymol import cmd
import os

# === Inputs ===
receptor_pdb = sys.argv[1]       # e.g., ABL1.renum.pdb
reference_pdb = sys.argv[2]      # e.g., reference_with_ligand.pdb
ligand_resn = sys.argv[3]        # e.g., JIN
output_file = sys.argv[4]        # e.g., centers/ABL1.txt

# === Load structures ===
cmd.load(receptor_pdb, "receptor")
cmd.load(reference_pdb, "reference")

# Extract ligand from reference and move it to new object so it follows alignment
ligand_sel = f"reference and resn {ligand_resn} and chain A"
cmd.create("ligand", ligand_sel)

# Align receptor to reference (transform is applied to "ligand" as well)
cmd.align("receptor and polymer", "reference and polymer")

# Now ligand is in receptor's frame — extract center
model = cmd.get_model("ligand")
if not model.atom:
    print(f"[ERROR] Ligand {ligand_resn} in chain A not found.")
    sys.exit(1)

coords = [atom.coord for atom in model.atom]
center = [sum(c[i] for c in coords) / len(coords) for i in range(3)]

# Output directory
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Save center
with open(output_file, "w") as f:
    f.write(f"{center[0]:.3f},{center[1]:.3f},{center[2]:.3f}\n")
