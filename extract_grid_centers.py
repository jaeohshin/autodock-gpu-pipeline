#!/usr/bin/env python3
"""
Extracts grid box centers for each aligned kinase conformation
based on the cognate ligand in the reference structure.

Usage:
    pymol -cq extract_grid_centers.py -- receptors/kinase_001.pdb references/ref_with_ligand.pdb LIG centers/kinase_001.txt
"""

import sys
import os
from pymol import cmd

# Input arguments
receptor_pdb   = sys.argv[1]  # ABL1.renum.pdb
reference_pdb  = sys.argv[2]  # reference_with_ligand.pdb
ligand_resn    = sys.argv[3]  # JIN
output_file    = sys.argv[4]  # centers/ABL1.txt

# Load receptor (target of docking)
cmd.load(receptor_pdb, "receptor")

# Load reference structure with ligand
cmd.load(reference_pdb, "reference")

# Align reference to receptor (i.e., transform the ligand!)
cmd.align("reference and polymer", "receptor and polymer")

# Extract transformed ligand coordinates from the aligned reference
ligand_sel = f"reference and resn {ligand_resn} and chain A"
model = cmd.get_model(ligand_sel)

if not model.atom:
    print(f"[ERROR] Ligand {ligand_resn} in chain A not found in reference.")
    sys.exit(1)

coords = [atom.coord for atom in model.atom]
center = [sum(c[i] for c in coords) / len(coords) for i in range(3)]

# Ensure directory exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Save center
with open(output_file, "w") as f:
    f.write(f"{center[0]:.3f},{center[1]:.3f},{center[2]:.3f}\n")
