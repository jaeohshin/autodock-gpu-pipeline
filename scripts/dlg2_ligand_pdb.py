"""
dlg2_ligand_pdb.py

Author: Jaeoh Shin (KIAS)
Date: 2025/05/26

Purpose:
--------
Automatically extract the main ligand from crystal structure,
detect which chain it is in contact with, align the crystal complex 
to the generated receptor model, and export the ligand into the model's
coordinate frame for cognate docking.

Usage:
------
pymol -cq transform_ligand_auto_detect.py -- crystal_complex.pdb generated_model.pdb transformed_ligand.pdb

"""

import sys
from pymol import cmd

# === Input arguments ===
crystal_pdb   = sys.argv[1]
model_pdb     = sys.argv[2]
output_ligand = sys.argv[3]

# === Load structures ===
cmd.load(crystal_pdb, "crystal_full")
cmd.load(model_pdb, "model_full")

# === Fix missing chain ID in model if necessary ===
model = cmd.get_model("model_full")
if not any(atom.chain.strip() for atom in model.atom):
    print("[INFO] No chain ID found in model structure. Assigning chain A.")
    for atom in model.atom:
        atom.chain = "A"
    cmd.delete("model_full")
    cmd.load_model(model, "model_full", state=1)

# === Identify first polymer chain in crystal structure ===
chains = []
for atom in cmd.get_model("crystal_full and polymer").atom:
    if atom.chain not in chains:
        chains.append(atom.chain)

if not chains:
    print("[ERROR] No polymer chains found in crystal structure.")
    sys.exit(1)

first_chain = chains[0]
print(f"[INFO] First protein chain detected: {first_chain}")

# === Detect main ligand (largest HETATM group) ===
het_counts = {}
for atom in cmd.get_model("crystal_full and not polymer").atom:
    key = (atom.resn, atom.resi, atom.chain)
    het_counts[key] = het_counts.get(key, 0) + 1

if not het_counts:
    print("[ERROR] No ligand or HETATM entries found.")
    sys.exit(1)

ligand_resn, ligand_resi, ligand_chain = max(het_counts.items(), key=lambda x: x[1])[0]
print(f"[INFO] Auto-detected ligand: resn {ligand_resn}, resi {ligand_resi}, chain {ligand_chain}")

# === Create chain-only objects for alignment ===
cmd.create("crystal_chain", f"crystal_full and chain {first_chain} and polymer")
cmd.create("model_chain",   f"model_full and chain {first_chain} and polymer")

if not cmd.count_atoms("crystal_chain") or not cmd.count_atoms("model_chain"):
    print("[ERROR] One of the alignment chains is empty.")
    sys.exit(1)

# === Perform alignment and apply transform to full crystal object ===
alignment_result = cmd.align("crystal_chain", "model_chain")
if alignment_result is None:
    print("[ERROR] Alignment failed.")
    sys.exit(1)

print(f"[INFO] Alignment RMSD: {alignment_result[0]:.3f} Å")

# === Apply same transform to the entire crystal complex (incl. ligand) ===
matrix = cmd.get_object_matrix("crystal_chain")
if matrix is not None:
    cmd.transform_object("crystal_full", matrix)
else:
    print("[ERROR] Failed to retrieve transformation matrix.")
    sys.exit(1)

# === Select ligand atoms that are in contact with the protein ===
cmd.select("ligand_all", f"crystal_full and resn {ligand_resn} and resi {ligand_resi} and chain {ligand_chain}")
cmd.select("receptor_chain", f"crystal_full and chain {first_chain} and polymer")
cmd.select("ligand_in_contact", "ligand_all within 5 of receptor_chain")

lig_model = cmd.get_model("ligand_in_contact")
if not lig_model.atom:
    print("[ERROR] Ligand is not in contact with receptor chain after alignment.")
    sys.exit(1)

print(f"[INFO] Found {len(lig_model.atom)} ligand atoms in contact with chain {first_chain}")

# === Save output ===
cmd.save(output_ligand, "ligand_in_contact")
print(f"[INFO] Saved transformed ligand to: {output_ligand}")
