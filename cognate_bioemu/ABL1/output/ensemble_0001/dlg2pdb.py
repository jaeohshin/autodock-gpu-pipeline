from meeko.dlg_parser import DLGParser
from meeko import PDBQTMolecule
from openbabel import pybel

# === CONFIGURATION ===
dlg_path = "ligand_docked.dlg"  # change to your dlg path
output_pdb = "ligand_docked_best.pdb"

# === PARSE DLG ===
parser = DLGParser()
parser.parse(dlg_path)

# Get PDBQTMolecule of the first pose (best scoring one)
molecule = parser.molecule
mol = molecule.get_pybel_mol(conformer=0)

# Write to PDB
mol.write("pdb", output_pdb, overwrite=True)

print(f"[INFO] Best docked pose saved to: {output_pdb}")

