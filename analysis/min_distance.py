import os
from Bio.PDB import PDBParser
import numpy as np
from scipy.spatial.distance import pdist

p = PDBParser(QUIET=True)

for i in range(1, 11):
    fname = f"receptor_{i:04d}.pdb"
    if not os.path.exists(fname):
        print(f"[{fname}] ❌ File not found.")
        continue

    try:
        structure = p.get_structure("receptor", fname)
        coords = np.array([atom.coord for atom in structure.get_atoms()])
        if len(coords) < 2:
            print(f"[{fname}] ⚠️ Not enough atoms to compute distance.")
            continue

        min_dist = pdist(coords).min()
        print(f"[{fname}] ✅ Min distance: {min_dist:.3f} Å")

    except Exception as e:
        print(f"[{fname}] ❌ Error: {e}")

