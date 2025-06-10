reference_pdb = "../receptor/original.pdb"
ligand_pdb = "ligand.pdb"
output_script_path = "align_ligand2bioemu.pml"

script_lines = [
    "bg_color white\n",
    f"load {reference_pdb}, ref\n",
]

# Generate transformation blocks
for i in range(1, 51):  # adjust upper limit as needed
    receptor_pdb = f"../receptor/receptor_{i:04d}.pdb"
    aligned_ligand_pdb = f"ligand_aligned_{i:04d}.pdb"

    script_lines.extend([
        f"\nload {receptor_pdb}, mob\n",
        f"load {ligand_pdb}, lig\n",
        f"align ref and polymer and chain A, mob and polymer\n",
        f"matrix_copy ref, lig\n",
        f"save {aligned_ligand_pdb}, lig\n",
        f"delete mob\n",
        f"delete lig\n",  # 🛠️ ensures clean reload each loop
    ])

script_lines.append("quit\n")

with open(output_script_path, "w") as f:
    f.writelines(script_lines)

print(f"[INFO] PyMOL script written to: {output_script_path}")
