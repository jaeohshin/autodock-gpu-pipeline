# Load input files
load reference_with_ligand.pdb, ref
load frame_re_100.pdb, frame
load ligand_poses.pdb, poses

# Visual setup
show cartoon, ref and polymer
show cartoon, frame
select ligand_crystal, ref and not polymer
show sticks, ligand_crystal
show sticks, poses
hide spheres, ligand_crystal

# Color settings
color green, ref and chain A
color red, ligand_crystal
color blue, frame
color orange, poses

# Align reference protein (ref) to model (frame)
# Align chain A only (polymer part)
align ref and polymer, frame and polymer

# Apply the alignment transform to entire ref structure (protein + ligand)
# (ref was transformed already by align — poses stay fixed)
# But for completeness, you can reapply matrix to ref explicitly:
matrix = get_object_matrix ref
transform_object ref, matrix

# Final visualization settings
hide everything, hydro
bg_color white
set ray_opaque_background, 0
