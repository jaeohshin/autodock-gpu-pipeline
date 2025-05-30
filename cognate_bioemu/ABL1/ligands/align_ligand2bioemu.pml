bg_color white
load ../receptor/original.pdb, ref

load ../receptor/receptor_0001.pdb, mob
load ligand.pdb, lig
align ref and polymer and chain A, mob and polymer
matrix_copy ref, lig
save ligand_aligned_0001.pdb, lig
delete mob
delete lig

load ../receptor/receptor_0002.pdb, mob
load ligand.pdb, lig
align ref and polymer and chain A, mob and polymer
matrix_copy ref, lig
save ligand_aligned_0002.pdb, lig
delete mob
delete lig

load ../receptor/receptor_0003.pdb, mob
load ligand.pdb, lig
align ref and polymer and chain A, mob and polymer
matrix_copy ref, lig
save ligand_aligned_0003.pdb, lig
delete mob
delete lig

load ../receptor/receptor_0004.pdb, mob
load ligand.pdb, lig
align ref and polymer and chain A, mob and polymer
matrix_copy ref, lig
save ligand_aligned_0004.pdb, lig
delete mob
delete lig

load ../receptor/receptor_0005.pdb, mob
load ligand.pdb, lig
align ref and polymer and chain A, mob and polymer
matrix_copy ref, lig
save ligand_aligned_0005.pdb, lig
delete mob
delete lig

load ../receptor/receptor_0006.pdb, mob
load ligand.pdb, lig
align ref and polymer and chain A, mob and polymer
matrix_copy ref, lig
save ligand_aligned_0006.pdb, lig
delete mob
delete lig

load ../receptor/receptor_0007.pdb, mob
load ligand.pdb, lig
align ref and polymer and chain A, mob and polymer
matrix_copy ref, lig
save ligand_aligned_0007.pdb, lig
delete mob
delete lig

load ../receptor/receptor_0008.pdb, mob
load ligand.pdb, lig
align ref and polymer and chain A, mob and polymer
matrix_copy ref, lig
save ligand_aligned_0008.pdb, lig
delete mob
delete lig

load ../receptor/receptor_0009.pdb, mob
load ligand.pdb, lig
align ref and polymer and chain A, mob and polymer
matrix_copy ref, lig
save ligand_aligned_0009.pdb, lig
delete mob
delete lig

load ../receptor/receptor_0010.pdb, mob
load ligand.pdb, lig
align ref and polymer and chain A, mob and polymer
matrix_copy ref, lig
save ligand_aligned_0010.pdb, lig
delete mob
delete lig
quit
