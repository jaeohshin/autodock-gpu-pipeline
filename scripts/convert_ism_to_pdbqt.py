import os
import sys
import gzip
import shutil
from multiprocessing import Pool, cpu_count
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

ROOT = "/store/jaeohshin/work/dock/dude_raw"
OUT_ROOT = "/store/jaeohshin/work/dock/virtual_screening/input/ligands_ism_pdbqt"

print("[INFO] Starting parallel Meeko-based ISM to PDBQT conversion...")

class CleanStderr:
    def write(self, msg):
        if ("get_symmetries_for_rmsd" not in msg and "too symmetric?" not in msg):
            sys.__stderr__.write(msg)
    def flush(self): pass

sys.stderr = CleanStderr()

def convert_ligand(args):
    i, smiles, chembl_id, out_dir, type_ = args
    newname = f"{type_}_{i:05d}_REMARKName={chembl_id}.pdbqt"
    out_path = os.path.join(out_dir, newname)

    if os.path.exists(out_path):
        return newname, None  # Already exists

    try:
        # RDKit build & 3-D embed
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("SMILES parsing failed")
        Chem.SanitizeMol(mol)
        mol.SetProp("_Name", chembl_id)
        mol = Chem.AddHs(mol) 
        if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
            raise ValueError("3D embedding failed")

        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol)
        else:
            AllChem.UFFOptimizeMolecule(mol)

        # RDKit debug
        
        print(f"[DEBUG] RDKit atoms for {chembl_id}:")
        for atom in mol.GetAtoms():
            print(f"  {atom.GetIdx():>2}: {atom.GetSymbol():<2} | "
                  f"Hybridization: {atom.GetHybridization()} | "
                  f"FormalCharge: {atom.GetFormalCharge()}")
            
        # Meeko preparation
        prep = MoleculePreparation()
        prep.add_hydrogens = False
        prep.ignore_protonation_states = False
        prep.generate_tautomers = False
        prep.rigid_macrocycles = True
        
        molsetups = prep.prepare(mol)  # <-- ADD THIS LINE!


        if not molsetups:
            raise ValueError("Meeko failed to generate molecule setup")
        molsetup = molsetups[0]
        

        writer = PDBQTWriterLegacy()
        pdbqt_str = writer.write_string(molsetup)[0]

        with open(out_path, "w") as out_f:
            out_f.write(pdbqt_str)

        return newname, None  # Success
    except Exception as e:
        return chembl_id, str(e)  # Error

# === Main loop ===
for target in sorted(os.listdir(ROOT)):
    target_path = os.path.join(ROOT, target)
    if not os.path.isdir(target_path):
        continue
    print(f"[TARGET] {target}")

    for type_ in ["actives", "decoys"]:
        ism_path = os.path.join(target_path, f"{type_}_final.ism")
        out_dir = os.path.join(OUT_ROOT, target, type_)
        list_out = os.path.join(OUT_ROOT, target, f"ligands_{type_}.list")
        error_log = os.path.join(OUT_ROOT, target, f"ligands_{type_}_errors.log")
        os.makedirs(out_dir, exist_ok=True)

        if not os.path.isfile(ism_path):
            print(f"  ⚠️ Missing ISM file: {ism_path}")
            continue

        jobs = []
        with open(ism_path) as f_in:
            for i, line in enumerate(f_in, 1):
                parts = line.strip().split()
                if len(parts) < 2:
                    continue
                smiles = parts[0]
                raw_id = parts[-1]
                chembl_id = raw_id if raw_id.startswith("CHEMBL") else f"CHEMBL{raw_id}"
                jobs.append((i, smiles, chembl_id, out_dir, type_))

        print(f"  → Converting {len(jobs)} ligands in parallel...")

        with Pool(processes=cpu_count()) as pool, \
             open(list_out, "w") as list_f, \
             open(error_log, "w") as err_f:
            for result_name, error in pool.imap_unordered(convert_ligand, jobs):
                if error:
                    err_f.write(f"{result_name}\t{error}\n")
                    print(f"    [ERROR] {result_name}: {error}")
                else:
                    list_f.write(f"{result_name}\n")
                    print(f"    ✓ {result_name}")

print("[DONE] Meeko-based ISM to PDBQT conversion complete.")

