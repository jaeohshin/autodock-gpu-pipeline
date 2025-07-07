#!/usr/bin/env python3

import sys
import warnings
import os
import gzip
import shutil
from rdkit import Chem
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy


ROOT = "/data/work/dock_local/dude_raw"
OUT_ROOT = "/data/work/dock_local/virtual_screening/input/ligands_sdf2meeko"

print("[INFO] Starting Meeko-based SDF to PDBQT conversion...")

class CleanStderr:
    def write(self, msg):
        if ("get_symmetries_for_rmsd" not in msg and 
            "too symmetric?" not in msg):
            sys.__stderr__.write(msg)
    def flush(self):
        pass

sys.stderr = CleanStderr()



for target in sorted(os.listdir(ROOT)):
    target_path = os.path.join(ROOT, target)
    if not os.path.isdir(target_path):
        continue
    print(f"[TARGET] {target}")

    for type_ in ["actives", "decoys"]:
        sdf_gz = os.path.join(target_path, f"{type_}_final.sdf.gz")
        sdf = sdf_gz[:-3]
        out_dir = os.path.join(OUT_ROOT, target, type_)
        list_out = os.path.join(OUT_ROOT, target, f"ligands_{type_}.list")
        error_log = os.path.join(OUT_ROOT, target, f"ligands_{type_}_errors.log")
        os.makedirs(out_dir, exist_ok=True)

        # Decompress if needed
        if os.path.isfile(sdf_gz):
            print(f"  → Decompressing {sdf_gz}")
            with gzip.open(sdf_gz, 'rb') as f_in:
                with open(sdf, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

        if not os.path.isfile(sdf):
            print(f"  ⚠️ Missing SDF file: {sdf}")
            continue

        print(f"  → Reading and converting ligands from {sdf}")
        suppl = Chem.SDMolSupplier(sdf, removeHs=False)

        prep = MoleculePreparation()
        prep.add_hydrogens = False
        prep.ignore_protonation_states = True
        prep.generate_tautomers = False
        prep.rigid_macrocycles = True

        with open(list_out, 'w') as list_f, open(error_log, 'w') as err_f:
            for i, mol in enumerate(suppl, start=1):
                if mol is None:
                    continue

                chembl_id = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else f"ligand_{i:05d}"
                mol.SetProp("_Name", chembl_id)
                newname = f"{type_}_{i:05d}_REMARKName={chembl_id}.pdbqt"
                out_path = os.path.join(out_dir, newname)

                try:
                    # Recast to strip query features and avoid HasQuery errors
                    mol_clean = Chem.Mol(mol)
                    mol_setups = prep.prepare(mol_clean)
                    writer = PDBQTWriterLegacy()
                    pdbqt_string = writer.write_string(mol_setups[0])[0]
                    with open(out_path, "w") as f:
                        f.write(pdbqt_string)
                    

                    list_f.write(f"{newname}\n")
                except Exception as e:
                    print(f"[ERROR] {chembl_id}: {e}")
                    err_f.write(f"{chembl_id}\t{e}\n")

print("[DONE] Meeko-based conversion finished.")
