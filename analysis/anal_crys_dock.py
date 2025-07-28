"""
Analyze the docking output file (.dlg) for a single crystal structure.

How to Run:

python script.py --kinase abl1
"""

import os
import re
import json
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

def extract_best_energy_from_dlg(filepath):
    energies = []
    with open(filepath, 'r') as f:
        for line in f:
            if "Estimated Free Energy of Binding" in line:
                match = re.search(r"=\s*(-?\d+\.\d+)", line)
                if match:
                    energies.append(float(match.group(1)))
    return min(energies) if energies else None

def label_ligand(filename):
    if filename.startswith("actives_"):
        return "active"
    elif filename.startswith("decoys_"):
        return "decoy"
    return "unknown"


def collect_best_scores(base_dir):
    from collections import defaultdict

    energy_dict = defaultdict(list)
    label_dict = {}

    if not os.path.isdir(base_dir):
        raise FileNotFoundError(f"[ERROR] Directory not found: {base_dir}")

    for fname in tqdm(os.listdir(base_dir), desc="Parsing ligands"):
        if not fname.endswith(".dlg"):
            continue
        match = re.search(r'REMARKName=(ZINC\d+|CHEMBL\d+)', fname)
        if not match:
            continue  # Skip malformed files
        ligand_id = match.group(1)  # Could be ZINC or CHEMBL
        fpath = os.path.join(base_dir, fname)
        energy = extract_best_energy_from_dlg(fpath)
        if energy is not None:
            energy_dict[ligand_id].append(energy)
            if ligand_id not in label_dict:
                label_dict[ligand_id] = label_ligand(fname)
    
    records = []
    for zinc_id, energies in energy_dict.items():
        best_energy = min(energies)
        label = label_dict.get(zinc_id, "unknown")
        records.append((zinc_id, best_energy, label))

    df = pd.DataFrame(records, columns=["Ligand", "Energy", "Label"])
    df["Best_Receptor"] = "crystal"
    return df



"""
def collect_best_scores(base_dir):
    records = []
    if not os.path.isdir(base_dir):
        raise FileNotFoundError(f"[ERROR] Directory not found: {base_dir}")

    for fname in tqdm(os.listdir(base_dir), desc="Parsing ligands"):
        if not fname.endswith(".dlg"):
            continue
        ligand = fname[:-4]
        fpath = os.path.join(base_dir, fname)
        energy = extract_best_energy_from_dlg(fpath)
        if energy is not None:
            records.append((ligand, energy, label_ligand(fname)))

    df = pd.DataFrame(records, columns=["Ligand", "Energy", "Label"])
    df["Best_Receptor"] = "crystal"  # fixed value
    return df
"""

def calculate_ef(df, top_percent):
    df_sorted = df.sort_values("Energy")
    total = len(df)
    actives_total = (df["Label"] == "active").sum()

    if total == 0 or actives_total == 0:
        print(f"[WARN] Not enough data to compute EF@{top_percent}%")
        return float("nan")

    top_n = max(1, int(total * top_percent / 100))
    actives_top = (df_sorted.head(top_n)["Label"] == "active").sum()
    ef = (actives_top / top_n) / (actives_total / total)
    return ef

def plot_score_histogram(df, out_path):
    if df.empty:
        print("[WARN] DataFrame is empty. Skipping histogram.")
        return

    if df["Energy"].nunique() < 2:
        print("[WARN] Not enough unique Energy values to plot KDE. Plotting histogram only.")
        kde_flag = False
    else:
        kde_flag = True

    plt.figure(figsize=(8, 4))
    sns.histplot(data=df, x="Energy", hue="Label", bins=50, kde=kde_flag)
    plt.title("Binding Energy Distribution")
    plt.xlabel("Best Binding Energy (kcal/mol)")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"[INFO] Histogram saved to {out_path}")

def per_receptor_ef(df, top_percent=1):
    top_n = max(1, int(len(df) * top_percent / 100))
    top_df = df.sort_values("Energy").head(top_n)
    return top_df[top_df["Label"] == "active"]["Best_Receptor"].value_counts()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--kinase", required=True, help="Kinase target name (e.g., abl1)")
    args = parser.parse_args()

    kinase = args.kinase.lower()
    docking_dir = f"../vs_crystal/docking_output/{kinase}"
    output_dir = f"../vs_crystal/results/{kinase}"
    os.makedirs(output_dir, exist_ok=True)

    df = collect_best_scores(docking_dir)
    print(f"[INFO] {kinase.upper()}  — Ligands: {len(df)}, Actives: {(df['Label']=='active').sum()}, Decoys: {(df['Label']=='decoy').sum()}")

    df.to_csv(os.path.join(output_dir, "crystal_docking_summary.csv"), index=False)

    ef_records = []
    for pct in [1, 5, 10]:
        ef = calculate_ef(df, pct)
        ef_records.append({"Top%": pct, "EF": ef})
        print(f"EF@{pct}% = {ef:.3f}")
    pd.DataFrame(ef_records).to_csv(os.path.join(output_dir, "ef_summary.csv"), index=False)

    plot_score_histogram(df, out_path=os.path.join(output_dir, "binding_energy_hist.png"))

    per_ef1 = per_receptor_ef(df, top_percent=1)
    per_ef1.to_csv(os.path.join(output_dir, "per_receptor_ef1.csv"))
