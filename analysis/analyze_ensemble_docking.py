import os
import re
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Pool, cpu_count
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

def parse_receptor_dir(args):
    receptor_path, receptor_name = args
    records = []
    for fname in os.listdir(receptor_path):
        if not fname.endswith(".dlg"):
            continue
        ligand = fname[:-4]
        fpath = os.path.join(receptor_path, fname)
        energy = extract_best_energy_from_dlg(fpath)
        if energy is not None:
            records.append((ligand, energy, label_ligand(fname), receptor_name))
    return records

def collect_best_scores(base_dir, max_receptor=100):
    receptor_args = []
    for i in range(1, max_receptor + 1):
        receptor_name = f"receptor_{i:04d}"
        receptor_path = os.path.join(base_dir, receptor_name)
        if os.path.isdir(receptor_path):
            receptor_args.append((receptor_path, receptor_name))

    print(f"[INFO] Parsing {len(receptor_args)} receptors with multiprocessing...")
    best_scores = {}

    with Pool(processes=4) as pool:
        for result in tqdm(pool.imap(parse_receptor_dir, receptor_args), total=len(receptor_args), desc="Receptors"):
            for ligand, energy, label, receptor in result:
                if (ligand not in best_scores) or (energy < best_scores[ligand][0]):
                    best_scores[ligand] = (energy, label, receptor)

    df = pd.DataFrame([(lig, en, lbl, rec) for lig, (en, lbl, rec) in best_scores.items()],
                      columns=["Ligand", "Energy", "Label", "Best_Receptor"])
    return df

def calculate_ef(df, top_percent):
    df_sorted = df.sort_values("Energy")
    total = len(df)
    actives_total = (df["Label"] == "active").sum()
    top_n = max(1, int(total * top_percent / 100))
    top_df = df_sorted.head(top_n)
    actives_top = (top_df["Label"] == "active").sum()
    ef = (actives_top / top_n) / (actives_total / total)
    return ef


def plot_score_histogram(df, out_path):
    plt.figure(figsize=(8, 4))
    sns.histplot(data=df, x="Energy", hue="Label", bins=50, kde=True)
    plt.title("Binding Energy Distribution")
    plt.xlabel("Best Binding Energy (kcal/mol)")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"[INFO] Histogram saved to {out_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--kinase", required=True, help="Kinase target name (e.g., abl1)")
    parser.add_argument("--max_receptor", type=int, default=100,
                        help="Max receptor conformation index to include (default: 100)")
    args = parser.parse_args()

    kinase = args.kinase.lower()
    max_r = args.max_receptor

    docking_dir = f"../virtual_screening/docking_output/{kinase}"
    output_dir = f"../virtual_screening/results/{kinase}_upto_{max_r:04d}"
    os.makedirs(output_dir, exist_ok=True)

    df = collect_best_scores(docking_dir, max_receptor=max_r)
    print(f"[INFO] {kinase.upper()} up to receptor_{max_r:04d} — Ligands: {len(df)}, Actives: {(df['Label']=='active').sum()}, Decoys: {(df['Label']=='decoy').sum()}")

    df.to_csv(os.path.join(output_dir, "ensemble_docking_summary.csv"), index=False)

    ef_records = []
    for pct in [1, 5, 10]:
        ef = calculate_ef(df, pct)
        ef_records.append({"Top%": pct, "EF": ef})
        print(f"EF@{pct}% = {ef:.3f}")
    pd.DataFrame(ef_records).to_csv(os.path.join(output_dir, "ef_summary.csv"), index=False)

    plot_score_histogram(df, out_path=os.path.join(output_dir, "binding_energy_hist.png"))
    
    
    # --- Per-receptor EF@1% contribution ---
    def per_receptor_ef(df, top_percent=1):
        top_n = max(1, int(len(df) * top_percent / 100))
        top_df = df.sort_values("Energy").head(top_n)
        return top_df[top_df["Label"] == "active"]["Best_Receptor"].value_counts()

    per_ef1 = per_receptor_ef(df, top_percent=1)
    per_ef1.to_csv(os.path.join(output_dir, "per_receptor_ef1.csv"))

    print("\n[INFO] Top 10 receptors contributing to EF@1%:")
    print(per_ef1.head(10))


