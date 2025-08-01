"""
Analyze the docking output file (.dlg).
How to Run:

python script.py --kinase abl1 --max_receptor 50

"""
import os
import re
import json
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Pool
from tqdm import tqdm
import shutil

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


# === 1. 먼저 parse_receptor_dir 정의 ===
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


# --- 꼭 상단에 위치시켜야 합니다 ---
def safe_parse(args):
    try:
        return parse_receptor_dir(args)
    except Exception as e:
        print(f"[ERROR] Failed to parse {args}: {e}")
        return []


# --- collect_best_scores 내부 ---
def collect_best_scores(base_dir, max_receptor=50, checkpoint_file=None):
    receptor_args = []
    for i in range(1, max_receptor + 1):
        receptor_name = f"receptor_{i:03d}"  # ← 혹시 0001이면 04d로
        receptor_path = os.path.join(base_dir, receptor_name)
        if os.path.isdir(receptor_path):
            receptor_args.append((receptor_path, receptor_name))

    print(f"[INFO] Parsing {len(receptor_args)} receptors with multiprocessing...")

    best_scores = {}

    with Pool(processes=4) as pool:
        for result in tqdm(pool.imap(safe_parse, receptor_args), total=len(receptor_args), desc="Receptors"):
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
    sns.histplot(data=df, x="Energy", hue="Label", bins=50, kde=True, stat="density", common_norm=False)
    plt.title("Normalized Binding Energy Distribution")
    plt.xlabel("Best Binding Energy (kcal/mol)")
    plt.ylabel("Density")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"[INFO] Histogram saved to {out_path}")

def per_receptor_ef(df, top_percent=1):
    top_n = max(1, int(len(df) * top_percent / 100))
    top_df = df.sort_values("Energy").head(top_n)
    return top_df[top_df["Label"] == "active"]["Best_Receptor"].value_counts()



def find_and_copy_top_receptors(per_ef1, kinase, output_dir, pdb_base_dir=None, copy_files=True):
    """
    Find and optionally copy the top receptor PDB files based on EF@1% performance.
    
    Args:
        per_ef1: Series with receptor performance counts
        kinase: kinase name
        output_dir: directory to save results
        pdb_base_dir: base directory containing PDB files (if None, will use default path)
        copy_files: whether to copy files or just create a list
    """
    
    # Get top 10 receptors
    top_receptors = per_ef1.head(10).index.tolist()
    
    # Create subdirectory for top receptors
    top_receptors_dir = os.path.join(output_dir, "top_10_receptors")
    if copy_files:
        os.makedirs(top_receptors_dir, exist_ok=True)
    
    # Use default PDB directory structure if not provided
    if pdb_base_dir is None:
        pdb_base_dir = f"/store/jaeohshin/work/dock/virtual_screening/input/receptors/{kinase}"
    
    # Check if the directory exists
    if not os.path.exists(pdb_base_dir):
        print(f"[ERROR] PDB directory not found: {pdb_base_dir}")
        return None, [], []
    
    # Find and process receptor files
    receptor_info = []
    found_files = []
    missing_files = []
    
    for receptor in top_receptors:
        count = per_ef1[receptor]
        
        # Direct path to PDB file (based on your structure)
        pdb_file = os.path.join(pdb_base_dir, f"{receptor}.pdb")
        
        if os.path.exists(pdb_file):
            found_files.append(pdb_file)
            
            if copy_files:
                # Copy to top receptors directory with descriptive name
                dest_file = os.path.join(top_receptors_dir, f"{receptor}_EF1_{count}hits.pdb")
                shutil.copy2(pdb_file, dest_file)
                print(f"[INFO] Copied {receptor} (EF@1%: {count} hits) -> {dest_file}")
            
            receptor_info.append({
                "Receptor": receptor,
                "EF1_Hits": count,
                "Source_File": pdb_file,
                "Status": "Found"
            })
        else:
            missing_files.append(receptor)
            receptor_info.append({
                "Receptor": receptor,
                "EF1_Hits": count,
                "Source_File": "NOT_FOUND",
                "Status": "Missing"
            })
            print(f"[WARN] PDB file not found: {pdb_file}")
    
    # Save receptor information
    receptor_df = pd.DataFrame(receptor_info)
    receptor_df.to_csv(os.path.join(output_dir, "top_10_receptors_info.csv"), index=False)
    
    # Print summary
    print(f"\n[INFO] Top 10 Receptor Summary:")
    print(f"Found: {len(found_files)} PDB files")
    print(f"Missing: {len(missing_files)} PDB files")
    
    if missing_files:
        print(f"Missing receptors: {missing_files}")
        print(f"Searched in: {pdb_base_dir}")
    
    return receptor_df, found_files, missing_files


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--kinase", required=True, help="Kinase target name (e.g., abl1)")
    parser.add_argument("--max_receptor", type=int, default=50,
                        help="Max receptor conformation index to include (default: 100)")
    args = parser.parse_args()

    kinase = args.kinase.lower()
    max_r = args.max_receptor

    docking_dir = f"../virtual_screening/docking_output/{kinase}"
    output_dir = f"../virtual_screening/results/{kinase}_upto_{max_r:04d}"
    os.makedirs(output_dir, exist_ok=True)
    checkpoint_path = os.path.join(output_dir, "checkpoint.json")

    df = collect_best_scores(docking_dir, max_receptor=max_r, checkpoint_file=checkpoint_path)
    print(f"[INFO] {kinase.upper()} up to receptor_{max_r:04d} — Ligands: {len(df)}, Actives: {(df['Label']=='active').sum()}, Decoys: {(df['Label']=='decoy').sum()}")

    df.to_csv(os.path.join(output_dir, "ensemble_docking_summary.csv"), index=False)

    ef_records = []
    for pct in [1, 5, 10]:
        ef = calculate_ef(df, pct)
        ef_records.append({"Top%": pct, "EF": ef})
        print(f"EF@{pct}% = {ef:.3f}")
    pd.DataFrame(ef_records).to_csv(os.path.join(output_dir, "ef_summary.csv"), index=False)

    plot_score_histogram(df, out_path=os.path.join(output_dir, "binding_energy_hist.png"))


    per_ef1 = per_receptor_ef(df, top_percent=1)
    per_ef1.to_csv(os.path.join(output_dir, "per_receptor_ef1.csv"))

    print("\n[INFO] Top 10 receptors contributing to EF@1%:")
    print(per_ef1.head(10))
        # NEW: Find and copy top receptor PDB files
    receptor_df, found_files, missing_files = find_and_copy_top_receptors(
        per_ef1=per_ef1,
        kinase=kinase,
        output_dir=output_dir,
        pdb_base_dir=None,  # Will auto-detect, or specify your PDB directory
        copy_files=True     # Set to False if you only want the file list
    )

