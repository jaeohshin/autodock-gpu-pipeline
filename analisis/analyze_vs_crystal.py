"""
Combined Virtual Screening Analysis for Crystal Structures
Analyzes docking results for all kinases and generates summary statistics.

Usage:
    python analyze_vs_crystal.py                    # Process missing kinases only
    python analyze_vs_crystal.py --force            # Recompute all kinases
    python analyze_vs_crystal.py --kinase abl1      # Process specific kinase only

Author: Jaeoh
"""

import os
import re
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from collections import defaultdict

# ======================== Configuration ========================
BASE_DIR = "/store/jaeohshin/work/dock"
KINASE_FILE = os.path.join(BASE_DIR, "kinase.txt")
DOCKING_BASE = os.path.join(BASE_DIR, "vs_crystal/docking_output")
RESULTS_BASE = os.path.join(BASE_DIR, "vs_crystal/results")

# ======================== Helper Functions ========================
def extract_best_energy_from_dlg(filepath):
    """Extract minimum binding energy from .dlg file"""
    energies = []
    with open(filepath, 'r') as f:
        for line in f:
            if "Estimated Free Energy of Binding" in line:
                match = re.search(r"=\s*(-?\d+\.\d+)", line)
                if match:
                    energies.append(float(match.group(1)))
    return min(energies) if energies else None

def label_ligand(filename):
    """Determine if ligand is active or decoy based on filename"""
    if filename.startswith("actives_"):
        return "active"
    elif filename.startswith("decoys_"):
        return "decoy"
    return "unknown"

def collect_best_scores(docking_dir):
    """Parse all .dlg files in directory and extract best scores per ligand"""
    energy_dict = defaultdict(list)
    label_dict = {}

    if not os.path.isdir(docking_dir):
        raise FileNotFoundError(f"[ERROR] Directory not found: {docking_dir}")

    dlg_files = [f for f in os.listdir(docking_dir) if f.endswith(".dlg")]
    
    for fname in tqdm(dlg_files, desc="Parsing ligands", leave=False):
        match = re.search(r'REMARKName=(ZINC\d+|CHEMBL\d+|C\d+|\d+)', fname)
        if not match:
            continue
        ligand_id = match.group(1)
        fpath = os.path.join(docking_dir, fname)
        energy = extract_best_energy_from_dlg(fpath)
        if energy is not None:
            energy_dict[ligand_id].append(energy)
            if ligand_id not in label_dict:
                label_dict[ligand_id] = label_ligand(fname)
    
    records = []
    for ligand_id, energies in energy_dict.items():
        best_energy = min(energies)
        label = label_dict.get(ligand_id, "unknown")
        records.append((ligand_id, best_energy, label))

    df = pd.DataFrame(records, columns=["Ligand", "Energy", "Label"])
    df["Best_Receptor"] = "crystal"
    return df

def calculate_ef(df, top_percent):
    """Calculate enrichment factor at given percentage"""
    df_sorted = df.sort_values("Energy")
    total = len(df)
    actives_total = (df["Label"] == "active").sum()

    if total == 0 or actives_total == 0:
        return float("nan")

    top_n = max(1, int(total * top_percent / 100))
    actives_top = (df_sorted.head(top_n)["Label"] == "active").sum()
    ef = (actives_top / top_n) / (actives_total / total)
    return ef

def plot_score_histogram(df, out_path):
    """Generate histogram of binding energy distribution"""
    if df.empty:
        print("    [WARN] DataFrame is empty. Skipping histogram.")
        return

    kde_flag = df["Energy"].nunique() >= 2
    
    plt.figure(figsize=(8, 4))
    sns.histplot(data=df, x="Energy", hue="Label", bins=50, kde=kde_flag)
    plt.title("Binding Energy Distribution")
    plt.xlabel("Best Binding Energy (kcal/mol)")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()

def per_receptor_ef(df, top_percent=1):
    """Calculate per-receptor contribution to EF"""
    top_n = max(1, int(len(df) * top_percent / 100))
    top_df = df.sort_values("Energy").head(top_n)
    return top_df[top_df["Label"] == "active"]["Best_Receptor"].value_counts()

def process_kinase(kinase_name, force=False):
    """Process single kinase: parse docking results and calculate EFs"""
    kinase_lower = kinase_name.lower()
    docking_dir = os.path.join(DOCKING_BASE, kinase_lower)
    output_dir = os.path.join(RESULTS_BASE, kinase_lower)
    
    os.makedirs(output_dir, exist_ok=True)
    
    ef_file = os.path.join(output_dir, "ef_summary.csv")
    
    # Skip if already processed (unless force flag)
    if os.path.exists(ef_file) and not force:
        print(f"  ✓ {kinase_name}: Already processed (use --force to recompute)")
        return True
    
    # Check if docking directory exists
    if not os.path.isdir(docking_dir):
        print(f"  ✗ {kinase_name}: Docking directory not found: {docking_dir}")
        return False
    
    try:
        # Parse docking results
        df = collect_best_scores(docking_dir)
        n_actives = (df['Label'] == 'active').sum()
        n_decoys = (df['Label'] == 'decoy').sum()
        
        print(f"  → {kinase_name}: Ligands={len(df)}, Actives={n_actives}, Decoys={n_decoys}")
        
        # Save detailed results
        df.to_csv(os.path.join(output_dir, "crystal_docking_summary.csv"), index=False)
        
        # Calculate enrichment factors
        ef_records = []
        for pct in [1, 5, 10]:
            ef = calculate_ef(df, pct)
            ef_records.append({"Top%": pct, "EF": ef})
        
        ef_df = pd.DataFrame(ef_records)
        ef_df.to_csv(ef_file, index=False)
        
        # Generate plots
        plot_score_histogram(df, out_path=os.path.join(output_dir, "binding_energy_hist.png"))
        
        # Per-receptor EF analysis
        per_ef1 = per_receptor_ef(df, top_percent=1)
        per_ef1.to_csv(os.path.join(output_dir, "per_receptor_ef1.csv"))
        
        # Print EF summary
        ef_str = ", ".join([f"EF@{row['Top%']}%={row['EF']:.2f}" for _, row in ef_df.iterrows()])
        print(f"  ✓ {kinase_name}: {ef_str}")
        
        return True
        
    except Exception as e:
        print(f"  ✗ {kinase_name}: ERROR - {str(e)}")
        return False

def aggregate_results():
    """Aggregate all kinase EF results into master summary"""
    print("\n" + "="*60)
    print("Aggregating results from all kinases...")
    print("="*60)
    
    # Read kinase list
    kinases = []
    with open(KINASE_FILE, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 1:
                kinases.append(parts[0])
    
    # Collect EF data
    ef_data = {
        'Kinase': [],
        'EF_1%': [],
        'EF_5%': [],
        'EF_10%': []
    }
    
    for kinase in kinases:
        kinase_lower = kinase.lower()
        ef_file = os.path.join(RESULTS_BASE, kinase_lower, "ef_summary.csv")
        
        if os.path.exists(ef_file):
            df = pd.read_csv(ef_file)
            ef_dict = {}
            for _, row in df.iterrows():
                top_pct = row['Top%']
                ef_value = row['EF']
                ef_dict[top_pct] = ef_value
            
            if 1 in ef_dict and 5 in ef_dict and 10 in ef_dict:
                ef_data['Kinase'].append(kinase)
                ef_data['EF_1%'].append(ef_dict[1])
                ef_data['EF_5%'].append(ef_dict[5])
                ef_data['EF_10%'].append(ef_dict[10])
            else:
                print(f"  ✗ {kinase}: Missing some EF values")
        else:
            print(f"  ✗ {kinase}: ef_summary.csv not found")
    
    # Create summary DataFrame
    df_summary = pd.DataFrame(ef_data)
    
    # Round to 3 decimal places
    for col in ['EF_1%', 'EF_5%', 'EF_10%']:
        df_summary[col] = df_summary[col].round(3)
    
    # Save master summary
    output_file = os.path.join(RESULTS_BASE, "all_kinases_ef_summary.csv")
    df_summary.to_csv(output_file, index=False)
    
    print(f"\n{'='*60}")
    print(f"Summary saved to: {output_file}")
    print(f"Total kinases processed: {len(df_summary)}")
    print(f"{'='*60}\n")
    print(df_summary.to_string(index=False))

# ======================== Main ========================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze virtual screening results for crystal structures")
    parser.add_argument("--kinase", help="Process specific kinase only (e.g., abl1)")
    parser.add_argument("--force", action="store_true", help="Force recomputation of existing results")
    args = parser.parse_args()
    
    # Read kinase list
    kinases = []
    with open(KINASE_FILE, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 1:
                kinases.append(parts[0])
    
    print("="*60)
    print(f"Virtual Screening Analysis - Crystal Structures")
    print(f"Total kinases: {len(kinases)}")
    print("="*60 + "\n")
    
    # Process kinases
    if args.kinase:
        # Single kinase mode
        kinase_upper = args.kinase.upper()
        if kinase_upper in kinases:
            process_kinase(kinase_upper, force=args.force)
        else:
            print(f"[ERROR] Kinase '{args.kinase}' not found in kinase.txt")
    else:
        # Process all kinases
        success_count = 0
        for kinase in tqdm(kinases, desc="Processing kinases"):
            if process_kinase(kinase, force=args.force):
                success_count += 1
        
        print(f"\n{'='*60}")
        print(f"Processing complete: {success_count}/{len(kinases)} kinases successful")
        print(f"{'='*60}\n")
    
    # Aggregate results
    aggregate_results()
