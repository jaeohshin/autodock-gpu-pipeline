"""
Modified version of your analysis code to use average scores and TM-score filtering.
Key changes marked with # MODIFIED comments.
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
import numpy as np

# MODIFIED: New function to parse TM-score data
def parse_tmscore_file(tmscore_file, min_tmscore=0.8):
    """
    Parse the ranked_report.txt file to get approved receptors based on TM-score threshold.
    
    Args:
        tmscore_file: Path to ranked_report.txt file
        min_tmscore: Minimum TM-score threshold (default: 0.8)
    
    Returns:
        dict: Dictionary with kinase -> list of approved receptor names
    """
    approved_receptors = {}
    current_kinase = None
    
    with open(tmscore_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            # Check if this line is a kinase name (no comma, no numbers)
            if ',' not in line and not line[0].isdigit():
                current_kinase = line.lower()
                approved_receptors[current_kinase] = []
            elif current_kinase and ',' in line:
                # Parse receptor line: "1, receptor_040.pdb, 0.8486"
                parts = [p.strip() for p in line.split(',')]
                if len(parts) >= 3:
                    try:
                        rank = int(parts[0])
                        receptor_file = parts[1]
                        tmscore = float(parts[2])
                        
                        if tmscore >= min_tmscore:
                            # Extract receptor name without .pdb extension
                            receptor_name = receptor_file.replace('.pdb', '')
                            approved_receptors[current_kinase].append(receptor_name)
                    except (ValueError, IndexError):
                        continue
    
    return approved_receptors

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

def safe_parse(args):
    try:
        return parse_receptor_dir(args)
    except Exception as e:
        print(f"[ERROR] Failed to parse {args}: {e}")
        return []

# MODIFIED: Updated to use only approved receptors
def collect_all_scores(base_dir, approved_receptors, kinase, checkpoint_file=None):
    """
    Collect ALL binding scores for each ligand across approved receptor conformations only.
    
    Args:
        base_dir: Base directory containing receptor folders
        approved_receptors: List of approved receptor names for this kinase
        kinase: Kinase name
        checkpoint_file: Optional checkpoint file
    
    Returns:
        DataFrame with all ligand-receptor pairs for approved receptors only
    """
    receptor_args = []
    
    # MODIFIED: Only process approved receptors
    for receptor_name in approved_receptors:
        receptor_path = os.path.join(base_dir, receptor_name)
        if os.path.isdir(receptor_path):
            receptor_args.append((receptor_path, receptor_name))
        else:
            print(f"[WARN] Approved receptor directory not found: {receptor_path}")

    print(f"[INFO] Parsing {len(receptor_args)} approved receptors (TM-score ≥ 0.8) with multiprocessing...")
    print(f"[INFO] Approved receptors: {approved_receptors}")

    all_records = []

    with Pool(processes=4) as pool:
        for result in tqdm(pool.imap(safe_parse, receptor_args), total=len(receptor_args), desc="Receptors"):
            all_records.extend(result)

    df_all = pd.DataFrame(all_records, columns=["Ligand", "Energy", "Label", "Receptor"])
    return df_all

def compute_average_scores(df_all):
    """
    Compute average binding energy for each ligand across all receptor conformations.
    Also tracks which receptors contributed and provides statistics.
    """
    # Group by ligand and compute statistics
    ligand_stats = df_all.groupby(['Ligand', 'Label'])['Energy'].agg([
        'mean',     # Average score
        'std',      # Standard deviation
        'min',      # Best (minimum) score
        'max',      # Worst (maximum) score  
        'count',    # Number of conformations
        'median'    # Median score
    ]).reset_index()
    
    # Find which receptor gave the best score for each ligand (for comparison)
    best_receptor = df_all.loc[df_all.groupby('Ligand')['Energy'].idxmin()]
    best_receptor = best_receptor[['Ligand', 'Receptor']].rename(columns={'Receptor': 'Best_Receptor'})
    
    # Merge the statistics with best receptor info
    result_df = ligand_stats.merge(best_receptor, on='Ligand', how='left')
    
    # Rename 'mean' to 'Energy' for compatibility with existing code
    result_df = result_df.rename(columns={'mean': 'Energy'})
    
    return result_df

def compute_scores_with_method(df_all, method='mean'):
    """
    Compute ligand scores using different aggregation methods.
    
    Args:
        df_all: DataFrame with all ligand-receptor scores
        method: 'mean', 'median', 'best' (min), 'worst' (max), 'trimmed_mean', 'weighted_mean'
    
    Returns:
        DataFrame with aggregated scores
    """
    if method == 'mean':
        agg_func = 'mean'
    elif method == 'median':
        agg_func = 'median'
    elif method == 'best':
        agg_func = 'min'
    elif method == 'worst':
        agg_func = 'max'
    elif method == 'trimmed_mean':
        # Remove top and bottom 10% outliers, then take mean
        def trimmed_mean(x):
            return x.quantile([0.1, 0.9]).mean() if len(x) > 2 else x.mean()
        agg_func = trimmed_mean
    else:
        agg_func = 'mean'  # default
    
    # Group by ligand and compute the aggregated score
    if method == 'trimmed_mean':
        ligand_scores = df_all.groupby(['Ligand', 'Label'])['Energy'].apply(agg_func).reset_index()
    else:
        ligand_scores = df_all.groupby(['Ligand', 'Label'])['Energy'].agg(agg_func).reset_index()
    
    # Add receptor count and best receptor info
    receptor_counts = df_all.groupby('Ligand')['Receptor'].count().reset_index()
    receptor_counts = receptor_counts.rename(columns={'Receptor': 'Receptor_Count'})
    
    best_receptor = df_all.loc[df_all.groupby('Ligand')['Energy'].idxmin()]
    best_receptor = best_receptor[['Ligand', 'Receptor']].rename(columns={'Receptor': 'Best_Receptor'})
    
    result_df = ligand_scores.merge(receptor_counts, on='Ligand', how='left')
    result_df = result_df.merge(best_receptor, on='Ligand', how='left')
    
    return result_df

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

def plot_score_comparison(df_all, output_dir):
    """
    Plot comparison of different scoring methods.
    """
    methods = ['best', 'mean', 'median', 'trimmed_mean']
    method_dfs = {}
    
    for method in methods:
        method_dfs[method] = compute_scores_with_method(df_all, method)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.ravel()
    
    for i, method in enumerate(methods):
        df_method = method_dfs[method]
        if not df_method.empty and df_method["Energy"].nunique() > 1:
            sns.histplot(data=df_method, x="Energy", hue="Label", bins=30, 
                        kde=True, stat="density", ax=axes[i])
            axes[i].set_title(f"{method.title()} Score Distribution")
            axes[i].set_xlabel("Binding Energy (kcal/mol)")
            axes[i].set_ylabel("Density")
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "scoring_methods_comparison.png"), dpi=150)
    plt.close()
    
    return method_dfs

def plot_score_histogram(df, out_path, method_name):
    if df.empty:
        print("[WARN] DataFrame is empty. Skipping histogram.")
        return

    plt.figure(figsize=(8, 4))
    sns.histplot(data=df, x="Energy", hue="Label", bins=50, kde=True, stat="density", common_norm=False)
    plt.title(f"Binding Energy Distribution ({method_name.title()} Scores, TM-score ≥ 0.8)")  # MODIFIED: Updated title
    plt.xlabel(f"{method_name.title()} Binding Energy (kcal/mol)")  # MODIFIED: Updated label
    plt.ylabel("Density")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"[INFO] Histogram saved to {out_path}")

def per_receptor_contribution(df_all, df_avg, top_percent=1):
    """
    Analyze which receptors contribute most to top-scoring ligands when using average scores.
    """
    # Get top ligands based on average scores
    top_n = max(1, int(len(df_avg) * top_percent / 100))
    top_ligands = df_avg.sort_values("Energy").head(top_n)["Ligand"].tolist()
    
    # Find all receptor contributions for these top ligands
    top_ligand_data = df_all[df_all["Ligand"].isin(top_ligands)]
    
    # Count receptor frequency in top ligands
    receptor_counts = top_ligand_data["Receptor"].value_counts()
    
    return receptor_counts

def find_and_copy_top_receptors(per_ef1, kinase, output_dir, pdb_base_dir=None, copy_files=True):
    """Find and optionally copy the top receptor PDB files based on performance."""
    
    top_receptors = per_ef1.head(10).index.tolist()
    top_receptors_dir = os.path.join(output_dir, "top_10_receptors")
    if copy_files:
        os.makedirs(top_receptors_dir, exist_ok=True)
    
    if pdb_base_dir is None:
        pdb_base_dir = f"/store/jaeohshin/work/dock/virtual_screening/input/receptors/{kinase}"
    
    if not os.path.exists(pdb_base_dir):
        print(f"[ERROR] PDB directory not found: {pdb_base_dir}")
        return None, [], []
    
    receptor_info = []
    found_files = []
    missing_files = []
    
    for receptor in top_receptors:
        count = per_ef1[receptor]
        pdb_file = os.path.join(pdb_base_dir, f"{receptor}.pdb")
        
        if os.path.exists(pdb_file):
            found_files.append(pdb_file)
            
            if copy_files:
                dest_file = os.path.join(top_receptors_dir, f"{receptor}_EF1_{count}hits.pdb")
                shutil.copy2(pdb_file, dest_file)
                print(f"[INFO] Copied {receptor} (EF@1%: {count} hits) -> {dest_file}")
            
            receptor_info.append({
                "Receptor": receptor,
                "EF1_Contribution": count,
                "Source_File": pdb_file,
                "Status": "Found"
            })
        else:
            missing_files.append(receptor)
            receptor_info.append({
                "Receptor": receptor,
                "EF1_Contribution": count,
                "Source_File": "NOT_FOUND",
                "Status": "Missing"
            })
            print(f"[WARN] PDB file not found: {pdb_file}")
    
    receptor_df = pd.DataFrame(receptor_info)
    receptor_df.to_csv(os.path.join(output_dir, "top_10_receptors_info.csv"), index=False)
    
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
    parser.add_argument("--tmscore_file", required=True, help="Path to ranked_report.txt file with TM-scores")  # MODIFIED: Added TM-score file
    parser.add_argument("--min_tmscore", type=float, default=0.8, help="Minimum TM-score threshold (default: 0.8)")  # MODIFIED: Added TM-score threshold
    parser.add_argument("--method", choices=['best', 'mean', 'median', 'trimmed_mean'], 
                        default='mean', help="Scoring method to use (default: mean)")
    args = parser.parse_args()

    kinase = args.kinase.lower()
    scoring_method = args.method
    tmscore_file = args.tmscore_file
    min_tmscore = args.min_tmscore

    # MODIFIED: Parse TM-score file to get approved receptors
    print(f"[INFO] Parsing TM-score file: {tmscore_file}")
    print(f"[INFO] Using minimum TM-score threshold: {min_tmscore}")
    
    approved_receptors_dict = parse_tmscore_file(tmscore_file, min_tmscore)
    
    if kinase not in approved_receptors_dict:
        print(f"[ERROR] Kinase '{kinase}' not found in TM-score file")
        print(f"Available kinases: {list(approved_receptors_dict.keys())}")
        exit(1)
    
    approved_receptors = approved_receptors_dict[kinase]
    
    if not approved_receptors:
        print(f"[ERROR] No receptors found for {kinase} with TM-score ≥ {min_tmscore}")
        exit(1)
    
    print(f"[INFO] Found {len(approved_receptors)} approved receptors for {kinase}")

    docking_dir = f"/store/jaeohshin/work/dock/virtual_screening/docking_output/{kinase}"
    output_dir = f"../virtual_screening/results/{kinase}_tmscore{min_tmscore}_{scoring_method}"  # MODIFIED: Updated output dir name
    os.makedirs(output_dir, exist_ok=True)

    # MODIFIED: Collect scores only from approved receptors
    print(f"[INFO] Using {scoring_method} scoring method with TM-score ≥ {min_tmscore} filter")
    df_all = collect_all_scores(docking_dir, approved_receptors, kinase)
    
    if df_all.empty:
        print(f"[ERROR] No docking data found for approved receptors")
        exit(1)
    
    # Compute scores using specified method
    df = compute_scores_with_method(df_all, method=scoring_method)
    
    print(f"[INFO] {kinase.upper()} with {len(approved_receptors)} approved receptors (TM-score ≥ {min_tmscore}) using {scoring_method} — Ligands: {len(df)}, Actives: {(df['Label']=='active').sum()}, Decoys: {(df['Label']=='decoy').sum()}")

    # Save detailed results
    df.to_csv(os.path.join(output_dir, f"ensemble_docking_summary_{scoring_method}.csv"), index=False)
    df_all.to_csv(os.path.join(output_dir, "all_scores_detailed.csv"), index=False)
    
    # MODIFIED: Save approved receptors info
    approved_df = pd.DataFrame({
        'Receptor': approved_receptors,
        'TM_Score_Threshold': min_tmscore,
        'Kinase': kinase
    })
    approved_df.to_csv(os.path.join(output_dir, "approved_receptors.csv"), index=False)

    # Calculate EF values
    ef_records = []
    for pct in [1, 5, 10]:
        ef = calculate_ef(df, pct)
        ef_records.append({"Top%": pct, "EF": ef, "Method": scoring_method, "TM_Score_Threshold": min_tmscore})  # MODIFIED: Track threshold
        print(f"EF@{pct}% ({scoring_method}, TM-score ≥ {min_tmscore}) = {ef:.3f}")
    pd.DataFrame(ef_records).to_csv(os.path.join(output_dir, f"ef_summary_{scoring_method}.csv"), index=False)

    # Generate plots
    plot_score_histogram(df, out_path=os.path.join(output_dir, f"binding_energy_hist_{scoring_method}.png"), method_name=scoring_method)
    
    # Generate method comparison plots
    method_comparison = plot_score_comparison(df_all, output_dir)
    
    # Analyze receptor contributions
    per_ef1 = per_receptor_contribution(df_all, df, top_percent=1)
    per_ef1.to_csv(os.path.join(output_dir, f"per_receptor_ef1_{scoring_method}.csv"))

    print(f"\n[INFO] Top 10 approved receptors contributing to EF@1% ({scoring_method} method):")
    print(per_ef1.head(10))
    
    # Find and copy top receptor PDB files
    receptor_df, found_files, missing_files = find_and_copy_top_receptors(
        per_ef1=per_ef1,
        kinase=kinase,
        output_dir=output_dir,
        pdb_base_dir=None,
        copy_files=True
    )
    
    # Save comparison of different methods
    if scoring_method == 'mean':  # Only do this once
        method_ef_comparison = []
        for method_name, method_df in method_comparison.items():
            for pct in [1, 5, 10]:
                ef = calculate_ef(method_df, pct)
                method_ef_comparison.append({
                    "Method": method_name,
                    "Top%": pct,
                    "EF": ef,
                    "TM_Score_Threshold": min_tmscore
                })
        
        comparison_df = pd.DataFrame(method_ef_comparison)
        comparison_df.to_csv(os.path.join(output_dir, "method_comparison_ef.csv"), index=False)
        
        print(f"\n[INFO] Method Comparison (EF@1%, TM-score ≥ {min_tmscore}):")
        ef1_comparison = comparison_df[comparison_df["Top%"] == 1]
        for _, row in ef1_comparison.iterrows():
            print(f"{row['Method']:12s}: EF@1% = {row['EF']:.3f}")
    
    print(f"\n[INFO] Analysis complete! Results saved to: {output_dir}")
    print(f"[INFO] Used {len(approved_receptors)} receptors with TM-score ≥ {min_tmscore}")