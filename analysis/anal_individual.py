"""
Modified version to calculate EF for each individual receptor.
Key additions: per-receptor EF calculation and analysis functions.
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
    
    # Only process approved receptors
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

def calculate_ef(df, top_percent):
    """Calculate Enrichment Factor for a given top percentage."""
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

# NEW: Function to calculate EF for each individual receptor
def calculate_per_receptor_ef(df_all, top_percentages=[1, 5, 10]):
    """
    Calculate EF for each individual receptor.
    
    Args:
        df_all: DataFrame with all ligand-receptor scores
        top_percentages: List of top percentages to calculate EF for
    
    Returns:
        DataFrame with EF values for each receptor
    """
    per_receptor_ef = []
    
    receptors = df_all['Receptor'].unique()
    print(f"[INFO] Calculating EF for {len(receptors)} individual receptors...")
    
    for receptor in tqdm(receptors, desc="Computing per-receptor EF"):
        receptor_data = df_all[df_all['Receptor'] == receptor]
        
        # Basic statistics
        total_ligands = len(receptor_data)
        total_actives = (receptor_data['Label'] == 'active').sum()
        total_decoys = (receptor_data['Label'] == 'decoy').sum()
        
        # Calculate EF for each percentage
        ef_results = {'Receptor': receptor, 
                     'Total_Ligands': total_ligands,
                     'Total_Actives': total_actives, 
                     'Total_Decoys': total_decoys}
        
        for pct in top_percentages:
            ef_value = calculate_ef(receptor_data, pct)
            ef_results[f'EF@{pct}%'] = ef_value
        
        # Additional metrics
        if not receptor_data.empty:
            ef_results['Best_Energy'] = receptor_data['Energy'].min()
            ef_results['Worst_Energy'] = receptor_data['Energy'].max()
            ef_results['Mean_Energy'] = receptor_data['Energy'].mean()
            ef_results['Std_Energy'] = receptor_data['Energy'].std()
            
            # Calculate active vs decoy energy statistics
            active_data = receptor_data[receptor_data['Label'] == 'active']
            decoy_data = receptor_data[receptor_data['Label'] == 'decoy']
            
            if not active_data.empty:
                ef_results['Active_Mean_Energy'] = active_data['Energy'].mean()
                ef_results['Active_Best_Energy'] = active_data['Energy'].min()
            else:
                ef_results['Active_Mean_Energy'] = float('nan')
                ef_results['Active_Best_Energy'] = float('nan')
                
            if not decoy_data.empty:
                ef_results['Decoy_Mean_Energy'] = decoy_data['Energy'].mean()
                ef_results['Decoy_Best_Energy'] = decoy_data['Energy'].min()
            else:
                ef_results['Decoy_Mean_Energy'] = float('nan')
                ef_results['Decoy_Best_Energy'] = float('nan')
        
        per_receptor_ef.append(ef_results)
    
    return pd.DataFrame(per_receptor_ef)

# NEW: Function to analyze and visualize per-receptor EF results
def analyze_per_receptor_ef(per_receptor_ef_df, output_dir, top_n=10):
    """
    Analyze and visualize per-receptor EF results.
    
    Args:
        per_receptor_ef_df: DataFrame with per-receptor EF values
        output_dir: Output directory for plots and files
        top_n: Number of top receptors to highlight
    """
    
    # Sort by EF@1% and get top performers
    ef1_col = 'EF@1%'
    top_receptors = per_receptor_ef_df.nlargest(top_n, ef1_col)
    
    print(f"\n[INFO] Top {top_n} Receptors by EF@1%:")
    print(top_receptors[['Receptor', ef1_col, 'EF@5%', 'EF@10%', 'Total_Actives', 'Total_Ligands']].to_string(index=False))
    
    # Save detailed results
    per_receptor_ef_df.to_csv(os.path.join(output_dir, "per_receptor_ef_detailed.csv"), index=False)
    top_receptors.to_csv(os.path.join(output_dir, f"top_{top_n}_receptors_by_ef1.csv"), index=False)
    
    # Create visualizations
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # 1. EF@1% distribution
    axes[0,0].hist(per_receptor_ef_df[ef1_col].dropna(), bins=20, alpha=0.7, edgecolor='black')
    axes[0,0].axvline(per_receptor_ef_df[ef1_col].mean(), color='red', linestyle='--', 
                      label=f'Mean: {per_receptor_ef_df[ef1_col].mean():.2f}')
    axes[0,0].set_xlabel('EF@1%')
    axes[0,0].set_ylabel('Number of Receptors')
    axes[0,0].set_title('Distribution of EF@1% Across Receptors')
    axes[0,0].legend()
    
    # 2. EF comparison across percentages (top receptors only)
    top_receptors_melted = top_receptors[['Receptor', 'EF@1%', 'EF@5%', 'EF@10%']].melt(
        id_vars='Receptor', var_name='EF_Type', value_name='EF_Value')
    sns.barplot(data=top_receptors_melted, x='Receptor', y='EF_Value', hue='EF_Type', ax=axes[0,1])
    axes[0,1].set_title(f'Top {top_n} Receptors: EF Comparison')
    axes[0,1].tick_params(axis='x', rotation=45)
    
    # 3. Correlation between EF@1% and EF@5%
    valid_data = per_receptor_ef_df.dropna(subset=['EF@1%', 'EF@5%'])
    axes[1,0].scatter(valid_data['EF@1%'], valid_data['EF@5%'], alpha=0.6)
    axes[1,0].set_xlabel('EF@1%')
    axes[1,0].set_ylabel('EF@5%')
    axes[1,0].set_title('EF@1% vs EF@5% Correlation')
    
    # Add correlation coefficient
    if len(valid_data) > 1:
        corr_coef = valid_data['EF@1%'].corr(valid_data['EF@5%'])
        axes[1,0].text(0.05, 0.95, f'R = {corr_coef:.3f}', transform=axes[1,0].transAxes, 
                       bbox=dict(boxstyle="round", facecolor='wheat', alpha=0.5))
    
    # 4. Active recovery rate (best energy vs mean energy)
    if 'Active_Best_Energy' in per_receptor_ef_df.columns and 'Active_Mean_Energy' in per_receptor_ef_df.columns:
        valid_energy = per_receptor_ef_df.dropna(subset=['Active_Best_Energy', 'Active_Mean_Energy'])
        axes[1,1].scatter(valid_energy['Active_Best_Energy'], valid_energy['Active_Mean_Energy'], alpha=0.6)
        axes[1,1].set_xlabel('Best Active Energy (kcal/mol)')
        axes[1,1].set_ylabel('Mean Active Energy (kcal/mol)')
        axes[1,1].set_title('Active Ligand Energy Distribution')
    else:
        axes[1,1].text(0.5, 0.5, 'Active energy data not available', 
                       ha='center', va='center', transform=axes[1,1].transAxes)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "per_receptor_ef_analysis.png"), dpi=150, bbox_inches='tight')
    plt.close()
    
    # Create a summary statistics table
    summary_stats = {
        'Total_Receptors': len(per_receptor_ef_df),
        'Mean_EF@1%': per_receptor_ef_df['EF@1%'].mean(),
        'Median_EF@1%': per_receptor_ef_df['EF@1%'].median(),
        'Std_EF@1%': per_receptor_ef_df['EF@1%'].std(),
        'Max_EF@1%': per_receptor_ef_df['EF@1%'].max(),
        'Min_EF@1%': per_receptor_ef_df['EF@1%'].min(),
        'Mean_EF@5%': per_receptor_ef_df['EF@5%'].mean(),
        'Mean_EF@10%': per_receptor_ef_df['EF@10%'].mean(),
    }
    
    summary_df = pd.DataFrame([summary_stats])
    summary_df.to_csv(os.path.join(output_dir, "per_receptor_ef_summary_stats.csv"), index=False)
    
    print(f"\n[INFO] Per-receptor EF analysis summary:")
    for key, value in summary_stats.items():
        if pd.notna(value):
            print(f"{key}: {value:.3f}")
    
    return top_receptors, summary_stats

# NEW: Function to compare per-receptor EF with ensemble EF
def compare_individual_vs_ensemble_ef(df_all, per_receptor_ef_df, output_dir):
    """
    Compare individual receptor EF performance with ensemble methods.
    """
    from scipy import stats
    
    # Calculate ensemble EF using different methods
    ensemble_methods = ['best', 'mean', 'median']
    ensemble_results = []
    
    for method in ensemble_methods:
        if method == 'best':
            # Take best (minimum) energy for each ligand across all receptors
            ensemble_scores = df_all.groupby(['Ligand', 'Label'])['Energy'].min().reset_index()
        elif method == 'mean':
            ensemble_scores = df_all.groupby(['Ligand', 'Label'])['Energy'].mean().reset_index()
        elif method == 'median':
            ensemble_scores = df_all.groupby(['Ligand', 'Label'])['Energy'].median().reset_index()
        
        # Calculate EF for ensemble
        for pct in [1, 5, 10]:
            ef = calculate_ef(ensemble_scores, pct)
            ensemble_results.append({
                'Method': f'Ensemble_{method}',
                'EF_Type': f'EF@{pct}%',
                'EF_Value': ef
            })
    
    ensemble_df = pd.DataFrame(ensemble_results)
    
    # Compare with individual receptor performance
    comparison_data = []
    
    # Add individual receptor results
    for _, row in per_receptor_ef_df.iterrows():
        for pct in [1, 5, 10]:
            comparison_data.append({
                'Method': 'Individual_Receptor',
                'Receptor': row['Receptor'],
                'EF_Type': f'EF@{pct}%',
                'EF_Value': row[f'EF@{pct}%']
            })
    
    # Add ensemble results
    for _, row in ensemble_df.iterrows():
        comparison_data.append({
            'Method': row['Method'],
            'Receptor': 'Ensemble',
            'EF_Type': row['EF_Type'],
            'EF_Value': row['EF_Value']
        })
    
    comparison_df = pd.DataFrame(comparison_data)
    comparison_df.to_csv(os.path.join(output_dir, "individual_vs_ensemble_ef_comparison.csv"), index=False)
    
    # Create comparison plot
    plt.figure(figsize=(12, 8))
    
    # Plot individual receptor EF@1% distribution
    individual_ef1 = comparison_df[(comparison_df['Method'] == 'Individual_Receptor') & 
                                   (comparison_df['EF_Type'] == 'EF@1%')]
    plt.hist(individual_ef1['EF_Value'].dropna(), bins=20, alpha=0.6, 
             label='Individual Receptors', color='lightblue', edgecolor='black')
    
    # Add ensemble method lines
    ensemble_ef1 = comparison_df[(comparison_df['Method'].str.startswith('Ensemble')) & 
                                 (comparison_df['EF_Type'] == 'EF@1%')]
    
    colors = ['red', 'green', 'orange']
    for i, method in enumerate(['Ensemble_best', 'Ensemble_mean', 'Ensemble_median']):
        method_data = ensemble_ef1[ensemble_ef1['Method'] == method]
        if not method_data.empty:
            ef_value = method_data['EF_Value'].iloc[0]
            plt.axvline(ef_value, color=colors[i], linestyle='--', linewidth=2,
                       label=f'{method.replace("Ensemble_", "Ensemble ")} (EF@1% = {ef_value:.2f})')
    
    plt.xlabel('EF@1%')
    plt.ylabel('Number of Receptors')
    plt.title('Individual Receptor EF@1% vs Ensemble Methods')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(output_dir, "individual_vs_ensemble_ef_comparison.png"), 
                dpi=150, bbox_inches='tight')
    plt.close()
    
    # Statistical summary
    individual_ef1_values = individual_ef1['EF_Value'].dropna()
    stats_summary = {
        'Individual_Receptors_Count': len(individual_ef1_values),
        'Individual_Mean_EF@1%': individual_ef1_values.mean(),
        'Individual_Median_EF@1%': individual_ef1_values.median(),
        'Individual_Std_EF@1%': individual_ef1_values.std(),
        'Individual_Max_EF@1%': individual_ef1_values.max(),
        'Percentile_90_EF@1%': individual_ef1_values.quantile(0.9),
        'Percentile_75_EF@1%': individual_ef1_values.quantile(0.75),
        'Receptors_Better_Than_Ensemble_Best': (individual_ef1_values > 
                                                ensemble_ef1[ensemble_ef1['Method'] == 'Ensemble_best']['EF_Value'].iloc[0]).sum(),
        'Receptors_Better_Than_Ensemble_Mean': (individual_ef1_values > 
                                                ensemble_ef1[ensemble_ef1['Method'] == 'Ensemble_mean']['EF_Value'].iloc[0]).sum(),
    }
    
    # Add ensemble values to summary
    for _, row in ensemble_ef1.iterrows():
        stats_summary[f"{row['Method']}_EF@1%"] = row['EF_Value']
    
    stats_df = pd.DataFrame([stats_summary])
    stats_df.to_csv(os.path.join(output_dir, "individual_vs_ensemble_stats.csv"), index=False)
    
    print(f"\n[INFO] Individual vs Ensemble EF@1% Comparison:")
    print(f"Individual receptors - Mean: {individual_ef1_values.mean():.3f}, "
          f"Max: {individual_ef1_values.max():.3f}, "
          f"90th percentile: {individual_ef1_values.quantile(0.9):.3f}")
    
    for _, row in ensemble_ef1.iterrows():
        print(f"{row['Method']}: {row['EF_Value']:.3f}")
    
    return comparison_df, stats_summary

# Keep original ensemble functions for compatibility
def compute_scores_with_method(df_all, method='mean'):
    """Compute ligand scores using different aggregation methods."""
    if method == 'mean':
        agg_func = 'mean'
    elif method == 'median':
        agg_func = 'median'
    elif method == 'best':
        agg_func = 'min'
    elif method == 'worst':
        agg_func = 'max'
    else:
        agg_func = 'mean'
    
    ligand_scores = df_all.groupby(['Ligand', 'Label'])['Energy'].agg(agg_func).reset_index()
    
    receptor_counts = df_all.groupby('Ligand')['Receptor'].count().reset_index()
    receptor_counts = receptor_counts.rename(columns={'Receptor': 'Receptor_Count'})
    
    best_receptor = df_all.loc[df_all.groupby('Ligand')['Energy'].idxmin()]
    best_receptor = best_receptor[['Ligand', 'Receptor']].rename(columns={'Receptor': 'Best_Receptor'})
    
    result_df = ligand_scores.merge(receptor_counts, on='Ligand', how='left')
    result_df = result_df.merge(best_receptor, on='Ligand', how='left')
    
    return result_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--kinase", required=True, help="Kinase target name (e.g., abl1)")
    parser.add_argument("--tmscore_file", required=True, help="Path to ranked_report.txt file with TM-scores")
    parser.add_argument("--min_tmscore", type=float, default=0.8, help="Minimum TM-score threshold (default: 0.8)")
    parser.add_argument("--per_receptor", action='store_true', help="Calculate EF for each individual receptor")
    parser.add_argument("--top_n", type=int, default=10, help="Number of top receptors to analyze (default: 10)")
    args = parser.parse_args()

    kinase = args.kinase.lower()
    tmscore_file = args.tmscore_file
    min_tmscore = args.min_tmscore
    calculate_per_receptor = args.per_receptor
    top_n = args.top_n

    # Parse TM-score file to get approved receptors
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
    output_dir = f"../virtual_screening/results/{kinase}_per_receptor_ef_tmscore{min_tmscore}"
    os.makedirs(output_dir, exist_ok=True)

    # Collect scores from approved receptors
    print(f"[INFO] Collecting docking scores from approved receptors...")
    df_all = collect_all_scores(docking_dir, approved_receptors, kinase)
    
    if df_all.empty:
        print(f"[ERROR] No docking data found for approved receptors")
        exit(1)
    
    print(f"[INFO] {kinase.upper()} — Total records: {len(df_all)}, "
          f"Receptors: {df_all['Receptor'].nunique()}, "
          f"Ligands: {df_all['Ligand'].nunique()}, "
          f"Actives: {(df_all['Label']=='active').sum()}, "
          f"Decoys: {(df_all['Label']=='decoy').sum()}")

    # Save all collected data
    df_all.to_csv(os.path.join(output_dir, "all_scores_detailed.csv"), index=False)

    if calculate_per_receptor:
        # NEW: Calculate per-receptor EF
        per_receptor_ef_df = calculate_per_receptor_ef(df_all, top_percentages=[1, 5, 10])
        
        # Analyze and visualize results
        top_receptors, summary_stats = analyze_per_receptor_ef(per_receptor_ef_df, output_dir, top_n)
        
        # Compare with ensemble methods
        comparison_df, comparison_stats = compare_individual_vs_ensemble_ef(df_all, per_receptor_ef_df, output_dir)
        
        print(f"\n[INFO] Per-receptor EF analysis complete!")
        print(f"[INFO] Results saved to: {output_dir}")
        print(f"[INFO] Best individual receptor EF@1%: {per_receptor_ef_df['EF@1%'].max():.3f}")
        print(f"[INFO] Mean individual receptor EF@1%: {per_receptor_ef_df['EF@1%'].mean():.3f}")
        
    else:
        # Original ensemble analysis
        print(f"[INFO] Use --per_receptor flag to calculate EF for individual receptors")
        print(f"[INFO] Running basic ensemble analysis...")
        
        # Calculate ensemble EF for comparison
        df_ensemble = compute_scores_with_method(df_all, method='mean')
        
        ef_records = []
        for pct in [1, 5, 10]:
            ef = calculate_ef(df_ensemble, pct)
            ef_records.append({"Top%": pct, "EF": ef, "Method": "ensemble_mean"})
            print(f"Ensemble EF@{pct}% (mean method) = {ef:.3f}")
        
        pd.DataFrame(ef_records).to_csv(os.path.join(output_dir, "ensemble_ef_summary.csv"), index=False)