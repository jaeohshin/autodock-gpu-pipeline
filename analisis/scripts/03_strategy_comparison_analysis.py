#!/usr/bin/env python3
"""
Kinase Virtual Screening Analysis Data Exporter
===============================================

This script processes raw docking results and exports clean, organized 
analysis data files for publication-quality figure generation.

Outputs:
1. kinase_metrics_summary.csv - Main metrics for each kinase-strategy combination
2. strategy_comparison_data.csv - Data formatted specifically for box plots
3. kinase_detailed_stats.csv - Additional statistics and metadata

Author: Jaeoh Shin
Date: November 12, 2025
"""

import os
import sys
import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, auc
from datetime import datetime

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input file
CSV_FILE = "../data/output.csv"
# Output directory
OUTPUT_DIR = "../data"

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# LOAD DATA
# ============================================================================

def load_and_validate_data(csv_file):
    """Load and validate the docking results data"""
    print(f"Loading docking results from: {csv_file}")
    
    if not os.path.exists(csv_file):
        print(f"ERROR: CSV file '{csv_file}' not found!")
        csv_files = [f for f in os.listdir('.') if f.endswith('.csv')]
        print("\nAvailable CSV files:")
        for f in csv_files:
            print(f"  - {f}")
        sys.exit(1)
    
    df = pd.read_csv(csv_file)
    
    # Validate required columns
    required_columns = ['kinase', 'compound_id', 'compound_type', 'docking_score', 'structure_id']
    missing_columns = [col for col in required_columns if col not in df.columns]
    
    if missing_columns:
        print(f"ERROR: Missing required columns: {missing_columns}")
        print(f"Available columns: {list(df.columns)}")
        sys.exit(1)
    
    print(f"✓ Loaded {len(df):,} docking results")
    print(f"✓ Kinases: {df['kinase'].nunique()}")
    print(f"✓ Compounds: {df['compound_id'].nunique()}")
    print(f"✓ Structures: {df['structure_id'].nunique()}")
    print(f"✓ Actives: {(df['compound_type'] == 'active').sum():,}")
    print(f"✓ Decoys: {(df['compound_type'] == 'decoy').sum():,}")
    
    return df

# ============================================================================
# ANALYSIS FUNCTIONS
# ============================================================================

def calculate_roc_metrics(actives, decoys):
    """Calculate ROC AUC and related metrics"""
    if len(actives) == 0 or len(decoys) == 0:
        return np.nan, np.nan, np.nan
    
    try:
        y_true = np.concatenate([np.ones(len(actives)), np.zeros(len(decoys))])
        y_scores = -np.concatenate([actives, decoys])  # Negative because lower scores are better
        
        fpr, tpr, thresholds = roc_curve(y_true, y_scores)
        roc_auc = auc(fpr, tpr)
        
        return roc_auc, fpr, tpr
    except Exception as e:
        print(f"Warning: Error calculating ROC metrics: {e}")
        return np.nan, np.nan, np.nan

def calculate_enrichment_factors(compound_scores, percentages=[1, 5, 10]):
    """Calculate enrichment factors at specified percentages"""
    n_actives = (compound_scores['compound_type'] == 'active').sum()
    n_total = len(compound_scores)
    
    if n_actives == 0 or n_total == 0:
        return {f'EF{p}%': np.nan for p in percentages}
    
    # Sort by docking score (lower is better)
    sorted_compounds = compound_scores.sort_values('docking_score')
    
    ef_results = {}
    for pct in percentages:
        top_n = max(1, int(pct/100 * n_total))
        actives_in_top = (sorted_compounds.head(top_n)['compound_type'] == 'active').sum()
        
        # EF = (actives_in_top / top_n) / (total_actives / total_compounds)
        expected_rate = n_actives / n_total
        observed_rate = actives_in_top / top_n
        ef = observed_rate / expected_rate if expected_rate > 0 else 0
        
        ef_results[f'EF{pct}%'] = ef
    
    return ef_results

def calculate_strategy_metrics(df, strategy='best_score'):
    """Calculate metrics for a specific strategy"""
    print(f"\nCalculating metrics for strategy: {strategy}")
    
    results = []
    
    for kinase in sorted(df['kinase'].unique()):
        kinase_data = df[df['kinase'] == kinase]
        
        # Apply strategy to get compound-level scores
        if strategy == 'best_score':
            # Select best (lowest) scoring conformation per compound
            compound_scores = kinase_data.groupby(['compound_id', 'compound_type'])['docking_score'].min().reset_index()
        elif strategy == 'consensus':
            # Use average score across conformations
            compound_scores = kinase_data.groupby(['compound_id', 'compound_type'])['docking_score'].mean().reset_index()
        else:
            raise ValueError(f"Unknown strategy: {strategy}")
        
        # Get actives and decoys
        actives = compound_scores[compound_scores['compound_type'] == 'active']['docking_score'].values
        decoys = compound_scores[compound_scores['compound_type'] == 'decoy']['docking_score'].values
        
        # Calculate ROC AUC
        roc_auc, fpr, tpr = calculate_roc_metrics(actives, decoys)
        
        # Calculate enrichment factors
        ef_metrics = calculate_enrichment_factors(compound_scores, [1, 5, 10])
        
        # Calculate additional statistics
        active_stats = {
            'active_mean_score': np.mean(actives) if len(actives) > 0 else np.nan,
            'active_median_score': np.median(actives) if len(actives) > 0 else np.nan,
            'active_std_score': np.std(actives) if len(actives) > 0 else np.nan,
        }
        
        decoy_stats = {
            'decoy_mean_score': np.mean(decoys) if len(decoys) > 0 else np.nan,
            'decoy_median_score': np.median(decoys) if len(decoys) > 0 else np.nan,
            'decoy_std_score': np.std(decoys) if len(decoys) > 0 else np.nan,
        }
        
        # Score separation metrics
        score_separation = {
            'mean_score_difference': active_stats['active_mean_score'] - decoy_stats['decoy_mean_score'],
            'median_score_difference': active_stats['active_median_score'] - decoy_stats['decoy_median_score'],
        }
        
        # Compile results
        result = {
            'kinase': kinase,
            'strategy': strategy,
            'AUC': roc_auc,
            **ef_metrics,
            **active_stats,
            **decoy_stats,
            **score_separation,
            'n_actives': len(actives),
            'n_decoys': len(decoys),
            'n_compounds': len(compound_scores),
            'n_structures': kinase_data['structure_id'].nunique(),
            'n_conformations': len(kinase_data),
        }
        
        results.append(result)
        
        # Progress indicator
        if len(results) % 10 == 0:
            print(f"  Processed {len(results)} kinases...")
    
    return pd.DataFrame(results)

# ============================================================================
# DATA EXPORT FUNCTIONS
# ============================================================================

def export_main_metrics(all_metrics_df, output_dir):
    """Export main metrics summary"""
    filename = os.path.join(output_dir, "kinase_metrics_summary.csv")
    
    # Select key columns for the main summary
    key_columns = [
        'kinase', 'strategy', 'AUC', 'EF1%', 'EF5%', 'EF10%',
        'n_actives', 'n_decoys', 'n_compounds', 'n_structures'
    ]
    
    summary_df = all_metrics_df[key_columns].copy()
    summary_df.to_csv(filename, index=False, float_format='%.4f')
    
    print(f"✓ Exported main metrics to: {filename}")
    return filename

def export_strategy_comparison_data(all_metrics_df, output_dir):
    """Export data formatted specifically for strategy comparison box plots"""
    filename = os.path.join(output_dir, "strategy_comparison_data.csv")
    
    # Reshape data for plotting (long format)
    plot_data = []
    
    for strategy in all_metrics_df['strategy'].unique():
        strategy_data = all_metrics_df[all_metrics_df['strategy'] == strategy]
        
        for _, row in strategy_data.iterrows():
            # Add AUC data point
            plot_data.append({
                'kinase': row['kinase'],
                'strategy': strategy,
                'metric': 'AUC',
                'value': row['AUC']
            })
            
            # Add EF data points
            for ef_col in ['EF1%', 'EF5%', 'EF10%']:
                plot_data.append({
                    'kinase': row['kinase'],
                    'strategy': strategy,
                    'metric': ef_col,
                    'value': row[ef_col]
                })
    
    plot_df = pd.DataFrame(plot_data)
    plot_df.to_csv(filename, index=False, float_format='%.4f')
    
    print(f"✓ Exported strategy comparison data to: {filename}")
    return filename

def export_detailed_stats(all_metrics_df, output_dir):
    """Export detailed statistics and metadata"""
    filename = os.path.join(output_dir, "kinase_detailed_stats.csv")
    
    all_metrics_df.to_csv(filename, index=False, float_format='%.4f')
    
    print(f"✓ Exported detailed statistics to: {filename}")
    return filename

def export_summary_statistics(all_metrics_df, output_dir):
    """Export overall summary statistics"""
    filename = os.path.join(output_dir, "analysis_summary.txt")
    
    with open(filename, 'w') as f:
        f.write("Kinase Virtual Screening Analysis Summary\n")
        f.write("=" * 50 + "\n")
        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Overall statistics
        f.write("Dataset Overview:\n")
        f.write(f"  Total kinases analyzed: {all_metrics_df['kinase'].nunique()}\n")
        f.write(f"  Strategies compared: {list(all_metrics_df['strategy'].unique())}\n\n")
        
        # Strategy-wise statistics
        for strategy in all_metrics_df['strategy'].unique():
            strategy_data = all_metrics_df[all_metrics_df['strategy'] == strategy]
            
            f.write(f"Strategy: {strategy}\n")
            f.write("-" * 20 + "\n")
            f.write(f"  Mean AUC: {strategy_data['AUC'].mean():.4f} ± {strategy_data['AUC'].std():.4f}\n")
            f.write(f"  Median AUC: {strategy_data['AUC'].median():.4f}\n")
            f.write(f"  AUC Range: {strategy_data['AUC'].min():.4f} - {strategy_data['AUC'].max():.4f}\n")
            
            for ef_col in ['EF1%', 'EF5%', 'EF10%']:
                f.write(f"  Mean {ef_col}: {strategy_data[ef_col].mean():.2f} ± {strategy_data[ef_col].std():.2f}\n")
            
            f.write(f"  Kinases with AUC > 0.7: {(strategy_data['AUC'] > 0.7).sum()}\n")
            f.write(f"  Kinases with AUC > 0.8: {(strategy_data['AUC'] > 0.8).sum()}\n")
            f.write("\n")
    
    print(f"✓ Exported analysis summary to: {filename}")
    return filename

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main analysis and export pipeline"""
    print("Kinase Virtual Screening Analysis Data Exporter")
    print("=" * 50)
    
    # Load data
    df = load_and_validate_data(CSV_FILE)
    
    # Calculate metrics for both strategies
    strategies = ['best_score', 'consensus']
    all_metrics = []
    
    for strategy in strategies:
        strategy_metrics = calculate_strategy_metrics(df, strategy)
        all_metrics.append(strategy_metrics)
    
    # Combine all metrics
    combined_metrics = pd.concat(all_metrics, ignore_index=True)
    
    print(f"\n✓ Analysis complete! Calculated metrics for {len(combined_metrics)} kinase-strategy combinations")
    
    # Export all data files
    print(f"\nExporting analysis results to: {OUTPUT_DIR}/")
    print("-" * 30)
    
    export_main_metrics(combined_metrics, OUTPUT_DIR)
    export_strategy_comparison_data(combined_metrics, OUTPUT_DIR)
    export_detailed_stats(combined_metrics, OUTPUT_DIR)
    export_summary_statistics(combined_metrics, OUTPUT_DIR)
    
    print("\n" + "=" * 50)
    print("Analysis data export completed successfully!")
    print(f"\nOutput files created in '{OUTPUT_DIR}/':")
    print("  1. kinase_metrics_summary.csv - Main metrics for plotting")
    print("  2. strategy_comparison_data.csv - Data for box plots") 
    print("  3. kinase_detailed_stats.csv - Complete statistics")
    print("  4. analysis_summary.txt - Human-readable summary")
    print("\nReady for Step 2: Publication figure generation!")

if __name__ == "__main__":
    main()
