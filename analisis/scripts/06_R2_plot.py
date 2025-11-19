#!/usr/bin/env python3
"""
SIMPLIFIED: All-in-One Conformational State R² Analysis
Combines three steps into one simple script:
1. Merge conformational states with docking results
2. Calculate performance metrics (EF1%, AUC)
3. Calculate R² per kinase

Usage:
    python3 R2_analysis.py
"""

import pandas as pd
import numpy as np
from sklearn.metrics import r2_score, roc_auc_score
from scipy import stats
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION - Edit these file paths
# ============================================================================
DOCKING_FILE = '../data/output.csv'
CONFORMATIONAL_STATE_FILE = '../data/conf_state_structure_mapping.csv'
OUTPUT_PREFIX = 'results'

# ============================================================================
# STEP 1: MERGE DATA
# ============================================================================
def load_and_merge_data():
    """Load docking results and conformational states, then merge"""
    print("="*80)
    print("STEP 1: LOADING AND MERGING DATA")
    print("="*80)
    
    # Load files
    print(f"\nLoading docking results: {DOCKING_FILE}")
    docking_df = pd.read_csv(DOCKING_FILE)
    print(f"  Loaded {len(docking_df):,} docking records")
    print(f"  Kinases: {docking_df['kinase'].nunique()}")
    print(f"  Structures: {docking_df['structure_id'].nunique()}")
    
    print(f"\nLoading conformational states: {CONFORMATIONAL_STATE_FILE}")
    conf_df = pd.read_csv(CONFORMATIONAL_STATE_FILE)
    print(f"  Loaded {len(conf_df):,} structure mappings")
    print(f"  Kinases: {conf_df['kinase'].nunique()}")
    
    # VALIDATION: Check if all structures have the same number of compounds per kinase
    print("\n" + "="*80)
    print("DATA VALIDATION: Checking compound counts per structure")
    print("="*80)
    
    validation_issues = []
    
    for kinase in sorted(docking_df['kinase'].unique()):
        kinase_data = docking_df[docking_df['kinase'] == kinase]
        
        # Count compounds per structure
        compounds_per_structure = kinase_data.groupby('structure_id')['compound_id'].nunique()
        
        # Check if all structures have the same number of compounds
        if compounds_per_structure.nunique() > 1:
            min_compounds = compounds_per_structure.min()
            max_compounds = compounds_per_structure.max()
            
            # Find structures with incomplete data
            incomplete_structures = compounds_per_structure[compounds_per_structure < max_compounds]
            
            validation_issues.append({
                'kinase': kinase,
                'expected_compounds': max_compounds,
                'incomplete_structures': len(incomplete_structures),
                'incomplete_ids': incomplete_structures.index.tolist()
            })
            
            print(f"\n⚠ WARNING: {kinase}")
            print(f"  Expected compounds per structure: {max_compounds}")
            print(f"  Incomplete structures found: {len(incomplete_structures)}")
            print(f"  Structure IDs with incomplete data: {incomplete_structures.index.tolist()[:10]}" + 
                  ("..." if len(incomplete_structures) > 10 else ""))
            
            # Show compound counts for incomplete structures
            for struct_id, count in incomplete_structures.head(5).items():
                print(f"    Structure {struct_id}: {count} compounds (missing {max_compounds - count})")
    
    if not validation_issues:
        print("\n✓ All structures have consistent compound counts")
    else:
        print(f"\n{'='*80}")
        print(f"SUMMARY: Found {len(validation_issues)} kinase(s) with incomplete docking data")
        print(f"{'='*80}")
        print("NOTE: These structures will be automatically excluded from analysis")
        print("      (structures with n_actives = 0 or n_decoys = 0)")
    
    # Merge
    print("\n" + "="*80)
    print("MERGING DATASETS")
    print("="*80)
    merged_df = docking_df.merge(conf_df[['kinase', 'structure_id', 'conformational_state']], 
                                  on=['kinase', 'structure_id'], 
                                  how='left')
    
    print(f"  Merged: {len(merged_df):,} records")
    print(f"  With conformational state: {merged_df['conformational_state'].notna().sum():,}")
    print(f"  Missing conformational state: {merged_df['conformational_state'].isna().sum():,}")
    
    return merged_df, validation_issues

# ============================================================================
# STEP 2: CALCULATE PERFORMANCE METRICS
# ============================================================================
def calculate_enrichment_factor(scores, labels, percentage=1):
    """Calculate enrichment factor at given percentage"""
    if len(scores) == 0:
        return 0.0
    
    df = pd.DataFrame({'score': scores, 'label': labels})
    df_sorted = df.sort_values('score')  # Lower docking score is better
    
    n_top = max(1, int(len(df_sorted) * percentage / 100))
    top_compounds = df_sorted.head(n_top)
    
    actives_in_top = (top_compounds['label'] == 'active').sum()
    total_actives = (df['label'] == 'active').sum()
    
    if total_actives == 0:
        return 0.0
    
    expected_actives = total_actives * percentage / 100
    if expected_actives == 0:
        return 0.0
    
    ef = actives_in_top / expected_actives
    return ef

def calculate_auc(scores, labels):
    """Calculate AUC"""
    if len(set(labels)) < 2:
        return 0.5
    
    y_true = [1 if label == 'active' else 0 for label in labels]
    y_scores = [-score for score in scores]  # Invert because lower is better
    
    try:
        auc = roc_auc_score(y_true, y_scores)
        return auc
    except:
        return 0.5

def calculate_structure_performance(merged_df):
    """Calculate EF and AUC for each structure"""
    print("\n" + "="*80)
    print("STEP 2: CALCULATING PERFORMANCE METRICS")
    print("="*80)
    
    performance_results = []
    grouped = merged_df.groupby(['kinase', 'structure_id'])
    
    total = len(grouped)
    print(f"\nProcessing {total} structures...")
    
    for i, ((kinase, structure_id), group) in enumerate(grouped, 1):
        if i % 100 == 0 or i == total:
            print(f"  Progress: {i}/{total}")
        
        conf_state = group['conformational_state'].iloc[0]
        scores = group['docking_score'].values
        labels = group['compound_type'].values
        
        # Calculate metrics
        ef1 = calculate_enrichment_factor(scores, labels, 1)
        ef5 = calculate_enrichment_factor(scores, labels, 5)
        ef10 = calculate_enrichment_factor(scores, labels, 10)
        auc = calculate_auc(scores, labels)
        
        # Count compounds
        n_actives = (labels == 'active').sum()
        n_decoys = (labels == 'decoy').sum()
        
        performance_results.append({
            'kinase': kinase,
            'structure_id': structure_id,
            'conformational_state': conf_state,
            'ef1': ef1,
            'ef5': ef5,
            'ef10': ef10,
            'auc': auc,
            'n_actives': n_actives,
            'n_decoys': n_decoys,
            'n_total': len(labels)
        })
    
    performance_df = pd.DataFrame(performance_results)
    
    # FILTER OUT INCOMPLETE STRUCTURES
    print(f"\nInitial structures: {len(performance_df)}")
    
    incomplete = performance_df[(performance_df['n_actives'] == 0) | (performance_df['n_decoys'] == 0)]
    if len(incomplete) > 0:
        print(f"\n⚠ Found {len(incomplete)} structures with incomplete data:")
        for kinase in incomplete['kinase'].unique():
            kinase_incomplete = incomplete[incomplete['kinase'] == kinase]
            print(f"  {kinase}: structures {kinase_incomplete['structure_id'].tolist()}")
        
        performance_df = performance_df[(performance_df['n_actives'] > 0) & (performance_df['n_decoys'] > 0)]
        print(f"\n✓ Filtered to {len(performance_df)} complete structures")
    else:
        print(f"\n✓ All structures have complete data")
    
    print(f"Final: {len(performance_df)} structures analyzed")
    
    return performance_df

# ============================================================================
# STEP 3: CALCULATE R² PER KINASE
# ============================================================================
def calculate_r2_by_kinase(performance_df):
    """Calculate R² for each kinase"""
    print("\n" + "="*80)
    print("STEP 3: CALCULATING R² BY KINASE")
    print("="*80)
    
    # Remove structures without conformational state
    df = performance_df.dropna(subset=['conformational_state']).copy()
    print(f"\nStructures with conformational state data: {len(df)}")
    
    results = []
    
    for kinase in sorted(df['kinase'].unique()):
        kinase_data = df[df['kinase'] == kinase].copy()
        
        # Need at least 2 states and 5 structures
        n_states = kinase_data['conformational_state'].nunique()
        n_structures = len(kinase_data)
        
        if n_states < 2 or n_structures < 5:
            print(f"  {kinase}: Skipped (n_structures={n_structures}, n_states={n_states})")
            continue
        
        # Calculate R²
        state_means = kinase_data.groupby('conformational_state')['ef1'].mean()
        kinase_data['predicted_ef1'] = kinase_data['conformational_state'].map(state_means)
        r2 = r2_score(kinase_data['ef1'], kinase_data['predicted_ef1'])
        
        # ANOVA F-test
        groups = [group['ef1'].values for _, group in kinase_data.groupby('conformational_state')]
        f_stat, p_value = stats.f_oneway(*groups)
        
        results.append({
            'kinase': kinase,
            'r2': r2,
            'var_explained_pct': r2 * 100,
            'n_structures': n_structures,
            'n_states': n_states,
            'f_stat': f_stat,
            'p_value': p_value
        })
        
        print(f"  {kinase}: R²={r2:.5f}, n_structures={n_structures}, n_states={n_states}")
    
    return pd.DataFrame(results)

# ============================================================================
# VISUALIZATION AND SUMMARY
# ============================================================================
def plot_r2_by_kinase(results_df, output_file):
    """Create publication-quality bar plot of R² by kinase"""
    
    if len(results_df) == 0:
        print("\nNo R² results to plot!")
        return
    
    # Sort by R²
    results_sorted = results_df.sort_values('r2', ascending=True)
    
    # Publication settings
    plt.rcParams['font.family'] = 'Arial'  # or 'Helvetica'
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.linewidth'] = 1.5
    
    # Create figure
    fig, ax = plt.subplots(figsize=(6, max(4, len(results_sorted) * 0.25)))  # 더 compact
    
    # Color mapping
    colors = []
    for r2_val in results_sorted['r2']:
        if r2_val > 0.3:
            colors.append('#2ca02c')  # green
        elif r2_val > 0.1:
            colors.append('#ff7f0e')  # orange
        else:
            colors.append('#d62728')  # red
    
    # Create bars
    bars = ax.barh(results_sorted['kinase'], results_sorted['r2'], 
                   color=colors, edgecolor='black', linewidth=0.5)
    
    # Labels - cleaner for publication
    ax.set_xlabel('Coefficient of Determination (R²)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Kinase', fontsize=11, fontweight='bold')
    ax.set_title('Conformational State Predictive Power', 
                 fontsize=12, fontweight='bold', pad=10)
    
    # Grid and limits
    ax.grid(True, alpha=0.3, axis='x', linestyle='--', linewidth=0.5)
    ax.set_xlim(0, 0.41)
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4])
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#2ca02c', edgecolor='black', label='Highest (R² > 0.3)'),
        Patch(facecolor='#ff7f0e', edgecolor='black', label='Moderate (0.1 < R² ≤ 0.3)'),
        Patch(facecolor='#d62728', edgecolor='black', label='Low (R² ≤ 0.1)')
    ]
    ax.legend(handles=legend_elements, loc='lower right', frameon=True, 
              fontsize=9, edgecolor='black')
    
    # Remove top and right spines (cleaner look)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    # Save both formats
    plt.savefig(f'../data/{output_file}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'../data/{output_file}.pdf', dpi=300, bbox_inches='tight', 
                transparent=True)  # transparent background for PDF
    
    print(f"\n✓ Publication figure saved:")
    print(f"  PNG: ../data/{output_file}.png")
    print(f"  PDF: ../data/{output_file}.pdf")
    
    plt.close()

def print_summary(results_df):
    """Print detailed summary"""
    
    if len(results_df) == 0:
        print("\nNo results to summarize!")
        return
    
    print("\n" + "="*80)
    print("FINAL RESULTS: R² BY KINASE")
    print("="*80)
    
    # Sort by R²
    results_sorted = results_df.sort_values('r2', ascending=False)
    
    print(f"\n{'Kinase':<10} {'R²':<8} {'Var%':<8} {'N':<6} {'States':<8} {'p-value':<12} {'Sig':<5}")
    print("─"*70)
    
    for _, row in results_sorted.iterrows():
        sig = "***" if row['p_value'] < 0.001 else "**" if row['p_value'] < 0.01 else "*" if row['p_value'] < 0.05 else "ns"
        print(f"{row['kinase']:<10} {row['r2']:<8.5f} {row['var_explained_pct']:<8.1f} "
              f"{row['n_structures']:<6} {row['n_states']:<8} {row['p_value']:<12.2e} {sig:<5}")
    
    # Category summary
    high = results_df[results_df['r2'] > 0.3]
    medium = results_df[(results_df['r2'] > 0.1) & (results_df['r2'] <= 0.3)]
    low = results_df[results_df['r2'] <= 0.1]
    
    print("\n" + "="*80)
    print("SUMMARY BY CATEGORY")
    print("="*80)
    print(f"HIGH (R² > 0.3):   {len(high):2d} kinases ({len(high)/len(results_df)*100:5.1f}%)")
    if len(high) > 0:
        print(f"  Kinases: {', '.join(high['kinase'].tolist())}")
    
    print(f"\nMEDIUM (0.1-0.3):  {len(medium):2d} kinases ({len(medium)/len(results_df)*100:5.1f}%)")
    if len(medium) > 0:
        print(f"  Kinases: {', '.join(medium['kinase'].tolist())}")
    
    print(f"\nLOW (R² < 0.1):    {len(low):2d} kinases ({len(low)/len(results_df)*100:5.1f}%)")
    if len(low) > 0:
        print(f"  Kinases: {', '.join(low['kinase'].tolist())}")
    
    print(f"\nOverall Statistics:")
    print(f"  Mean R²:   {results_df['r2'].mean():.5f}")
    print(f"  Median R²: {results_df['r2'].median():.5f}")
    print(f"  Min R²:    {results_df['r2'].min():.5f} ({results_df.loc[results_df['r2'].idxmin(), 'kinase']})")
    print(f"  Max R²:    {results_df['r2'].max():.5f} ({results_df.loc[results_df['r2'].idxmax(), 'kinase']})")

def show_top_states_per_kinase(performance_df, results_df):
    """Show best performing conformational states for top kinases"""
    
    if len(results_df) == 0:
        return
    
    print("\n" + "="*80)
    print("BEST CONFORMATIONAL STATES BY KINASE (Top R² kinases)")
    print("="*80)
    
    # Get top 5 kinases by R²
    top_kinases = results_df.nlargest(5, 'r2')['kinase'].tolist()
    
    df = performance_df.dropna(subset=['conformational_state'])
    
    for kinase in top_kinases:
        kinase_data = df[df['kinase'] == kinase]
        
        if len(kinase_data) == 0:
            continue
        
        # Calculate mean EF1% by conformational state
        state_performance = kinase_data.groupby('conformational_state').agg({
            'ef1': ['mean', 'count']
        })
        
        state_performance.columns = ['mean_ef1', 'n_structures']
        state_performance = state_performance.sort_values('mean_ef1', ascending=False)
        
        r2_val = results_df[results_df['kinase'] == kinase]['r2'].iloc[0]
        
        print(f"\n{kinase} (R² = {r2_val:.5f}):")
        print(f"  {'State':<25} {'Mean EF1%':<12} {'N Structures'}")
        print("  " + "─"*50)
        
        for state, row in state_performance.head(5).iterrows():
            print(f"  {state:<25} {row['mean_ef1']:<12.5f} {int(row['n_structures'])}")

# ============================================================================
# MAIN FUNCTION
# ============================================================================
def main():
    """Run complete analysis pipeline"""
    
    print("\n" + "="*80)
    print("SIMPLIFIED CONFORMATIONAL STATE R² ANALYSIS")
    print("="*80)
    
    try:
        # Step 1: Load and merge data
        merged_df, validation_issues = load_and_merge_data()
        
        # Step 2: Calculate performance metrics
        performance_df = calculate_structure_performance(merged_df)
        
        # Save performance data with 5 decimal places
        performance_file = f'../data/{OUTPUT_PREFIX}_performance.csv'
        performance_df.to_csv(performance_file, index=False, float_format='%.5f')
        print(f"\nPerformance data saved: {performance_file}")
        
        # Step 3: Calculate R² by kinase
        results_df = calculate_r2_by_kinase(performance_df)
        
        if len(results_df) > 0:
            # Save R² results with 5 decimal places
            r2_file = f'../data/{OUTPUT_PREFIX}_r2.csv'
            results_df.to_csv(r2_file, index=False, float_format='%.5f')
            print(f"R² results saved: {r2_file}")
            
            # Create visualization
            plot_file = f'../data/{OUTPUT_PREFIX}_r2_plot'
            plot_r2_by_kinase(results_df, plot_file)
            
            # Print summaries
            print_summary(results_df)
            show_top_states_per_kinase(performance_df, results_df)
            
            print("\n" + "="*80)
            print("ANALYSIS COMPLETE!")
            print("="*80)
            print(f"\nOutput files:")
            print(f"  1. {performance_file}")
            print(f"  2. {r2_file}")
            print(f"  3. {plot_file}.png")
            print(f"  4. {plot_file}.pdf")
            
            if validation_issues:
                print(f"\n⚠ NOTE: {len(validation_issues)} kinase(s) had incomplete docking data")
                print("  Incomplete structures were automatically excluded from analysis")
                print("  Run docking to completion for complete results")
        else:
            print("\nWARNING: No valid R² results generated!")
            print("This might happen if kinases don't have enough structures or conformational states.")
        
    except FileNotFoundError as e:
        print(f"\nERROR: File not found!")
        print(f"  {e}")
        print(f"\nMake sure these files exist:")
        print(f"  - {DOCKING_FILE}")
        print(f"  - {CONFORMATIONAL_STATE_FILE}")
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
