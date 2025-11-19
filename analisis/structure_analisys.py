#!/usr/bin/env python3
"""
Structure-Performance Correlation Analysis
Merges docking results with conformational state data and analyzes:
1. Performance differences between conformational states
2. Within-state performance variation
3. Individual structure performance metrics
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

def load_and_merge_data(mapping_file, docking_file):
    """Load and merge conformational state mapping with docking results"""
    print("Loading data files...")
    
    # Load mapping file
    mapping_df = pd.read_csv(mapping_file)
    print(f"Mapping file shape: {mapping_df.shape}")
    print(f"Mapping columns: {mapping_df.columns.tolist()}")
    print(f"Unique kinases in mapping: {mapping_df['kinase'].nunique()}")
    print(f"Total structures in mapping: {len(mapping_df)}")
    
    # Load docking results
    docking_df = pd.read_csv(docking_file)
    print(f"\nDocking file shape: {docking_df.shape}")
    print(f"Docking columns: {docking_df.columns.tolist()}")
    print(f"Unique kinases in docking: {docking_df['kinase'].nunique()}")
    print(f"Unique structures in docking: {docking_df['structure_id'].nunique()}")
    
    # Check data types
    print(f"\nData types - Mapping structure_id: {mapping_df['structure_id'].dtype}")
    print(f"Data types - Docking structure_id: {docking_df['structure_id'].dtype}")
    
    # Merge datasets
    print("\nMerging datasets...")
    merged_df = docking_df.merge(mapping_df, on=['kinase', 'structure_id'], how='left')
    
    print(f"Merged shape: {merged_df.shape}")
    print(f"Missing conformational states: {merged_df['conformational_state'].isna().sum()}")
    
    # Show conformational state distribution
    print(f"\nConformational state distribution:")
    state_counts = merged_df['conformational_state'].value_counts()
    print(state_counts)
    
    return merged_df, mapping_df, docking_df

def calculate_enrichment_factor(scores, labels, percentage=1):
    """Calculate enrichment factor at given percentage"""
    if len(scores) == 0:
        return 0.0
    
    # Create dataframe and sort by score (lower is better for docking)
    df = pd.DataFrame({'score': scores, 'label': labels})
    df_sorted = df.sort_values('score')
    
    # Calculate number of compounds in top percentage
    n_top = max(1, int(len(df_sorted) * percentage / 100))
    top_compounds = df_sorted.head(n_top)
    
    # Count actives in top compounds
    actives_in_top = (top_compounds['label'] == 'active').sum()
    total_actives = (df['label'] == 'active').sum()
    
    if total_actives == 0:
        return 0.0
    
    # Calculate enrichment factor
    expected_actives = total_actives * percentage / 100
    if expected_actives == 0:
        return 0.0
    
    ef = actives_in_top / expected_actives
    return ef

def calculate_auc(scores, labels):
    """Calculate AUC (Area Under Curve)"""
    if len(set(labels)) < 2:  # Need both actives and decoys
        return 0.5
    
    # Convert labels to binary (1 for active, 0 for decoy)
    y_true = [1 if label == 'active' else 0 for label in labels]
    
    # For docking scores, lower is better, so we need to invert
    y_scores = [-score for score in scores]
    
    try:
        auc = roc_auc_score(y_true, y_scores)
        return auc
    except:
        return 0.5

def calculate_structure_performance(merged_df):
    """Calculate EF and AUC for each individual structure"""
    print("\nCalculating performance metrics for each structure...")
    
    performance_results = []
    
    # Group by kinase and structure_id
    grouped = merged_df.groupby(['kinase', 'structure_id'])
    
    total_structures = len(grouped)
    processed = 0
    
    for (kinase, structure_id), group in grouped:
        processed += 1
        if processed % 100 == 0:
            print(f"Processing structure {processed}/{total_structures}")
        
        # Get conformational state (should be same for all compounds in this structure)
        conf_state = group['conformational_state'].iloc[0]
        
        # Get scores and labels
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
        n_total = len(labels)
        
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
            'n_total': n_total
        })
    
    performance_df = pd.DataFrame(performance_results)
    print(f"Calculated performance for {len(performance_df)} structures")
    
    return performance_df

def analyze_conformational_states(performance_df):
    """Analyze performance differences between conformational states"""
    print("\n=== CONFORMATIONAL STATE ANALYSIS ===")
    
    # Overall statistics by conformational state
    state_stats = performance_df.groupby('conformational_state').agg({
        'ef1': ['count', 'mean', 'std', 'min', 'max'],
        'ef5': ['mean', 'std'],
        'ef10': ['mean', 'std'],
        'auc': ['mean', 'std', 'min', 'max']
    }).round(3)
    
    print("\nPerformance by Conformational State:")
    print(state_stats)
    
    # Show top performing conformational states
    state_means = performance_df.groupby('conformational_state')[['ef1', 'ef5', 'ef10', 'auc']].mean().round(3)
    print(f"\nTop conformational states by EF1%:")
    top_ef1 = state_means.sort_values('ef1', ascending=False)
    print(top_ef1)
    
    return state_stats, state_means

def analyze_within_state_variation(performance_df):
    """Analyze performance variation within the same conformational state"""
    print("\n=== WITHIN-STATE VARIATION ANALYSIS ===")
    
    within_state_results = {}
    
    for state in performance_df['conformational_state'].unique():
        if pd.isna(state):
            continue
            
        state_data = performance_df[performance_df['conformational_state'] == state]
        
        if len(state_data) < 2:  # Need at least 2 structures to see variation
            continue
        
        # Calculate variation metrics
        ef1_range = state_data['ef1'].max() - state_data['ef1'].min()
        ef1_std = state_data['ef1'].std()
        auc_range = state_data['auc'].max() - state_data['auc'].min()
        auc_std = state_data['auc'].std()
        
        within_state_results[state] = {
            'n_structures': len(state_data),
            'ef1_mean': state_data['ef1'].mean(),
            'ef1_std': ef1_std,
            'ef1_range': ef1_range,
            'ef1_min': state_data['ef1'].min(),
            'ef1_max': state_data['ef1'].max(),
            'auc_mean': state_data['auc'].mean(),
            'auc_std': auc_std,
            'auc_range': auc_range,
            'auc_min': state_data['auc'].min(),
            'auc_max': state_data['auc'].max()
        }
        
        # Show examples of best and worst performers within this state
        best_structure = state_data.loc[state_data['ef1'].idxmax()]
        worst_structure = state_data.loc[state_data['ef1'].idxmin()]
        
        print(f"\nState: {state} ({len(state_data)} structures)")
        print(f"  EF1% range: {ef1_range:.2f} (std: {ef1_std:.2f})")
        print(f"  Best performer: {best_structure['kinase']}-{best_structure['structure_id']} (EF1%: {best_structure['ef1']:.2f})")
        print(f"  Worst performer: {worst_structure['kinase']}-{worst_structure['structure_id']} (EF1%: {worst_structure['ef1']:.2f})")
    
    return within_state_results

def identify_top_performers(performance_df, top_n=10):
    """Identify top performing individual structures"""
    print(f"\n=== TOP {top_n} PERFORMING STRUCTURES ===")
    
    # Sort by EF1% (primary) and AUC (secondary)
    top_structures = performance_df.nlargest(top_n, ['ef1', 'auc'])
    
    print("Top performers by EF1%:")
    for idx, row in top_structures.iterrows():
        print(f"  {row['kinase']}-{row['structure_id']:02d} ({row['conformational_state']}): "
              f"EF1%={row['ef1']:.2f}, AUC={row['auc']:.3f}")
    
    # Show conformational state distribution of top performers
    top_states = top_structures['conformational_state'].value_counts()
    print(f"\nConformational states among top {top_n} performers:")
    print(top_states)
    
    return top_structures

def create_visualizations(performance_df):
    """Create visualizations for the analysis"""
    print("\nCreating visualizations...")
    
    # Set up plotting style
    plt.style.use('default')
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Filter out missing conformational states
    states_with_data = performance_df.dropna(subset=['conformational_state'])
    
    if states_with_data.empty:
        print("WARNING: No structures with conformational state data found!")
        print("Cannot create conformational state visualizations.")
        
        # Create alternative visualizations
        # 1. Overall EF1% distribution
        ax1 = axes[0, 0]
        performance_df['ef1'].hist(bins=20, ax=ax1)
        ax1.set_title('Overall EF1% Distribution')
        ax1.set_xlabel('EF1%')
        ax1.set_ylabel('Frequency')
        
        # 2. Overall AUC distribution
        ax2 = axes[0, 1]
        performance_df['auc'].hist(bins=20, ax=ax2)
        ax2.set_title('Overall AUC Distribution')
        ax2.set_xlabel('AUC')
        ax2.set_ylabel('Frequency')
        
        # 3. EF1% vs AUC scatter plot
        ax3 = axes[1, 0]
        ax3.scatter(performance_df['auc'], performance_df['ef1'], alpha=0.6)
        ax3.set_xlabel('AUC')
        ax3.set_ylabel('EF1%')
        ax3.set_title('EF1% vs AUC (All Structures)')
        
        # 4. Performance by kinase
        ax4 = axes[1, 1]
        kinase_means = performance_df.groupby('kinase')['ef1'].mean().sort_values(ascending=False)
        if len(kinase_means) > 0:
            kinase_means.head(10).plot(kind='bar', ax=ax4)
            ax4.set_title('Top 10 Kinases by Mean EF1%')
            ax4.set_xlabel('Kinase')
            ax4.set_ylabel('Mean EF1%')
            ax4.tick_params(axis='x', rotation=45)
    else:
        # Original visualizations with conformational states
        # 1. EF1% distribution by conformational state
        ax1 = axes[0, 0]
        sns.boxplot(data=states_with_data, x='conformational_state', y='ef1', ax=ax1)
        ax1.set_title('EF1% Distribution by Conformational State')
        ax1.set_xlabel('Conformational State')
        ax1.set_ylabel('EF1%')
        ax1.tick_params(axis='x', rotation=45)
        
        # 2. AUC distribution by conformational state
        ax2 = axes[0, 1]
        sns.boxplot(data=states_with_data, x='conformational_state', y='auc', ax=ax2)
        ax2.set_title('AUC Distribution by Conformational State')
        ax2.set_xlabel('Conformational State')
        ax2.set_ylabel('AUC')
        ax2.tick_params(axis='x', rotation=45)
        
        # 3. EF1% vs AUC scatter plot colored by conformational state
        ax3 = axes[1, 0]
        for state in states_with_data['conformational_state'].unique():
            state_data = states_with_data[states_with_data['conformational_state'] == state]
            ax3.scatter(state_data['auc'], state_data['ef1'], label=state, alpha=0.7)
        ax3.set_xlabel('AUC')
        ax3.set_ylabel('EF1%')
        ax3.set_title('EF1% vs AUC by Conformational State')
        ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # 4. Number of structures per conformational state
        ax4 = axes[1, 1]
        state_counts = states_with_data['conformational_state'].value_counts()
        if len(state_counts) > 0:
            state_counts.plot(kind='bar', ax=ax4)
            ax4.set_title('Number of Structures per Conformational State')
            ax4.set_xlabel('Conformational State')
            ax4.set_ylabel('Count')
            ax4.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig('structure_performance_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

def main():
    """Main analysis pipeline"""
    print("=== STRUCTURE-PERFORMANCE CORRELATION ANALYSIS ===")
    
    # File paths - adjust these to your actual file locations
    mapping_file = 'conformational_state_counts_structure_mapping.csv'
    docking_file = 'docking_results_clean_20251031_143707.csv'
    
    # Step 1: Load and merge data
    merged_df, mapping_df, docking_df = load_and_merge_data(mapping_file, docking_file)
    
    # Step 2: Calculate individual structure performance
    performance_df = calculate_structure_performance(merged_df)
    
    # Save performance results
    performance_df.to_csv('individual_structure_performance.csv', index=False)
    print("\nSaved individual structure performance to 'individual_structure_performance.csv'")
    
    # Step 3: Analyze conformational states
    state_stats, state_means = analyze_conformational_states(performance_df)
    
    # Step 4: Analyze within-state variation
    within_state_results = analyze_within_state_variation(performance_df)
    
    # Step 5: Identify top performers
    top_structures = identify_top_performers(performance_df)
    
    # Step 6: Create visualizations
    create_visualizations(performance_df)
    
    # Step 7: Summary insights
    print("\n=== KEY INSIGHTS ===")
    print(f"1. Total structures analyzed: {len(performance_df)}")
    print(f"2. Conformational states found: {performance_df['conformational_state'].nunique()}")
    print(f"3. Best performing state (by mean EF1%): {state_means['ef1'].idxmax()} (EF1%: {state_means['ef1'].max():.2f})")
    print(f"4. Worst performing state (by mean EF1%): {state_means['ef1'].idxmin()} (EF1%: {state_means['ef1'].min():.2f})")
    
    # Check for within-state variation
    high_variation_states = []
    for state, results in within_state_results.items():
        if results['ef1_range'] > 2.0:  # Arbitrary threshold for "high variation"
            high_variation_states.append((state, results['ef1_range']))
    
    if high_variation_states:
        print(f"5. States with high within-state variation (EF1% range > 2.0):")
        for state, variation in sorted(high_variation_states, key=lambda x: x[1], reverse=True):
            print(f"   {state}: range = {variation:.2f}")
    else:
        print("5. No states show high within-state variation in EF1%")
    
    return performance_df, state_stats, within_state_results

if __name__ == "__main__":
    performance_df, state_stats, within_state_results = main()
