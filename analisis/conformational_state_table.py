#!/usr/bin/env python3
"""
Generate a table showing EF1% values for different conformational states
Similar to the table in the research paper
"""

import pandas as pd
import numpy as np

def create_conformational_state_table(performance_df, kinase_filter=None, output_file='conformational_state_performance_table.csv'):
    """
    Create a table showing EF1% performance for different conformational states
    
    Args:
        performance_df: DataFrame with individual structure performance
        kinase_filter: If specified, filter to only this kinase (e.g., 'FGFR1')
        output_file: Output CSV filename
    """
    
    # Filter to specific kinase if requested
    if kinase_filter:
        df = performance_df[performance_df['kinase'] == kinase_filter].copy()
        print(f"Filtering to kinase: {kinase_filter}")
        print(f"Structures for {kinase_filter}: {len(df)}")
    else:
        df = performance_df.copy()
        print("Using all kinases")
    
    # Remove rows with missing conformational states
    df = df.dropna(subset=['conformational_state'])
    
    if len(df) == 0:
        print("No data with conformational states found!")
        return None
    
    print(f"Total structures with conformational states: {len(df)}")
    
    # Calculate statistics for each conformational state
    state_stats = df.groupby('conformational_state').agg({
        'ef1': ['count', 'mean', 'std', 'min', 'max'],
        'auc': ['mean', 'std']
    }).round(3)
    
    # Flatten column names
    state_stats.columns = ['_'.join(col).strip() for col in state_stats.columns]
    
    # Create a clean summary table
    summary_table = []
    
    for state in state_stats.index:
        count = state_stats.loc[state, 'ef1_count']
        mean_ef1 = state_stats.loc[state, 'ef1_mean']
        std_ef1 = state_stats.loc[state, 'ef1_std']
        min_ef1 = state_stats.loc[state, 'ef1_min']
        max_ef1 = state_stats.loc[state, 'ef1_max']
        mean_auc = state_stats.loc[state, 'auc_mean']
        
        summary_table.append({
            'Structural State': state,
            'Count': int(count),
            'EF1% Mean': mean_ef1,
            'EF1% Std': std_ef1,
            'EF1% Min': min_ef1,
            'EF1% Max': max_ef1,
            'AUC Mean': mean_auc
        })
    
    # Convert to DataFrame and sort by EF1% mean (descending)
    table_df = pd.DataFrame(summary_table)
    table_df = table_df.sort_values('EF1% Mean', ascending=False)
    
    # Print formatted table
    print(f"\n{'='*60}")
    if kinase_filter:
        print(f"CONFORMATIONAL STATE PERFORMANCE TABLE FOR {kinase_filter}")
    else:
        print("CONFORMATIONAL STATE PERFORMANCE TABLE (ALL KINASES)")
    print(f"{'='*60}")
    
    print(f"{'Structural State':<20} {'Count':<6} {'EF1%':<8} {'Std':<6} {'Range':<10} {'AUC':<6}")
    print(f"{'-'*60}")
    
    for _, row in table_df.iterrows():
        state = row['Structural State']
        count = row['Count']
        ef1_mean = row['EF1% Mean']
        ef1_std = row['EF1% Std']
        ef1_range = f"{row['EF1% Min']:.1f}-{row['EF1% Max']:.1f}"
        auc_mean = row['AUC Mean']
        
        print(f"{state:<20} {count:<6} {ef1_mean:<8.1f} {ef1_std:<6.1f} {ef1_range:<10} {auc_mean:<6.3f}")
    
    # Save to CSV
    table_df.to_csv(output_file, index=False)
    print(f"\nTable saved to: {output_file}")
    
    # Additional analysis: Show which states have high variation
    print(f"\n{'='*40}")
    print("WITHIN-STATE VARIATION ANALYSIS")
    print(f"{'='*40}")
    
    high_variation_states = table_df[table_df['EF1% Std'] > 2.0]  # Arbitrary threshold
    if len(high_variation_states) > 0:
        print("States with high variation (Std > 2.0):")
        for _, row in high_variation_states.iterrows():
            state = row['Structural State']
            std_val = row['EF1% Std']
            count = row['Count']
            print(f"  {state}: {std_val:.1f} std ({count} structures)")
    else:
        print("No states show high within-state variation (Std > 2.0)")
    
    # Show best and worst individual performers for top states
    print(f"\n{'='*40}")
    print("BEST/WORST PERFORMERS IN TOP STATES")
    print(f"{'='*40}")
    
    top_states = table_df.head(3)['Structural State'].values
    
    for state in top_states:
        state_structures = df[df['conformational_state'] == state]
        if len(state_structures) > 1:
            best = state_structures.loc[state_structures['ef1'].idxmax()]
            worst = state_structures.loc[state_structures['ef1'].idxmin()]
            
            print(f"\n{state} ({len(state_structures)} structures):")
            print(f"  Best:  {best['kinase']}-{best['structure_id']:02d} (EF1%: {best['ef1']:.1f})")
            print(f"  Worst: {worst['kinase']}-{worst['structure_id']:02d} (EF1%: {worst['ef1']:.1f})")
    
    return table_df

def compare_with_literature(table_df, literature_values=None):
    """
    Compare your results with literature values (like the table in the figure)
    
    Args:
        table_df: Your conformational state performance table
        literature_values: Dict of {state: ef1_value} from literature
    """
    
    if literature_values is None:
        # Values from the figure you showed (for FGFR1)
        literature_values = {
            'DFGin-BLAminus': 2.9,
            'DFGin-BLAplus': 6.4,
            'DFGin-BLBminus': 0.7,
            'DFGin-BLBplus': 8.6,
            'DFGin-BLBtrans': 2.1,
            'DFGout-BBAminus': 3.6,
            'DFGout-Unassigned': 1.4
        }
    
    print(f"\n{'='*50}")
    print("COMPARISON WITH LITERATURE VALUES")
    print(f"{'='*50}")
    print(f"{'State':<20} {'Literature':<12} {'Your Data':<12} {'Difference':<12}")
    print(f"{'-'*56}")
    
    for state, lit_value in literature_values.items():
        your_row = table_df[table_df['Structural State'] == state]
        if len(your_row) > 0:
            your_value = your_row['EF1% Mean'].iloc[0]
            difference = your_value - lit_value
            print(f"{state:<20} {lit_value:<12.1f} {your_value:<12.1f} {difference:<12.1f}")
        else:
            print(f"{state:<20} {lit_value:<12.1f} {'N/A':<12} {'N/A':<12}")

def main():
    """Main function to generate conformational state performance table"""
    
    # Load the performance data (from your previous analysis)
    try:
        performance_df = pd.read_csv('individual_structure_performance.csv')
        print(f"Loaded performance data: {len(performance_df)} structures")
    except FileNotFoundError:
        print("Error: 'individual_structure_performance.csv' not found!")
        print("Please run the structure analysis script first.")
        return
    
    # Option 1: Create table for all kinases
    print("\n" + "="*60)
    print("OPTION 1: ALL KINASES COMBINED")
    print("="*60)
    
    all_kinases_table = create_conformational_state_table(
        performance_df, 
        kinase_filter=None,
        output_file='conformational_state_table_all_kinases.csv'
    )
    
    # Option 2: Create table for specific kinase (like FGFR1 in the paper)
    available_kinases = performance_df['kinase'].unique()
    print(f"\nAvailable kinases: {sorted(available_kinases)}")
    
    # Try FGFR1 first (to match the paper), then any kinase with good data
    target_kinases = ['FGFR1', 'ABL1', 'AKT1', 'AKT2']  # Priority order
    
    for kinase in target_kinases:
        if kinase in available_kinases:
            kinase_data = performance_df[performance_df['kinase'] == kinase]
            kinase_with_states = kinase_data.dropna(subset=['conformational_state'])
            
            if len(kinase_with_states) > 10:  # Need reasonable amount of data
                print(f"\n" + "="*60)
                print(f"OPTION 2: {kinase} ONLY (like the paper)")
                print("="*60)
                
                kinase_table = create_conformational_state_table(
                    performance_df,
                    kinase_filter=kinase,
                    output_file=f'conformational_state_table_{kinase.lower()}.csv'
                )
                
                # Compare with literature if it's FGFR1
                if kinase == 'FGFR1' and kinase_table is not None:
                    compare_with_literature(kinase_table)
                
                break
    
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print("Created tables showing EF1% performance by conformational state")
    print("This replicates the analysis style from the research paper")
    print("Key findings will show which conformational states are optimal for VS")

if __name__ == "__main__":
    main()
