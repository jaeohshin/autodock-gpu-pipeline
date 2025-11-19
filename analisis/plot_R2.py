#!/usr/bin/env python3
"""
Simplified R-squared analysis per kinase
"""

import pandas as pd
import numpy as np
from sklearn.metrics import r2_score
from scipy import stats
import matplotlib.pyplot as plt

def calculate_r2_by_kinase(performance_df):
    """Calculate R² for each kinase"""
    
    # Remove structures without conformational state data
    df = performance_df.dropna(subset=['conformational_state']).copy()
    
    results = []
    
    for kinase in sorted(df['kinase'].unique()):
        kinase_data = df[df['kinase'] == kinase].copy()
        
        # Need at least 2 states and 5 structures
        n_states = kinase_data['conformational_state'].nunique()
        n_structures = len(kinase_data)
        
        if n_states < 2 or n_structures < 5:
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
    
    return pd.DataFrame(results)

def plot_r2_by_kinase(results_df):
    """Create bar plot of R² by kinase"""
    
    # Sort by R²
    results_sorted = results_df.sort_values('r2', ascending=True)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create bars with color coding
    bars = ax.barh(results_sorted['kinase'], results_sorted['r2'])
    
    for i, bar in enumerate(bars):
        r2_val = results_sorted.iloc[i]['r2']
        if r2_val > 0.3:
            bar.set_color('darkgreen')
        elif r2_val > 0.1:
            bar.set_color('orange')
        else:
            bar.set_color('red')
    
    ax.set_xlabel('R-squared', fontsize=13)
    ax.set_ylabel('Kinase', fontsize=14)
    ax.set_title('R² by Kinase: Conformational State Predictive Power\n(Green: High >0.3, Orange: Medium 0.1-0.3, Red: Low <0.1)', 
                 fontsize=14, pad=20)
    ax.grid(True, alpha=0.3, axis='x')
    ax.tick_params(axis='y', labelsize=13)

    plt.tight_layout()
    plt.savefig('r2_by_kinase.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('r2_by_kinase.png', dpi=300, bbox_inches='tight')
    plt.show()

def print_summary(results_df):
    """Print summary statistics"""
    
    print(f"\n{'='*80}")
    print(f"R² BY KINASE SUMMARY")
    print(f"{'='*80}\n")
    
    # Sort by R²
    results_sorted = results_df.sort_values('r2', ascending=False)
    
    print(f"{'Kinase':<10} {'R²':<8} {'Var%':<8} {'N':<6} {'States':<8} {'p-value':<12} {'Sig':<5}")
    print(f"{'─'*70}")
    
    for _, row in results_sorted.iterrows():
        sig = "***" if row['p_value'] < 0.001 else "**" if row['p_value'] < 0.01 else "*" if row['p_value'] < 0.05 else "ns"
        print(f"{row['kinase']:<10} {row['r2']:<8.3f} {row['var_explained_pct']:<8.1f} {row['n_structures']:<6} {row['n_states']:<8} {row['p_value']:<12.2e} {sig:<5}")
    
    # Category summary
    high = results_df[results_df['r2'] > 0.3]
    medium = results_df[(results_df['r2'] > 0.1) & (results_df['r2'] <= 0.3)]
    low = results_df[results_df['r2'] <= 0.1]
    
    print(f"\n{'='*80}")
    print(f"CATEGORIES")
    print(f"{'='*80}")
    print(f"HIGH (R² > 0.3):   {len(high):2d} kinases ({len(high)/len(results_df)*100:5.1f}%)")
    print(f"MEDIUM (0.1-0.3):  {len(medium):2d} kinases ({len(medium)/len(results_df)*100:5.1f}%)")
    print(f"LOW (R² < 0.1):    {len(low):2d} kinases ({len(low)/len(results_df)*100:5.1f}%)")
    
    print(f"\nMean R²: {results_df['r2'].mean():.3f}")
    print(f"Median R²: {results_df['r2'].median():.3f}")

def main():
    """Main analysis"""
    
    # Load data
    try:
        df = pd.read_csv('structure_performance_with_states.csv')
        print(f"Loaded {len(df)} structures")
    except FileNotFoundError:
        print("Error: 'structure_performance_with_states.csv' not found!")
        return
    
    # Calculate R² by kinase
    results = calculate_r2_by_kinase(df)
    
    if len(results) > 0:
        # Print summary
        print_summary(results)
        
        # Create plot
        plot_r2_by_kinase(results)
        
        # Save results
        results.to_csv('r2_by_kinase.csv', index=False)
        print(f"\nResults saved to 'r2_by_kinase.csv'")
    else:
        print("No valid kinases found for analysis!")

if __name__ == "__main__":
    main()