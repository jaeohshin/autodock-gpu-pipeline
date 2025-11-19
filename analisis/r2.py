#!/usr/bin/env python3
"""
Calculate R-squared analysis for specific kinases individually for user-defined EF percentages.
"""

import pandas as pd
import numpy as np
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
from scipy import stats
import os

# --- Utility Functions (Kept for Context) ---

def calculate_enrichment_factor(scores, labels, percentage=1):
    """Calculate enrichment factor at given percentage (Not used in R2 analysis, but kept)"""
    if len(scores) == 0: return 0.0
    df = pd.DataFrame({'score': scores, 'label': labels})
    df_sorted = df.sort_values('score')
    n_top = max(1, int(len(df_sorted) * percentage / 100))
    top_compounds = df_sorted.head(n_top)
    actives_in_top = (top_compounds['label'] == 'active').sum()
    total_actives = (df['label'] == 'active').sum()
    if total_actives == 0: return 0.0
    expected_actives = total_actives * percentage / 100
    if expected_actives == 0: return 0.0
    ef = actives_in_top / expected_actives
    return ef

# --- Core R-squared Calculation Function (Metric-Agnostic) ---

def calculate_per_kinase_r2_for_metric(performance_df, ef_metric, target_kinases=None):
    """
    Calculate R-squared for each kinase individually for a specific EF metric.
    """
    
    if ef_metric not in performance_df.columns:
        raise ValueError(f"Error: Required column '{ef_metric}' not found in the DataFrame. Check if your EF percentages match available columns.")
        
    df = performance_df.dropna(subset=['conformational_state']).copy()
    
    if target_kinases:
        df = df[df['kinase'].isin(target_kinases)]
        print(f"Analyzing kinases for {ef_metric}: {target_kinases}")
    else:
        print(f"Analyzing all kinases for {ef_metric}: {sorted(df['kinase'].unique())}")
    
    results = []
    
    print(f"\n{'='*80}")
    print(f"PER-KINASE R-SQUARED ANALYSIS for {ef_metric.upper()}")
    print(f"{'='*80}")
    
    for kinase in sorted(df['kinase'].unique()):
        kinase_data = df[df['kinase'] == kinase].copy()
        
        n_states = kinase_data['conformational_state'].nunique()
        n_structures = len(kinase_data)
        
        if n_states < 2:
            # print(f"\n{kinase}: SKIPPED - Only {n_states} conformational state(s)")
            continue
        
        if n_structures < 5:
            # print(f"\n{kinase}: SKIPPED - Only {n_structures} structures")
            continue
        
        # Calculate R-squared for this kinase
        state_means = kinase_data.groupby('conformational_state')[ef_metric].mean()
        predicted_col = 'predicted_' + ef_metric
        kinase_data[predicted_col] = kinase_data['conformational_state'].map(state_means)
        
        r2 = r2_score(kinase_data[ef_metric], kinase_data[predicted_col])
        
        # Variance decomposition
        total_var = np.var(kinase_data[ef_metric])
        between_var = np.var(kinase_data[predicted_col])
        # within_var = total_var - between_var # within_var is more robustly calculated via ANOVA residuals, but this approximation is fine for r2
        
        # ANOVA F-test
        groups = [group[ef_metric].values for name, group in kinase_data.groupby('conformational_state')]
        if len(groups) > 1:
            f_stat, p_value = stats.f_oneway(*groups)
        else:
            f_stat, p_value = 0, 1
        
        results.append({
            'kinase': kinase,
            'r2': r2,
            'n_structures': n_structures,
            'n_states': n_states,
            'f_statistic': f_stat,
            'p_value': p_value,
            'total_variance': total_var,
            'between_variance': between_var,
            'mean_' + ef_metric: kinase_data[ef_metric].mean(),
        })
        
        # --- Print detailed results for this kinase (Simplified for user input loop) ---
        print(f"\n{kinase.upper()} ({ef_metric.upper()}):")
        print(f"Structures: {n_structures}, States: {n_states}, R-squared: {r2:.4f}, p-value: {p_value:.3e}")
        
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        print(f"\n{'='*80}")
        print(f"SUMMARY TABLE for {ef_metric.upper()}")
        print(f"{'='*80}")
        
        results_df_sorted = results_df.sort_values('r2', ascending=False)
        
        print(f"{'Kinase':<8} {'R²':<8} {'Var Exp':<8} {'N Struct':<8} {'N States':<8} {'p-value':<10}")
        print(f"{'─'*70}")
        
        for _, row in results_df_sorted.iterrows():
            significance = "***" if row['p_value'] < 0.001 else "**" if row['p_value'] < 0.01 else "*" if row['p_value'] < 0.05 else ""
            print(f"{row['kinase']:<8} {row['r2']:<8.3f} {row['r2']*100:<8.1f}% {row['n_structures']:<8} {row['n_states']:<8} {row['p_value']:<10.2e}{significance}")
        
    return results_df

# --- Visualization Function (Metric-Agnostic) ---

def create_per_kinase_visualization(performance_df, results_df, ef_metric, target_kinases=None):
    """Create visualizations for per-kinase R-squared analysis"""
    
    df = performance_df.dropna(subset=['conformational_state']).copy()
    
    if target_kinases:
        df = df[df['kinase'].isin(target_kinases)]
        results_df = results_df[results_df['kinase'].isin(target_kinases)]
    
    n_kinases = len(results_df)
    if n_kinases == 0:
        print(f"No results to visualize for {ef_metric}.")
        return

    # Use a fixed grid for manageable plotting (max 6 plots)
    plot_kinases = results_df.sort_values('r2', ascending=False).head(6)
    n_plots = len(plot_kinases)

    cols = 3
    rows = int(np.ceil(n_plots / cols))
    
    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 5*rows))
    axes = axes.flatten()
    
    # Plot each kinase
    for i, (_, kinase_info) in enumerate(plot_kinases.iterrows()):
        kinase = kinase_info['kinase']
        r2 = kinase_info['r2']
        
        kinase_data = df[df['kinase'] == kinase].copy()
        
        state_means = kinase_data.groupby('conformational_state')[ef_metric].mean()
        predicted_col = 'predicted_' + ef_metric
        kinase_data[predicted_col] = kinase_data['conformational_state'].map(state_means)
        
        ax = axes[i]
        
        # Scatter plot with conformational state colors
        states = kinase_data['conformational_state'].unique()
        colors = plt.cm.Set3(np.linspace(0, 1, max(1, len(states))))
        
        for state, color in zip(states, colors):
            state_data = kinase_data[kinase_data['conformational_state'] == state]
            ax.scatter(state_data[predicted_col], state_data[ef_metric], 
                      label=state, alpha=0.7, s=50, color=color)
        
        # Perfect prediction line
        min_val = min(kinase_data[ef_metric].min(), kinase_data[predicted_col].min())
        max_val = max(kinase_data[ef_metric].max(), kinase_data[predicted_col].max())
        ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.8)
        
        ax.set_xlabel(f'Predicted {ef_metric.upper()}')
        ax.set_ylabel(f'Actual {ef_metric.upper()}')
        ax.set_title(f'{kinase} ({ef_metric.upper()})\nR² = {r2:.3f}')
        ax.grid(True, alpha=0.3)
        
    # Hide unused subplots
    for i in range(n_plots, len(axes)):
        axes[i].axis('off')
    
    filename = f'per_kinase_r2_analysis_top6_{ef_metric}.png'
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Top 6 kinase plots saved to '{filename}'")
    
# --- Main Function (UPDATED for User Input) ---

def main():
    """Main function for per-kinase R-squared analysis using user-defined EF percentages."""
    
    data_file = 'structure_performance_with_states.csv'
    try:
        performance_df = pd.read_csv(data_file)
        print(f"Loaded {len(performance_df)} structures from '{data_file}'")
    except FileNotFoundError:
        print(f"Error: '{data_file}' not found!")
        print("Please ensure the data file containing EF columns (e.g., 'ef1', 'ef5', 'ef10') exists.")
        return
    
    # 1. Get user input for EF percentages
    while True:
        user_input = input("\nEnter EF percentages to analyze (e.g., 1, 5, 10, or 0.5, 2): ")
        
        try:
            # Parse input string into a list of floats
            percentages = [float(p.strip()) for p in user_input.split(',')]
            
            # Convert to integers if they are whole numbers, otherwise keep as float
            metrics = [f'ef{int(p)}' if p == int(p) else f'ef{p}' for p in percentages]
            break
        except ValueError:
            print("Invalid input. Please enter numbers separated by commas.")
            continue
    
    print(f"Metrics to be analyzed: {', '.join(metrics)}")
    
    # Check if metrics exist in the DataFrame
    valid_metrics = [m for m in metrics if m in performance_df.columns]
    invalid_metrics = [m for m in metrics if m not in performance_df.columns]
    
    if not valid_metrics:
        print("\nError: None of the specified EF metrics were found in the data file.")
        print(f"Expected columns: {metrics}")
        return

    if invalid_metrics:
        print(f"\nWarning: Skipping non-existent columns: {', '.join(invalid_metrics)}")

    # Target kinases for detailed analysis (using a predefined list)
    target_kinases = ['ABL1', 'AKT1', 'AKT2', 'BRAF', 'CSF1R']
    
    all_results_combined = []
    
    for metric in valid_metrics:
        print(f"\n\n{'#'*80}")
        print(f"STARTING ANALYSIS FOR METRIC: {metric.upper()}")
        print(f"{'#'*80}")
        
        try:
            # 1. Analyze all kinases and save results
            all_results_df = calculate_per_kinase_r2_for_metric(performance_df, metric, target_kinases=None)
            
            if len(all_results_df) > 0:
                filename_all = f'all_kinases_r2_results_{metric}.csv'
                all_results_df.to_csv(filename_all, index=False)
                print(f"\nAll kinases results saved to '{filename_all}'")
                
                # Add a metric column for combination
                all_results_df['metric'] = metric
                all_results_combined.append(all_results_df)

                # 2. Create visualizations for a subset (using all_results_df head)
                # Filter down to the specific target kinases if they exist in the results
                target_results_df = all_results_df[all_results_df['kinase'].isin(target_kinases)]
                if len(target_results_df) > 0:
                    create_per_kinase_visualization(performance_df, target_results_df, metric, target_kinases)
                else:
                    print(f"Skipping detailed visualization: No target kinases found with enough data for {metric}.")

        except ValueError as e:
            print(e)
            print(f"Skipping {metric} analysis.")
        except Exception as e:
            print(f"An unexpected error occurred during {metric} analysis: {e}")
            
    if all_results_combined:
        final_df = pd.concat(all_results_combined)
        final_df.to_csv('combined_r2_results_summary.csv', index=False)
        print(f"\n\nFINAL SUMMARY: All results combined and saved to 'combined_r2_results_summary.csv'")

if __name__ == "__main__":
    main()
