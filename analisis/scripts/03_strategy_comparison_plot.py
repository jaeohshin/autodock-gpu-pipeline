#!/usr/bin/env python3
"""
Minimal Combined Figure Generator for Kinase Virtual Screening
=============================================================

Generates only the combined EF1% and AUC comparison figure.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from matplotlib import rcParams
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION AND STYLING
# ============================================================================

# Input data directory
DATA_DIR = "../data"

# Output directory for figures
FIGURE_DIR = "../data"
os.makedirs(FIGURE_DIR, exist_ok=True)

# Publication-quality styling parameters
STYLE_CONFIG = {
    'font.family': 'sans-serif',
    'font.size': 13,
    'axes.titlesize': 15,
    'axes.labelsize': 13,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
    'lines.linewidth': 1.5,
    'patch.linewidth': 1,
    'axes.linewidth': 1,
    'axes.grid': False,
    'axes.spines.top': False,
    'axes.spines.right': False,
}

# Apply styling
rcParams.update(STYLE_CONFIG)

# Color scheme for strategies
STRATEGY_COLORS = {
    'best_score': '#2E86AB',      # Professional blue
    'consensus': '#A23B72',       # Deep magenta
}

STRATEGY_LABELS = {
    'best_score': 'Best Score',
    'consensus': 'Consensus'
}

# ============================================================================
# DATA LOADING
# ============================================================================

def load_analysis_data():
    """Load the exported analysis data"""
    print("Loading analysis data...")
    
    # Load strategy comparison data (long format)
    comparison_file = os.path.join(DATA_DIR, "strategy_comparison_data.csv")
    if not os.path.exists(comparison_file):
        raise FileNotFoundError(f"Comparison data file not found: {comparison_file}")
    
    comparison_df = pd.read_csv(comparison_file)
    
    print(f"✓ Loaded comparison data: {len(comparison_df)} data points")
    
    return comparison_df

# ============================================================================
# FIGURE GENERATION
# ============================================================================

def create_combined_figure(comparison_df):
    """Create a combined figure with both EF1% and AUC comparisons"""
    print("Generating combined comparison figure...")
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # EF1% subplot (Panel A)
    ef_data = comparison_df[comparison_df['metric'] == 'EF1%'].copy()
    
    sns.boxplot(
        data=ef_data, x='strategy', y='value',
        palette=[STRATEGY_COLORS[s] for s in ef_data['strategy'].unique()],
        showfliers=False,
        ax=ax1, width=0.6, linewidth=1.5, notch=True
    )
    
    sns.stripplot(
        data=ef_data, x='strategy', y='value',
        color='black', alpha=0.5, size=6, ax=ax1, jitter=0.2
    )
    
    # AUC subplot (Panel B)
    auc_data = comparison_df[comparison_df['metric'] == 'AUC'].copy()
    
    sns.boxplot(
        data=auc_data, x='strategy', y='value',
        palette=[STRATEGY_COLORS[s] for s in auc_data['strategy'].unique()],
        ax=ax2, width=0.6, linewidth=1.5, notch=True
    )
    
    sns.stripplot(
        data=auc_data, x='strategy', y='value',
        color='black', alpha=0.5, size=6, ax=ax2, jitter=0.2
    )
    
    # Customize subplots
    strategies = list(auc_data['strategy'].unique())
    
    # EF1% subplot styling (Panel A)
    ax1.set_xlabel('Strategy', fontweight='bold')
    ax1.set_ylabel('EF1%', fontweight='bold')
    ax1.set_title('A) EF1% Performance', fontweight='bold', loc='left')
    ax1.set_xticklabels([STRATEGY_LABELS[s] for s in strategies])
    ax1.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.set_yticks([0, 5, 10, 15, 20])
    
    # AUC subplot styling (Panel B)
    ax2.set_xlabel('Strategy', fontweight='bold')
    ax2.set_ylabel('AUC', fontweight='bold')
    ax2.set_title('B) AUC Performance', fontweight='bold', loc='left')
    ax2.set_xticklabels([STRATEGY_LABELS[s] for s in strategies])
    ax2.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
    ax2.set_ylim(0.4, 0.9)
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save combined figure
    output_file = os.path.join(FIGURE_DIR, "Figure_Combined_Strategy_Comparison.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_file.replace('.png', '.pdf'), bbox_inches='tight')
    
    print(f"✓ Combined comparison figure saved: {output_file}")
    
    return fig

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main figure generation pipeline"""
    print("Combined Figure Generator")
    print("=" * 30)
    
    try:
        # Load data
        comparison_df = load_analysis_data()
        
        # Generate combined figure
        fig_combined = create_combined_figure(comparison_df)
        
        print("\n" + "=" * 30)
        print("Figure generation completed!")
        print(f"\nSaved: '{FIGURE_DIR}/Figure_Combined_Strategy_Comparison.pdf'")
        
        # Show figure
        plt.show()
        
    except Exception as e:
        print(f"Error generating figure: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
