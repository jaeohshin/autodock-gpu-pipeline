#!/usr/bin/env python3
"""
Generate publication-quality figure for FGFR1 structure performance variability
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from sklearn.metrics import roc_curve, auc

# Configuration
CSV_FILE = "../data/output.csv"
TARGET_KINASE = "ABL1"

print(f"Loading data for {TARGET_KINASE}...")

# Load data
df = pd.read_csv(CSV_FILE)
kinase_data = df[df['kinase'] == TARGET_KINASE].copy()

print(f"Found {len(kinase_data)} docking results")
print(f"Structures: {kinase_data['structure_id'].nunique()}")
print(f"Compounds: {kinase_data['compound_id'].nunique()}")

# Calculate metrics for each structure
structure_metrics = []

for struct_id in sorted(kinase_data['structure_id'].unique()):
    struct_data = kinase_data[kinase_data['structure_id'] == struct_id]
    
    actives = struct_data[struct_data['compound_type'] == 'active']['docking_score'].values
    decoys = struct_data[struct_data['compound_type'] == 'decoy']['docking_score'].values
    
    if len(actives) > 0 and len(decoys) > 0:
        # Calculate AUC
        y_true = np.concatenate([np.ones(len(actives)), np.zeros(len(decoys))])
        y_scores = -np.concatenate([actives, decoys])
        
        fpr, tpr, _ = roc_curve(y_true, y_scores)
        auc_score = auc(fpr, tpr)
        
        # Calculate EF1%
        struct_scores = struct_data[['compound_id', 'compound_type', 'docking_score']].copy()
        struct_scores = struct_scores.sort_values('docking_score')
        
        n_actives = len(actives)
        n_total = len(struct_scores)
        top_1_percent = max(1, int(0.01 * n_total))
        actives_in_top_1 = (struct_scores.head(top_1_percent)['compound_type'] == 'active').sum()
        ef_1 = (actives_in_top_1 / top_1_percent) / (n_actives / n_total) if n_actives > 0 else 0
        
        structure_metrics.append({
            'structure_id': struct_id,
            'AUC': auc_score,
            'EF1%': ef_1,
            'n_actives': n_actives,
            'n_decoys': len(decoys)
        })

struct_df = pd.DataFrame(structure_metrics)

# Calculate consensus (ensemble) metrics
compound_scores = kinase_data.groupby(['compound_id', 'compound_type'])['docking_score'].mean().reset_index()
consensus_actives = compound_scores[compound_scores['compound_type'] == 'active']['docking_score'].values
consensus_decoys = compound_scores[compound_scores['compound_type'] == 'decoy']['docking_score'].values

consensus_y_true = np.concatenate([np.ones(len(consensus_actives)), np.zeros(len(consensus_decoys))])
consensus_y_scores = -np.concatenate([consensus_actives, consensus_decoys])
consensus_fpr, consensus_tpr, _ = roc_curve(consensus_y_true, consensus_y_scores)
consensus_auc = auc(consensus_fpr, consensus_tpr)

# Calculate consensus EF1%
consensus_scores = compound_scores.sort_values('docking_score')
n_total_compounds = len(consensus_scores)
top_1_percent_compounds = max(1, int(0.01 * n_total_compounds))
consensus_actives_in_top = (consensus_scores.head(top_1_percent_compounds)['compound_type'] == 'active').sum()
consensus_ef1 = (consensus_actives_in_top / top_1_percent_compounds) / (len(consensus_actives) / n_total_compounds)

print(f"\nConsensus metrics:")
print(f"  AUC: {consensus_auc:.3f}")
print(f"  EF1%: {consensus_ef1:.1f}")

# Statistics
mean_auc = struct_df['AUC'].mean()
std_auc = struct_df['AUC'].std()
mean_ef1 = struct_df['EF1%'].mean()
std_ef1 = struct_df['EF1%'].std()

print(f"\nStructure statistics:")
print(f"  AUC: {mean_auc:.3f} ± {std_auc:.3f}")
print(f"  EF1%: {mean_ef1:.1f} ± {std_ef1:.1f}")

# Create publication-quality figure
# Simpler, cleaner version
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Plot 1: EF1%
# x_positions 없이 바로
ax1.scatter(range(len(struct_df)), struct_df['EF1%'], s=90, color='#ff7f0e',  # orange
           alpha=0.8, edgecolors='black', linewidth=0.5, label='Individual')
ax1.axhline(y=consensus_ef1, color='#d62728', linestyle='--', 
           linewidth=2, label=f'Consensus: {consensus_ef1:.1f}')
# Mean line은 legend에 안 넣고 참고용으로만
ax1.axhline(y=mean_ef1, color='gray', linestyle=':', 
           linewidth=1.5, alpha=0.5)  # label 없음


ax1.set_ylabel('Enrichment Factor at 1%', fontsize=12, fontweight='bold')
ax1.set_ylim([-0.1, 7.5])
ax1.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
ax1.legend(loc='upper right', fontsize=10, frameon=True, edgecolor='black')

# Stats는 작게 corner에
ax1.text(0.02, 0.02, f'μ = {mean_ef1:.1f}, σ = {std_ef1:.1f}', 
        transform=ax1.transAxes, fontsize=9, 
        verticalalignment='bottom', style='italic', color='gray')

# Plot 2: AUC
ax2.scatter(range(len(struct_df)), struct_df['AUC'], s=90, color='#2ca02c',  # green
           alpha=0.8, edgecolors='black', linewidth=0.5, label='Individual')
ax2.axhline(y=consensus_auc, color='#d62728', linestyle='--', 
           linewidth=2, label=f'Consensus: {consensus_auc:.3f}')
ax2.axhline(y=mean_auc, color='gray', linestyle=':', 
           linewidth=1.5, alpha=0.5)  # label 없음

ax2.set_ylabel('Area Under ROC Curve', fontsize=12, fontweight='bold')
ax2.set_xlabel(f'{TARGET_KINASE.upper()} Structure Index', fontsize=12, fontweight='bold')
ax2.set_ylim([0.4, 0.8])
ax2.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
ax2.legend(loc='lower right', fontsize=10, frameon=True, edgecolor='black')

# Stats corner
ax2.text(0.02, 0.02, f'μ = {mean_auc:.3f}, σ = {std_auc:.3f}', 
        transform=ax2.transAxes, fontsize=9, 
        verticalalignment='bottom', style='italic', color='gray')

plt.tight_layout()

# Save in multiple formats
plt.savefig(f'../data/{TARGET_KINASE}_structure_performance.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'../data/{TARGET_KINASE}_structure_performance.png', dpi=300, bbox_inches='tight')

print(f"\nFigures saved:")
print(f"  - {TARGET_KINASE}_structure_performance.pdf")
print(f"  - {TARGET_KINASE}_structure_performance.png")

# Save metrics table
struct_df.to_csv(f'../data/{TARGET_KINASE}_structure_metrics.csv', index=False)
print(f"  - {TARGET_KINASE}_structure_metrics.csv")

plt.show()
