import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set matplotlib parameters
plt.rcParams.update({
    'font.size': 16,
    'axes.titlesize': 20,
    'axes.labelsize': 20,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'legend.fontsize': 16,
    'figure.titlesize': 22
})

# Read BioEmu data
bioemu_df = pd.read_csv('individual_structure_performance.csv')

# Read your single crystal structure data
your_data = pd.read_csv('all_kinases_ef_summary.csv')

# Define mapping from DUD-E names to HGNC symbols
name_mapping = {
    'KPCB': 'PRKCB',
    'MK01': 'MAPK1',
    'MK10': 'MAPK10',
    'MK14': 'MAPK14',
    'MP2K1': 'MAP2K1',
    'MAPK2': 'MAPKAPK2',
    'VGFR2': 'KDR',
    'FAK1': 'PTK2',
    'TGFR1': 'TGFBR1'
}

# Apply name mapping
bioemu_df['Kinase'] = bioemu_df['kinase'].replace(name_mapping)
your_data['Kinase'] = your_data['Kinase'].replace(name_mapping)

# Create plots for each EF threshold
ef_columns = [
    ('ef1', 'EF_1%', '1\%'),
    ('ef5', 'EF_5%', '5\%'),
    ('ef10', 'EF_10%', '10\%')
]
colors = ["steelblue", "forestgreen", "darkorange"]

for i, (bioemu_col, crystal_col, label_pct) in enumerate(ef_columns):
    fig, ax = plt.subplots(1, 1, figsize=(13, 11))
    
    # Merge BioEmu data with crystal structure data
    # Each BioEmu structure will be matched with its kinase's crystal EF
    merged_data = pd.merge(
        bioemu_df[['Kinase', bioemu_col]],
        your_data[['Kinase', crystal_col]],
        on='Kinase',
        how='inner'
    )
    
    print(f"\n{crystal_col}: {len(merged_data)} BioEmu structures matched")
    
    # Get X and Y data
    x_data = merged_data[crystal_col]  # Crystal structure (one per kinase)
    y_data = merged_data[bioemu_col]   # BioEmu ensemble (50 per kinase)
    
    # Calculate axis limits
    x_min, x_max = x_data.min(), x_data.max()
    y_min, y_max = y_data.min(), y_data.max()
    
    x_range = x_max - x_min
    y_range = y_max - y_min
    x_padding = max(0.1 * x_range, 0.5)
    y_padding = max(0.1 * y_range, 0.5)
    
    x_lim = [0, x_max + x_padding]
    y_lim = [0, y_max + y_padding]
    
    # Make axes symmetric
    overall_min = min(x_lim[0], y_lim[0])
    overall_max = max(x_lim[1], y_lim[1])
    
    # Create scatter plot
    ax.scatter(
        x_data, 
        y_data,
        s=100,
        color=colors[i],
        alpha=0.5,
        edgecolor='black',
        linewidth=0.5
    )
    
    # Add diagonal line (y = x)
    ax.plot(
        [overall_min, overall_max], 
        [overall_min, overall_max],
        'r--', 
        linewidth=3, 
        alpha=0.4,
        label='y = x'
    )
    
    # Set axis limits
    ax.set_xlim(overall_min, overall_max)
    ax.set_ylim(overall_min, overall_max)
    
    # Calculate correlation
    correlation = np.corrcoef(x_data, y_data)[0, 1]
    
    # Calculate how many BioEmu structures are better than crystal
    better_count = (y_data > x_data).sum()
    total_count = len(y_data)
    better_pct = (better_count / total_count) * 100
    
    # Set labels
    ax.set_xlabel(f"EF$_{{{label_pct}}}$ (Crystal Structure)", fontsize=26, fontweight='bold')
    ax.set_ylabel(f"EF$_{{{label_pct}}}$ (BioEmu Ensemble)", fontsize=26, fontweight='bold')
    
    # Add correlation and improvement stats
    stats_text = f'r = {correlation:.3f}\n'
    stats_text += f'{better_count}/{total_count} ({better_pct:.1f}%)\nabove diagonal'
    
    ax.text(0.05, 0.95, stats_text,
            transform=ax.transAxes,
            fontsize=22,
            fontweight='bold',
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='black'))
    
    # Add grid
    ax.grid(True, alpha=0.6, linewidth=1.5)
    
    # Make tick labels larger and bolder
    ax.tick_params(axis='both', which='major', labelsize=28, width=2, length=8)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save as PDF
    filename_pdf = f'crystal_vs_bioemu_{crystal_col.lower().replace("%", "percent")}.pdf'
    plt.savefig(filename_pdf,
                bbox_inches='tight',
                facecolor='white',
                edgecolor='none',
                format='pdf')
    
    print(f"Plot saved as '{filename_pdf}'")
    plt.close()

# Print summary statistics
print("\n" + "="*70)
print("           CRYSTAL vs BIOEMU ENSEMBLE COMPARISON")
print("="*70)

for bioemu_col, crystal_col, _ in ef_columns:
    merged_data = pd.merge(
        bioemu_df[['Kinase', bioemu_col]],
        your_data[['Kinase', crystal_col]],
        on='Kinase',
        how='inner'
    )
    
    x_data = merged_data[crystal_col]
    y_data = merged_data[bioemu_col]
    
    correlation = np.corrcoef(x_data, y_data)[0, 1]
    better_count = (y_data > x_data).sum()
    total_count = len(y_data)
    better_pct = (better_count / total_count) * 100
    
    mean_improvement = (y_data - x_data).mean()
    
    print(f"\n{crystal_col}:")
    print(f"  Correlation: {correlation:.3f}")
    print(f"  BioEmu better than Crystal: {better_count}/{total_count} ({better_pct:.1f}%)")
    print(f"  Mean improvement: {mean_improvement:+.3f}")
    print(f"  Crystal range: {x_data.min():.3f} - {x_data.max():.3f}")
    print(f"  BioEmu range: {y_data.min():.3f} - {y_data.max():.3f}")

print("="*70)