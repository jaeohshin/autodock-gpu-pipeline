import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set matplotlib parameters for better text visibility
plt.rcParams.update({
    'font.size': 16,          # Base font size
    'axes.titlesize': 20,     # Title font size
    'axes.labelsize': 20,     # Axis label font size
    'xtick.labelsize': 20,    # X tick label size
    'ytick.labelsize': 20,    # Y tick label size
    'legend.fontsize': 16,    # Legend font size (fixed typo from 1618)
    'figure.titlesize': 22    # Figure title size
})

# Read CSV files
your_data = pd.read_csv('all_kinases_ef_summary.csv')
song_data = pd.read_csv('ef_crystal_song.csv')

# Define mapping from DUD-E names to HGNC symbols for plotting
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

# Apply name mapping to both datasets for consistency
your_data['Kinase_HGNC'] = your_data['Kinase'].replace(name_mapping)
song_data['Kinase_HGNC'] = song_data['Kinase'].replace(name_mapping)

# Merge datasets on kinase names
merged_df = pd.merge(
    your_data, 
    song_data, 
    on='Kinase_HGNC', 
    suffixes=('_You', '_Song'),
    how='inner'  # Only keep kinases present in both datasets
)

print(f"Found {len(merged_df)} common kinases between the two datasets")
print(f"Common kinases: {sorted(merged_df['Kinase_HGNC'].tolist())}\n")

# Prepare data for plotting
data = {
    "Kinase": merged_df['Kinase_HGNC'].tolist(),
    "EF@1% (You)": merged_df['EF_1%_You'].tolist(),
    "EF@1% (Paper)": merged_df['EF_1%_Song'].tolist(),
    "EF@5% (You)": merged_df['EF_5%_You'].tolist(),
    "EF@5% (Paper)": merged_df['EF_5%_Song'].tolist(),
    "EF@10% (You)": merged_df['EF_10%_You'].tolist(),
    "EF@10% (Paper)": merged_df['EF_10%_Song'].tolist()
}

df = pd.DataFrame(data)

# Set plot style
sns.set(style="whitegrid")

# Create separate figures for each threshold
thresholds = ["EF@1%", "EF@5%", "EF@10%"]
colors = ["steelblue", "forestgreen", "darkorange"]

for i, ef in enumerate(thresholds):
    # Create new figure for each threshold
    fig, ax = plt.subplots(1, 1, figsize=(13, 11))
    
    # Get data for current threshold
    x_data = df[f"{ef} (You)"]
    y_data = df[f"{ef} (Paper)"]
    
    # Calculate axis limits with some padding
    x_min, x_max = x_data.min(), x_data.max()
    y_min, y_max = y_data.min(), y_data.max()
    
    # Add padding to accommodate labels
    x_range = x_max - x_min
    y_range = y_max - y_min
    x_padding = max(0.1 * x_range, 0.2)  # Padding for labels
    y_padding = max(0.1 * y_range, 0.2)  # Padding for labels
    
    x_lim = [0, x_max + x_padding]
    y_lim = [0, y_max + y_padding]
    
    # Make axes limits symmetric around the diagonal for better comparison
    overall_min = min(x_lim[0], y_lim[0])
    overall_max = max(x_lim[1], y_lim[1])
    
    # Create scatter plot with larger markers
    sns.scatterplot(
        x=f"{ef} (You)", y=f"{ef} (Paper)", data=df, ax=ax,
        s=200, color=colors[i], edgecolor="black", linewidth=2, alpha=0.8
    )
    
    # Add diagonal line (y = x) over the appropriate range
    ax.plot([overall_min, overall_max], [overall_min, overall_max],
            'r--', label="y = x", linewidth=3, alpha=0.4)
    
    # Set axis limits
    ax.set_xlim(overall_min, overall_max)
    ax.set_ylim(overall_min, overall_max)
    
    # Add simple text labels for all kinases
    for idx, row in df.iterrows():
        x_val = row[f"{ef} (You)"]
        y_val = row[f"{ef} (Paper)"]
        kinase_name = row['Kinase']
        
        # Add simple text annotation slightly offset from the point
        ax.text(x_val + 0.1, y_val + 0.1, kinase_name, 
                fontsize=12,
                fontweight='bold',
                alpha=0.6,
                color='black',
                ha='left', va='bottom')
    
    # Calculate and display correlation coefficient
    correlation = np.corrcoef(x_data, y_data)[0, 1]
    
    # Set labels and title with larger fonts
    # Extract the percentage dynamically (e.g., "1%" from "EF@1%")
    threshold_pct = ef.split('@')[1]
    threshold_label = threshold_pct.replace('%', r'\%')

    ax.set_xlabel(f"EF$_{{{threshold_label}}}$ (This Study)", fontsize=26, fontweight='bold')
    ax.set_ylabel(f"EF$_{{{threshold_label}}}$ (Song et al.)", fontsize=26, fontweight='bold')
    ax.text(0.1, 0.9, f'r = {correlation:.3f}', 
            transform=ax.transAxes,
            fontsize=26,
            fontweight='bold',
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='black'))
    
    # Add grid for better readability
    ax.grid(True, alpha=0.6, linewidth=1.5)
    
    # Make tick labels larger and bolder
    ax.tick_params(axis='both', which='major', labelsize=28, width=2, length=8)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save as PDF
    filename_pdf = f'kinase_ef_{ef.replace("@", "_").replace("%", "percent")}_validation_simple.pdf'
    plt.savefig(filename_pdf, 
                bbox_inches='tight', 
                facecolor='white', 
                edgecolor='none',
                format='pdf')

    print(f"Plot saved as '{filename_pdf}'")
    
    # Optionally show the plot
    # plt.show()
    plt.close()

# Print summary statistics for each threshold
print("\n" + "="*60)
print("           KINASE ENHANCEMENT FACTOR VALIDATION")
print("                  SUMMARY STATISTICS")
print("="*60)

for ef in ["EF@1%", "EF@5%", "EF@10%"]:
    x_data = df[f"{ef} (You)"]
    y_data = df[f"{ef} (Paper)"]
    correlation = np.corrcoef(x_data, y_data)[0, 1]
    mean_abs_diff = np.mean(np.abs(x_data - y_data))
    
    print(f"\n{ef}:")
    print(f"  Correlation coefficient: {correlation:.3f}")
    print(f"  Mean absolute difference: {mean_abs_diff:.2f}")
    print(f"  Your EF range: {x_data.min():.1f} - {x_data.max():.1f}")
    print(f"  Paper EF range: {y_data.min():.1f} - {y_data.max():.1f}")
    
    # Find the kinase with the largest absolute difference
    max_diff_idx = np.argmax(np.abs(x_data - y_data))
    max_diff_kinase = df.iloc[max_diff_idx]['Kinase']
    max_diff_value = np.abs(x_data - y_data).iloc[max_diff_idx]
    print(f"  Largest discrepancy: {max_diff_kinase} ({max_diff_value:.1f})")

print("="*60)