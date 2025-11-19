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
    'legend.fontsize': 1618,    # Legend font size
    'figure.titlesize': 22    # Figure title size
})

# Create the dataset from the CSV content
data = {
    "Kinase": [
        "ABL1", "AKT1", "AKT2", "BRAF", "CDK2", "CSF1R", "EGFR", "FGFR1", "IGF1R", "JAK2",
        "KIT", "LCK", "MAPK2", "MK10", "MK14", "MET", "MP2K1", "PLK1", "KPCB", "ROCK1", "TGFR1", "WEE1"
    ],
    "EF@1% (You)": [
        9.9, 3.1, 19.7, 1.4, 10.8, 8.5, 7.4, 1.0, 12.9, 9.4,
        4.9, 9.8, 6.0, 3.9, 3.7, 8.4, 1.7, 1.9, 23.0, 2.1, 9.8, 54.4
    ],
    "EF@1% (Paper)": [
        7.6, 3.4, 22.8, 4.6, 5.9, 8.4, 5.3, 7.1, 5.4, 3.7,
        6.6, 8.8, 4.3, 4.7, 1.4, 7.3, 2.5, 4.6, 19.8, 3.0, 1.5, 47.7
    ],
    "EF@5% (You)": [
        3.1, 2.9, 8.8, 2.0, 4.8, 3.9, 3.4, 1.0, 5.2, 2.7,
        3.1, 4.9, 3.8, 3.5, 2.5, 6.4, 0.5, 3.4, 7.9, 2.6, 4.8, 15.7
    ],
    "EF@5% (Paper)": [
        3.2, 3.1, 8.2, 2.8, 3.2, 4.5, 2.4, 3.7, 4.1, 1.9,
        3.7, 4.4, 4.3, 4.6, 1.9, 6.7, 1.8, 1.5, 5.5, 3.6, 4.7, 13.7
    ],
    "EF@10% (You)": [
        2.6, 2.5, 5.2, 2.1, 3.4, 2.9, 2.6, 1.0, 3.3, 1.7,
        2.2, 3.4, 3.2, 3.0, 1.9, 4.9, 0.6, 3.0, 5.2, 2.6, 3.8, 8.0
    ],
    "EF@10% (Paper)": [
        2.8, 2.7, 5.3, 2.0, 2.5, 3.1, 2.0, 2.4, 3.0, 1.3,
        2.7, 3.0, 3.8, 3.3, 2.0, 4.3, 1.4, 1.5, 3.6, 3.3, 3.5, 7.2
    ]
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
    
    #x_lim = [x_min - x_padding, x_max + x_padding]
    #y_lim = [y_min - y_padding, y_max + y_padding]
    
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
    
#    ax.legend(fontsize=18, loc='upper right')
    
    # Add grid for better readability
    ax.grid(True, alpha=0.6, linewidth=1.5)
    
    # Make tick labels larger and bolder
    ax.tick_params(axis='both', which='major', labelsize=28, width=2, length=8)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save each figure as separate PNG
    filename = f'kinase_ef_{ef.replace("@", "_").replace("%", "percent")}_validation_simple.png'
    # plt.savefig(filename, 
    #             dpi=300, 
    #             bbox_inches='tight', 
    #             facecolor='white', 
    #             edgecolor='none',
    #             format='png')
    
    # Add PDF save right after
    filename_pdf = f'kinase_ef_{ef.replace("@", "_").replace("%", "percent")}_validation_simple.pdf'
    plt.savefig(filename_pdf, 
                bbox_inches='tight', 
                facecolor='white', 
                edgecolor='none',
                format='pdf')

    print(f"Plot saved as '{filename_pdf}'")
    
    # Show the plot
    #plt.show()

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