import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams

# Set publication-ready parameters
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
rcParams['font.size'] = 14
rcParams['axes.linewidth'] = 1.5
rcParams['xtick.major.width'] = 1.5
rcParams['ytick.major.width'] = 1.5
rcParams['xtick.minor.width'] = 1.0
rcParams['ytick.minor.width'] = 1.0
rcParams['lines.linewidth'] = 1.5
rcParams['patch.linewidth'] = 1.0
rcParams['figure.dpi'] = 300
rcParams['savefig.dpi'] = 300
rcParams['text.usetex'] = False  # Set to True if you have LaTeX installed

# Data
categories = [
    'DFGin-\nABAminus',
    'DFGin-\nBLAminus', 
    'DFGin-\nBLAplus',
    'DFGin-\nBLBminus',
    'DFGin-\nBLBplus',
    'DFGin-\nBLBtrans',
    'DFGin-\nUnassigned',
    'DFGinter-\nBABtrans',
    'DFGinter-\nUnassigned',
    'DFGout-\nBBAminus',
    'DFGout-\nUnassigned',
    'Unassigned-\nUnassigned'
]

# Data values (percentages)
all_human = [8, 55, 6, 5, 6, 8, 6, 0.5, 2, 8, 4, 2]
kinase_25 = [6, 30, 2, 3, 16, 10, 2, 0.5, 1, 12, 6, 3]
af2_standard = [4, 75, 2, 2, 8, 4, 3, 0.5, 1, 8, 4, 1]
alphafold3 = [5, 80, 1, 1, 10, 6, 2, 0.3, 1, 6, 4, 1]
bioemu = [10.6, 21.4, 4.8, 5.1, 5.4, 1.6, 27.9, 0.2, 6.1, 1.0, 8.0, 8.2]

# Create figure with publication dimensions (3.5" or 7" width for single/double column)
fig, ax = plt.subplots(figsize=(14, 8))  # 14x8 inches for double-column figure

# Set width of bars and positions
bar_width = 0.15
x = np.arange(len(categories))

# Publication-quality colors (colorblind-friendly palette)
colors = {
    'all_human': '#1f77b4',      # Blue
    'kinase_25': '#ff7f0e',      # Orange  
    'af2_standard': '#2ca02c',   # Green
    'alphafold3': '#d62728',     # Red
    'bioemu': '#9467bd'          # Purple
}

# Create bars with professional styling
bars1 = ax.bar(x - 2*bar_width, all_human, bar_width,
               label='All Human Kinase Structures', 
               color=colors['all_human'], 
               edgecolor='black', linewidth=0.8, alpha=0.9)

bars2 = ax.bar(x - bar_width, kinase_25, bar_width,
               label='25 Kinase Structures', 
               color=colors['kinase_25'], 
               edgecolor='black', linewidth=0.8, alpha=0.9)

bars3 = ax.bar(x, af2_standard, bar_width,
               label='Standard AF2 Predicted Structures', 
               color=colors['af2_standard'], 
               edgecolor='black', linewidth=0.8, alpha=0.9)

bars4 = ax.bar(x + bar_width, alphafold3, bar_width,
               label='AlphaFold3 Predicted Structures', 
               color=colors['alphafold3'], 
               edgecolor='black', linewidth=0.8, alpha=0.9)

bars5 = ax.bar(x + 2*bar_width, bioemu, bar_width,
               label='BioEmu-generated Structures', 
               color=colors['bioemu'], 
               edgecolor='black', linewidth=0.8, alpha=0.9)

# Customize axes with publication standards
ax.set_xlabel('Kinase Structural State', fontsize=14, fontweight='bold', labelpad=10)
ax.set_ylabel('Distribution (%)', fontsize=14, fontweight='bold', labelpad=10)
#ax.set_title('Distribution of Kinase Structural States Across Different Structure Sources', 
#             fontsize=16, fontweight='bold', pad=20)

# Set x-axis
ax.set_xticks(x)
ax.set_xticklabels(categories, rotation=45, ha='right', fontsize=12)

# Set y-axis with proper range and ticks
ax.set_ylim(0, 90)
ax.set_yticks(np.arange(0, 90, 10))
ax.tick_params(axis='y', labelsize=11)
ax.tick_params(axis='x', labelsize=11)

# Add professional grid
ax.grid(True, alpha=0.3, axis='y', linestyle='-', linewidth=0.8)
ax.set_axisbelow(True)

# Professional legend placement
legend = ax.legend(loc='upper right', 
                   frameon=True, 
                   fancybox=False, 
                   shadow=False, 
                   fontsize=14,
                   framealpha=0.95,
                   edgecolor='black',
                   facecolor='white')
legend.get_frame().set_linewidth(1.0)

# Remove top and right spines for cleaner look
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.5)
ax.spines['bottom'].set_linewidth(1.5)

# Adjust layout to prevent clipping
plt.tight_layout()

# Save with publication quality settings
plt.savefig('kincore_kinase_structural_states_publication.pdf', 
            dpi=300, 
            bbox_inches='tight',
            facecolor='white',
            edgecolor='none',
            transparent=False,
            format='pdf')

# plt.savefig('kinase_structural_states_publication.png', 
#             dpi=300, 
#             bbox_inches='tight',
#             facecolor='white',
#             edgecolor='none',
#             transparent=False,
#             format='png')

# # Also save as EPS for LaTeX (vector format)
# plt.savefig('kinase_structural_states_publication.eps', 
#             dpi=300, 
#             bbox_inches='tight',
#             facecolor='white',
#             edgecolor='none',
#             transparent=False,
#             format='eps')

# Display the plot
plt.show()

# Print data summary for verification
print("\nData Summary for Publication:")
print("=" * 60)
for i, cat in enumerate(categories):
    clean_cat = cat.replace('\n', '-')
    print(f"{clean_cat:<25}: "
          f"Human={all_human[i]:4.1f}% | "
          f"25K={kinase_25[i]:4.1f}% | " 
          f"AF2={af2_standard[i]:4.1f}% | "
          f"AF3={alphafold3[i]:4.1f}% | "
          f"BioEmu={bioemu[i]:4.1f}%")

print(f"\nTotal bars plotted: {len(categories) * 5}")
print("Files saved: PDF, PNG, and EPS formats")