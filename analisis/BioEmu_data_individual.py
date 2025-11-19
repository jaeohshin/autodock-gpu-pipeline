import pandas as pd
import numpy as np

# Read BioEmu data
bioemu_df = pd.read_csv('individual_structure_performance.csv')

print("="*70)
print("                 BIOEMU DATA PREPROCESSING")
print("="*70)

# Check the data
print(f"\nTotal number of structures: {len(bioemu_df)}")
print(f"\nColumns in BioEmu data: {bioemu_df.columns.tolist()}")
print(f"\nFirst few rows:")
print(bioemu_df.head())

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
bioemu_df['kinase_HGNC'] = bioemu_df['kinase'].replace(name_mapping)

# Count structures per kinase
kinase_counts = bioemu_df.groupby('kinase_HGNC').size().sort_values(ascending=False)

print("\n" + "="*70)
print("             STRUCTURES PER KINASE (HGNC Names)")
print("="*70)
print(kinase_counts)

# Check for kinases with less than 50 structures
less_than_50 = kinase_counts[kinase_counts < 50]
if len(less_than_50) > 0:
    print(f"\n⚠️  Kinases with fewer than 50 structures:")
    for kinase, count in less_than_50.items():
        print(f"   {kinase}: {count} structures")
else:
    print("\n✓ All kinases have exactly 50 structures")

print(f"\nTotal number of unique kinases: {len(kinase_counts)}")
print(f"\nKinases present: {sorted(kinase_counts.index.tolist())}")

# Check if all 26 target kinases are present
target_kinases = [
    'ABL1', 'AKT1', 'AKT2', 'BRAF', 'CDK2', 'CSF1R', 'EGFR', 'PTK2', 
    'FGFR1', 'IGF1R', 'JAK2', 'PRKCB', 'KIT', 'LCK', 'MAPKAPK2', 
    'MET', 'MAPK1', 'MAPK10', 'MAPK14', 'MAP2K1', 'PLK1', 'ROCK1', 
    'SRC', 'TGFBR1', 'KDR', 'WEE1'
]

present_kinases = set(kinase_counts.index.tolist())
target_set = set(target_kinases)

missing_kinases = target_set - present_kinases
extra_kinases = present_kinases - target_set

print("\n" + "="*70)
print("             KINASE COVERAGE CHECK")
print("="*70)

if missing_kinases:
    print(f"\n⚠️  Missing kinases (not in BioEmu data): {sorted(missing_kinases)}")
else:
    print("\n✓ All 26 target kinases are present!")

if extra_kinases:
    print(f"\nExtra kinases (in BioEmu but not in target list): {sorted(extra_kinases)}")

# Show different aggregation statistics for each kinase
print("\n" + "="*70)
print("      AGGREGATION OPTIONS FOR EF_1% (Preview)")
print("="*70)
print(f"{'Kinase':<12} {'Count':<8} {'Mean':<10} {'Median':<10} {'Range (Min-Max)':<20}")
print("-"*70)

aggregation_stats = bioemu_df.groupby('kinase_HGNC')['ef1'].agg(['count', 'max', 'mean', 'median', 'min'])
for kinase, row in aggregation_stats.iterrows():
    range_str = f"{row['min']:.3f} - {row['max']:.3f}"
    print(f"{kinase:<12} {int(row['count']):<8} {row['mean']:<10.3f} {row['median']:<10.3f} {range_str:<20}")

# Show conformational state distribution
print("\n" + "="*70)
print("           CONFORMATIONAL STATE DISTRIBUTION")
print("="*70)
conf_states = bioemu_df['conformational_state'].value_counts()
print(conf_states)

# Prepare aggregated datasets for different methods
print("\n" + "="*70)
print("        CREATING AGGREGATED DATASETS")
print("="*70)

# Maximum EF aggregation
bioemu_max = bioemu_df.groupby('kinase_HGNC').agg({
    'ef1': 'max',
    'ef5': 'max',
    'ef10': 'max'
}).reset_index()
bioemu_max.columns = ['Kinase', 'EF_1%', 'EF_5%', 'EF_10%']
# Round to 3 decimal places
bioemu_max[['EF_1%', 'EF_5%', 'EF_10%']] = bioemu_max[['EF_1%', 'EF_5%', 'EF_10%']].round(3)

# Mean EF aggregation
bioemu_mean = bioemu_df.groupby('kinase_HGNC').agg({
    'ef1': 'mean',
    'ef5': 'mean',
    'ef10': 'mean'
}).reset_index()
bioemu_mean.columns = ['Kinase', 'EF_1%', 'EF_5%', 'EF_10%']
# Round to 3 decimal places
bioemu_mean[['EF_1%', 'EF_5%', 'EF_10%']] = bioemu_mean[['EF_1%', 'EF_5%', 'EF_10%']].round(3)

# Median EF aggregation
bioemu_median = bioemu_df.groupby('kinase_HGNC').agg({
    'ef1': 'median',
    'ef5': 'median',
    'ef10': 'median'
}).reset_index()
bioemu_median.columns = ['Kinase', 'EF_1%', 'EF_5%', 'EF_10%']
# Round to 3 decimal places
bioemu_median[['EF_1%', 'EF_5%', 'EF_10%']] = bioemu_median[['EF_1%', 'EF_5%', 'EF_10%']].round(3)

# Save aggregated datasets
bioemu_max.to_csv('bioemu_max_ef.csv', index=False)
bioemu_mean.to_csv('bioemu_mean_ef.csv', index=False)
bioemu_median.to_csv('bioemu_median_ef.csv', index=False)

print("\n✓ Created 3 aggregated CSV files:")
print("  1. bioemu_max_ef.csv    (maximum EF across 50 structures)")
print("  2. bioemu_mean_ef.csv   (mean EF across 50 structures)")
print("  3. bioemu_median_ef.csv (median EF across 50 structures)")

# Show comparison preview
print("\n" + "="*70)
print("   DETAILED COMPARISON: Mean, Median, Range (EF_1%)")
print("="*70)
print(f"{'Kinase':<12} {'Mean':<10} {'Median':<10} {'Min':<10} {'Max':<10} {'Range':<10}")
print("-"*70)

comparison = pd.merge(bioemu_max, bioemu_mean, on='Kinase', suffixes=('_max', '_mean'))
comparison = pd.merge(comparison, bioemu_median, on='Kinase')

# Get min values for range display
bioemu_min = bioemu_df.groupby('kinase_HGNC').agg({'ef1': 'min'}).reset_index()
bioemu_min.columns = ['Kinase', 'EF_1%_min']
comparison = pd.merge(comparison, bioemu_min, on='Kinase')

for _, row in comparison.iterrows():
    range_val = row['EF_1%_max'] - row['EF_1%_min']
    print(f"{row['Kinase']:<12} {row['EF_1%_mean']:<10.3f} {row['EF_1%']:<10.3f} {row['EF_1%_min']:<10.3f} {row['EF_1%_max']:<10.3f} {range_val:<10.3f}")

print("\n" + "="*70)
print("NEXT STEPS:")
print("="*70)
print("1. Review the aggregation statistics above")
print("2. Decide which aggregation method to use (max, mean, or median)")
print("3. Use the corresponding CSV file for plotting comparison with your data")
print("="*70)
