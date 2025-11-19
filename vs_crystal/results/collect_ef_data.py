import pandas as pd
import os

# Define the base directory path
base_dir = "/data/work/dock/vs_crystal/results"

# Read kinase names from kinase.txt
kinase_file = os.path.join(base_dir, "kinase.txt")
kinases = []

with open(kinase_file, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 1:
            kinases.append(parts[0])

print(f"Found {len(kinases)} kinases")

# Initialize a dictionary to store EF data
ef_data = {
    'Kinase': [],
    'EF_1%': [],
    'EF_5%': [],
    'EF_10%': []
}

# Read enhancement factors from each kinase folder
for kinase in kinases:
    kinase_lower = kinase.lower()
    ef_file = os.path.join(base_dir, kinase_lower, "ef_summary.csv")
    
    if os.path.exists(ef_file):
        df = pd.read_csv(ef_file)
        
        ef_dict = {}
        for _, row in df.iterrows():
            top_pct = row['Top%']
            ef_value = row['EF']
            ef_dict[top_pct] = ef_value
        
        # Only add if all three percentages are present
        if 1 in ef_dict and 5 in ef_dict and 10 in ef_dict:
            ef_data['Kinase'].append(kinase)
            ef_data['EF_1%'].append(ef_dict[1])
            ef_data['EF_5%'].append(ef_dict[5])
            ef_data['EF_10%'].append(ef_dict[10])
            print(f"✓ {kinase}: EF@1%={ef_dict[1]:.2f}, EF@5%={ef_dict[5]:.2f}, EF@10%={ef_dict[10]:.2f}")
        else:
            print(f"✗ {kinase}: Missing some EF values")
    else:
        print(f"✗ {kinase}: ef_summary.csv not found")

# Create DataFrame
df_summary = pd.DataFrame(ef_data)

columns_to_format = ['EF_1%', 'EF_5%', 'EF_10%']
for col in columns_to_format:
    # Use .round() to change the underlying data in the DataFrame
    df_summary[col] = df_summary[col].round(3)
    
# Save summary to CSV
output_file = os.path.join(base_dir, "all_kinases_ef_summary.csv")
df_summary.to_csv(output_file, index=False)

print(f"\n{'='*60}")
print(f"Summary saved to: {output_file}")
print(f"Total kinases processed: {len(df_summary)}")
print(f"{'='*60}\n")
print(df_summary.to_string(index=False))
