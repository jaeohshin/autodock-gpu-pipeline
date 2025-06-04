import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# === Configuration ===
input_file = "results/ABL1_tmscore_matrix.tsv"  # Update path if needed
output_file = "results/ABL1_tmscore_matrix.png"  # Output image

# === Load TSV ===
df = pd.read_csv(input_file, sep="\t")

# === Pivot to matrix ===
matrix = df.pivot(index="receptor_i", columns="receptor_j", values="TMscore")
matrix = matrix.apply(pd.to_numeric, errors="coerce")  # Handle 'NA'

# === Plot ===
plt.figure(figsize=(12, 10))
sns.heatmap(matrix, cmap="viridis", square=True, linewidths=0.4,
            cbar_kws={"label": "TM-score"}, annot=False)

plt.title("TM-score Similarity Matrix (ABL1)", fontsize=14)
plt.xlabel("Receptor J")
plt.ylabel("Receptor I")
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout()

# === Save ===
plt.savefig(output_file, dpi=300)
print(f"[INFO] Saved heatmap to {output_file}")
plt.show()
