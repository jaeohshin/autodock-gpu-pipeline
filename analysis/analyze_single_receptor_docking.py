import os
import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# === 설정 ===
kinase = "akt2"
receptor = "receptor_010"
dlg_dir = Path(f"/store/jaeohshin/work/dock/virtual_screening/docking_output/{kinase}/{receptor}")

# === 에너지 추출 함수 ===
def extract_best_energy_from_dlg(filepath):
    energies = []
    with open(filepath, 'r') as f:
        for line in f:
            if "Estimated Free Energy of Binding" in line:
                match = re.search(r"=\s*(-?\d+\.\d+)", line)
                if match:
                    energies.append(float(match.group(1)))
    return min(energies) if energies else None

# === 라벨링 함수 ===
def label_ligand(filename):
    if filename.startswith("actives_"):
        return "active"
    elif filename.startswith("decoys_"):
        return "decoy"
    return "unknown"

# === 데이터 수집 ===
records = []
for dlg_file in dlg_dir.glob("*.dlg"):
    ligand = dlg_file.stem
    energy = extract_best_energy_from_dlg(dlg_file)
    label = label_ligand(ligand)
    if energy is not None:
        records.append((ligand, energy, label))

df = pd.DataFrame(records, columns=["Ligand", "Energy", "Label"])
print(f"[INFO] {receptor} — total ligands: {len(df)} | actives: {(df['Label']=='active').sum()} | decoys: {(df['Label']=='decoy').sum()}")

# === 저장 (선택사항) ===
df.to_csv(f"{receptor}_all_ligands.csv", index=False)

# === 시각화 ===
plt.figure(figsize=(8, 4))
sns.histplot(data=df, x="Energy", hue="Label", bins=50, kde=True, stat="density", common_norm=False)
plt.title(f"{receptor} — Normalized Binding Energy Distribution")
plt.xlabel("Binding Energy (kcal/mol)")
plt.ylabel("Density")
plt.tight_layout()
plt.savefig(f"{receptor}_energy_distribution.png", dpi=150)
plt.show()

