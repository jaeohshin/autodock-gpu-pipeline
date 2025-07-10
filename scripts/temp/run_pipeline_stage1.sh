#!/bin/bash
set -e
source ~/miniconda3/etc/profile.d/conda.sh

KINASES=("BRAF" "EGFR" "FGFR1" "IGF1R" "JAK2" "KIT" "PRKCB" "LCK" "MAPKAPKxi2" "MET" "MAPKAPK1" "MAPKAPK10" "PLK1" "TGFBR1" "WEE1")  # Replace with real names

for kinase in "${KINASES[@]}"; do
    echo "🔵 Processing: $kinase"

    # Step 1 – BioEmu
    echo "→ Step 1: BioEmu"
    conda activate bioemu
    python scripts/run_bioemu.py "$kinase"

    # Step 2 – Extract backbone
    echo "→ Step 2: Extract backbone"
    python scripts/extract_backbone.py "output/$kinase/bioemu"

    # Step 3 – Glycine mutation
    echo "→ Step 3: Add Gly"
    python scripts/add_glycine_sidechains.py "output/$kinase/bioemu"

    # Step 4 – Energy Minimization
    echo "→ Step 4: Minimization"
    conda activate gro
    bash scripts/minimize_backbones.sh "output/$kinase/bioemu/gly/"

    # Step 5 – MolProbity scoring (optional, insert here if needed)

    # Step 6 – Restore real sequence
    echo "→ Step 6: Restore sequence"
    python scripts/correct_resnames.py "output/$kinase/bioemu/gly/"

    # Step 7 – FlowPacker
    echo "→ Step 7: FlowPacker"
    conda activate flowpacker
    FLOWPACKER_DIR="/store/jaeohshin/tools/flowpacker"
    BASE_CONFIG="${FLOWPACKER_DIR}/config/inference/_base.yaml"
    KINASE_CONFIG="${FLOWPACKER_DIR}/config/inference/${kinase}.yaml"
    FINAL_PATH="/store/jaeohshin/work/kine/output/$kinase/final_backbones"

    cp "$BASE_CONFIG" "$KINASE_CONFIG"
    sed -i "s|^  test_path:.*|  test_path: '$FINAL_PATH/'|" "$KINASE_CONFIG"

    pushd "$FLOWPACKER_DIR" > /dev/null
    python scripts/sampler_pdb.py "$kinase" "$kinase"
    popd > /dev/null
    rm "$KINASE_CONFIG"

    echo "✅ Done: $kinase"
done

