#!/bin/bash
for kinase in braf fak1 fgfr1 igf1r jak2 kit mapk2 met mk01 wee1; do
    echo "Running analysis for $kinase..."
    python anal_ens_dock.py --kinase $kinase --max_receptor 50
done
