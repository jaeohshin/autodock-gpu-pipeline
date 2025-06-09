#!/bin/bash

RECEPTOR_DIR="$1"
OUT_FILE="$2"

# Convert OUT_FILE to absolute path
OUT_FILE=$(realpath "$OUT_FILE")

# Ensure output directory exists
mkdir -p "$(dirname "$OUT_FILE")"

# Move to receptor directory
cd "$RECEPTOR_DIR" || exit 1

# Generate receptor list
ls receptor_*.pdb > receptor.list

# Write header
echo -e "receptor_i\treceptor_j\tTMscore" > "$OUT_FILE"

# Pairwise TMscore
for i in $(cat receptor.list); do
    for j in $(cat receptor.list); do
        result=$(TMscore "$i" "$j")
        score=$(echo "$result" | grep -E "^TM-score\s*=" | awk '{print $3}')
        [[ -z "$score" ]] && score="NA"
        echo -e "${i}\t${j}\t${score}" >> "$OUT_FILE"
    done
done

