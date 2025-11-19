#!/bin/bash
#To validate BioEmu generated structures, calculate TMscore
#and collect top 20 by TM-score.

reference="ABL1.renum.pdb"
output_file="tmscore.txt"
offset=224  # because BioEmu starts at 1, and crystal starts at 225
top_dir="top20"
top_list="top20.txt"

mkdir -p "$top_dir"
echo -e "Structure\tTM-score\tAlignedResidues\tRMSD" > "$output_file"

for target in receptor_[0-9][0-9][0-9][0-9].pdb; do
    echo "[INFO] Processing $target..."

    # Create a temporary renumbered file
    tmpfile=$(mktemp)
    pdb_reres -$offset "$target" > "$tmpfile"

    # Run TMscore
    result=$(TMscore "$tmpfile" "$reference")

    # Clean up temporary file
    rm "$tmpfile"

    # Parse TM-score, RMSD, aligned residues
    tm_score=$(echo "$result" | grep -E "^TM-score\s*=" | awk '{print $3}')
    rmsd=$(echo "$result" | grep -E "RMSD of\s+the common residues" | sed -E 's/.*=\s*([0-9.]+).*/\1/')
    aligned=$(echo "$result" | grep "Number of residues in common=" | awk '{print $6}')

    # Fallbacks
    [[ -z "$tm_score" ]] && tm_score="NA"
    [[ -z "$rmsd" ]] && rmsd="NA"
    [[ -z "$aligned" ]] && aligned="NA"

    echo -e "${target}\t${tm_score}\t${aligned}\t${rmsd}" >> "$output_file"
done

# Extract top 20 by TM-score
tail -n +2 "$output_file" | sort -k2,2nr | head -n 20 | cut -f1 > "$top_list"

# Copy top 20 to receptor/top20
while read -r pdb; do
    cp "$pdb" "$top_dir/"
done < "$top_list"

echo "[INFO] Top 20 structures copied to $top_dir"

