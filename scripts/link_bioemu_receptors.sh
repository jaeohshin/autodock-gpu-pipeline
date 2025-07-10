for i in $(seq 1 50); do
    src="/data/work/flowpacker/samples/abl1/run_1/frame_re_${i}.pdb"
    dst="/data/work/dock/virtual_screening/input/receptors/ABL1/receptor_$(printf "%04d" $i).pdb"

    if [ -f "$src" ]; then
        ln -s "$src" "$dst"
    else
        echo "[WARN] Skipping missing: $src"
    fi
done
