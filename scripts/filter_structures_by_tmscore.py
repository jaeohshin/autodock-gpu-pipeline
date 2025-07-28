## Calculate TMscore of receptors
# Run
# python tmscore_calculate.py . ../receptor_crystal/


import subprocess
import glob
import os
import argparse
import csv

def run_tmscore(model_path, ref_path, tmscore_exec="TMscore"):
    cmd = [tmscore_exec, model_path, ref_path, "-seq"]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True)
    for line in result.stdout.splitlines():
        if line.strip().startswith("TM-score"):
            tokens = line.split()
            if '=' in tokens:
                eq_index = tokens.index('=')
                score_str = tokens[eq_index + 1]
                return float(score_str)
            else:
                return float(tokens[2])
    return None

def collect_tmscores(models_dir, ref_path, tmscore_exec="TMscore"):
    models = sorted(glob.glob(os.path.join(models_dir, "*.pdb")))
    results = []
    for model in models:
        score = run_tmscore(model, ref_path, tmscore_exec)
        if score is not None:
            results.append((os.path.basename(model), score))
    results.sort(key=lambda x: x[1], reverse=True)
    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Grouped TM-score ranking report to CSV.")
    parser.add_argument("parent_dir", help="Parent directory containing kinase subfolders.")
    parser.add_argument("crystal_dir", help="Directory with reference crystal subfolders.")
    parser.add_argument("--tmscore_exec", default="TMscore", help="Path to TMscore executable.")
    parser.add_argument("--output_csv", default="all_tmscores.csv", help="Output CSV file.")
    args = parser.parse_args()

    subdirs = [d for d in sorted(os.listdir(args.parent_dir)) if os.path.isdir(os.path.join(args.parent_dir, d))]
    all_results = []

    for subdir in subdirs:
        models_dir = os.path.join(args.parent_dir, subdir)
        ref_path = os.path.join(args.crystal_dir, subdir, "receptor.pdb")
        if not os.path.exists(ref_path):
            print(f"Skipping {subdir}: reference PDB not found.")
            continue
        kinase_results = collect_tmscores(models_dir, ref_path, args.tmscore_exec)
        if not kinase_results:
            print(f"{subdir}: no valid scores.")
            continue
        for model, score in kinase_results:
            all_results.append({
                "kinase": subdir,
                "receptor_name": model,
                "tmscore": score
            })

    # Write CSV output
    with open(args.output_csv, "w", newline="") as csvfile:
        fieldnames = ["kinase", "receptor_name", "tmscore"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_results)

    print(f"CSV report saved to {args.output_csv}")
