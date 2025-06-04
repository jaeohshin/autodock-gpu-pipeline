import os
import glob
from datetime import datetime

ROOT_DIR = "../cognate"
timestamp = datetime.now().strftime("%Y-%m-%d_%H%M")
OUT_CSV = os.path.join(ROOT_DIR, f"avg_reference_rmsd_{timestamp}.csv")

def extract_reference_rmsds(dlg_file):
    rmsds = []
    in_table = False
    with open(dlg_file, 'r') as f:
        for line in f:
            if "Rank | Sub-" in line and "Reference" in line:
                in_table = True
                continue
            if in_table:
                if "RANKING" in line:
                    parts = line.strip().split()
                    try:
                        rmsd = float(parts[5])  # 6th column
                        rmsds.append(rmsd)
                    except (IndexError, ValueError):
                        continue
                elif line.strip() == "":
                    break
    return rmsds

def main():
    dlg_paths = sorted(glob.glob(os.path.join(ROOT_DIR, "*/output/ligand_docked.dlg")))
    print(f"[INFO] Found {len(dlg_paths)} docking result files.")

    with open(OUT_CSV, "w") as fout:
        fout.write("Kinase,NumRuns,AvgReferenceRMSD,GoodPose_AvgRMSD,GoodPose_MinRMSD\n")
        for dlg_path in dlg_paths:
            kinase = os.path.basename(os.path.dirname(os.path.dirname(dlg_path)))
            try:
                rmsds = extract_reference_rmsds(dlg_path)
                if not rmsds:
                    print(f"[WARN] No RMSDs found in {kinase}")
                    continue
                avg_rmsd = sum(rmsds) / len(rmsds)
                min_rmsd = min(rmsds)
                good_avg = 1 if avg_rmsd < 2.0 else 0
                good_min = 1 if min_rmsd < 2.0 else 0
                fout.write(f"{kinase},{len(rmsds)},{avg_rmsd:.3f},{good_avg},{good_min}\n")
                print(f"[DONE] {kinase}: Avg={avg_rmsd:.3f}, Min={min_rmsd:.3f}, GoodAvg={good_avg}, GoodMin={good_min}")
            except Exception as e:
                print(f"[ERROR] {kinase}: {e}")

    print(f"[INFO] Results saved to: {OUT_CSV}")

if __name__ == "__main__":
    main()
