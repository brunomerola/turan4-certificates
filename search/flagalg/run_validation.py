"""Run the validation ladder (each configuration in its own process) and write results/validation_table.{txt,json}.

Usage (heavy rows take ~1 min; run under reproduce/run_capped.sh):  python run_validation.py
"""
import json
import os
import subprocess
import sys
import time

PY = sys.executable
HERE = os.path.dirname(os.path.abspath(__file__))
CONFIGS = [
    # label, args, literature value
    ("V1 Mantel", "--r 2 --p 3 --N 3 --sharp --target 1/2", "1/2"),
    ("V1 Mantel", "--r 2 --p 3 --N 4 --sharp --target 1/2", "1/2"),
    ("V1 Mantel", "--r 2 --p 3 --N 5 --sharp --target 1/2", "1/2"),
    ("V2 Turan K4", "--r 2 --p 4 --N 4 --sharp --target 1/3", "1/3"),
    ("V2 Turan K4", "--r 2 --p 4 --N 5 --sharp --target 1/3", "1/3"),
    ("V2 Turan K4", "--r 2 --p 4 --N 6 --sharp --target 1/3", "1/3"),
    ("V3 K4^3", "--r 3 --p 4 --N 5", "(none published)"),
    ("V3 K4^3", "--r 3 --p 4 --N 6", "0.438334 (1-0.561666, Razborov 2010)"),
    ("K4^3+ind E1", "--r 3 --p 4 --N 6 --allowed 1,2,4 --sharp --target 4/9", "4/9 (1-5/9, Razborov 2010)"),
    ("t(5,4)", "--r 4 --p 5 --N 6 --sharp --target 1/4", ">=0.2636 known (627/2379)"),
    ("t(6,4)", "--r 4 --p 6 --N 6 --sharp --target 1/10", ">=0.10324 known (51/494)"),
    ("no-SOS check", "--r 3 --p 4 --N 6 --types none", "T(6,4,3)/20 = 3/10"),
]


def main():
    rows = []
    for label, args, lit in CONFIGS:
        tag = args.replace("--", "").replace(" ", "_").replace("/", "-").replace(",", "")
        out = os.path.join(HERE, "results", f"val_{tag}.json")
        t = time.time()
        subprocess.run([PY, os.path.join(HERE, "turan.py")] + args.split() + ["--out", out], check=True,
                       stdout=subprocess.DEVNULL, cwd=HERE)
        wall = time.time() - t
        d = json.load(open(out))
        sh = d.get("sharp_rounding") or {}
        rows.append({
            "label": label, "r": d["r"], "p": d["p"], "N": d["N"], "n_admissible": d["n_admissible"],
            "blocks": [b["n_flags"] for b in d["blocks"]], "float_bound": d["sdp"]["b"],
            "solver_status": d["sdp"]["status"], "dual_obj": -d["sdp"]["obj_dual"],
            "exact_eig": d["eig_rounding"]["b_exact_float"],
            "exact_sharp": sh.get("b_exact") if sh.get("success") else None,
            "time_s": round(wall, 2), "peak_rss_mb": round(d["peak_rss_mb"], 1), "literature": lit})
        print(rows[-1], flush=True)
    with open(os.path.join(HERE, "results", "validation_table.json"), "w") as f:
        json.dump(rows, f, indent=1)
    hdr = (f"{'problem':14s} {'r':>2s} {'p':>2s} {'N':>2s} {'#adm':>5s} {'float bound':>15s} {'exact (eig)':>15s} "
           f"{'exact sharp':>11s} {'time s':>7s} {'RSS MB':>7s}  blocks / literature")
    lines = [hdr, "-" * len(hdr)]
    for r in rows:
        lines.append(f"{r['label']:14s} {r['r']:2d} {r['p']:2d} {r['N']:2d} {r['n_admissible']:5d} "
                     f"{r['float_bound']:15.12f} {r['exact_eig']:15.12f} {str(r['exact_sharp'] or '-'):>11s} "
                     f"{r['time_s']:7.2f} {r['peak_rss_mb']:7.1f}  {r['blocks']} / {r['literature']}")
    txt = "\n".join(lines)
    with open(os.path.join(HERE, "results", "validation_table.txt"), "w") as f:
        f.write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
