"""Exploration (float): frozen-cut LP value restricted to a subset of W (e.g. the alive graphs after facial
reduction), plus a few rounds of new eigenvector cuts on that subset.  Imports highspy (via lpcg7), no ortools.
Usage: python lp_subset.py results/support_W.npz [--rounds 0]
"""
import argparse
import json
import time

import numpy as np
import scipy.linalg

import paths
from common import build_model
from lp_resolve import lp_args
from lpcg7 import ALL_BLOCKS, NFIX, LPCG7


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("support")
    ap.add_argument("--rounds", type=int, default=0)
    ap.add_argument("--max-cuts", type=int, default=60)
    a = ap.parse_args()
    model = build_model(5, ALL_BLOCKS)
    S = np.load(a.support)["masks"]
    st = np.load(paths.STATE_CONV)
    cg = LPCG7(model, S, 0.25, lp_args())
    cg.add_graphs(S)
    items = []
    for k in range(len(model.keys)):
        items += [(k, v) for v in st[f"V{k}"]]
    cg.add_cuts(items)
    cg.set_stationarity(False)
    for r in range(a.rounds + 1):
        t = time.time()
        cg.h.run()
        sol = cg.h.getSolution()
        y = np.maximum(np.array(sol.col_value)[NFIX:], 0.0)
        b = float(np.array(sol.row_dual)[0])
        Ys = cg.moments(y)
        cand, mins = [], []
        for k, Y in enumerate(Ys):
            w, V = scipy.linalg.eigh(Y)
            mins.append(float(w[0]))
            cand += [(float(w[i]), k, V[:, i]) for i in range(min(12, w.size)) if w[i] < -1e-13]
        cand.sort(key=lambda z: z[0])
        print(json.dumps({"round": r, "graphs": len(cg.W), "cuts": len(cg.cuts), "b_lp": b,
                          "obj": float(y @ (cg.e / 35.0)), "support": int((y > 0).sum()), "min_eig": min(mins),
                          "lp_s": round(time.time() - t, 1)}), flush=True)
        if r < a.rounds:
            cg.add_cuts([(k, v / np.linalg.norm(v)) for _, k, v in cand[:a.max_cuts]])
    np.savez(a.support.replace(".npz", "_lpy.npz"), y=y, masks=np.array(cg.W, dtype=np.int64), b=b)


if __name__ == "__main__":
    main()
