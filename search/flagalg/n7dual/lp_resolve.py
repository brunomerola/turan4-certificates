"""Step 0 of the N = 7 dual-point search: re-solve the frozen-cut restricted LP of the converged p = 5 run (working set W, cuts V*) and
store its primal pseudo-density y (float) together with the per-key moment-matrix terms of every graph of W.

Float only (HiGHS); nothing here is certified.  Imports highspy (through lpcg7); never import ortools here.
Usage: python lp_resolve.py [--out results/lp_y.npz]
"""
import argparse
import json
import time
from types import SimpleNamespace

import numpy as np

import paths  # noqa: F401  (sys.path)
from common import build_model
from lpcg7 import ALL_BLOCKS, NFIX, LPCG7


def lp_args(**kw):
    a = dict(lp_tol=1e-9, lp_threads=1, lp_solver="choose", simplex_strategy=-1, row_scale=1, keep_cuts=10 ** 9,
             prune_age=10, tau_max=0.0, b_cap=2e-3, nchunks=8, hot_size=1000, rows_per_price=0, price_tol=1e-9)
    a.update(kw)
    return SimpleNamespace(**a)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", default=paths.STATE_CONV)
    ap.add_argument("--out", default=paths.RES + "/lp_y.npz")
    a = ap.parse_args()
    t0 = time.time()
    model = build_model(5, ALL_BLOCKS)
    st = np.load(a.state)
    W = st["W"].astype(np.int64)
    cg = LPCG7(model, W, 0.25, lp_args())
    cg.add_graphs(W)
    items = []
    for k in range(len(model.keys)):
        items += [(k, v) for v in st[f"V{k}"]]
    cg.add_cuts(items)
    cg.set_stationarity(False)
    print(f"# setup {time.time() - t0:.1f}s: graphs {len(cg.W)}, cuts {len(cg.cuts)}", flush=True)
    t = time.time()
    cg.h.run()
    sol = cg.h.getSolution()
    yall = np.array(sol.col_value)
    rd = np.array(sol.row_dual)
    y = np.maximum(yall[NFIX:], 0.0)
    b_lp = float(rd[0])
    lam = np.maximum(rd[NFIX:], 0.0)
    print(f"# LP {time.time() - t:.1f}s: b_lp = {b_lp:.14f}, sum y = {y.sum():.15f}, y.d = "
          f"{float(y @ (cg.e / 35.0)):.14f}, support {int((y > 0).sum())}, z = {yall[0]:.3e}", flush=True)
    Ys = cg.moments(y)
    rep = []
    for k, Y in enumerate(Ys):
        w = np.linalg.eigvalsh(Y)
        T = cg.T[k]
        occ = np.zeros(Y.shape[0], dtype=bool)
        P = cg.P[k]
        for g in range(P.shape[0]):
            occ[P[g, T[1]]] = True
            occ[P[g, T[2]]] = True
        diag0 = np.abs(np.diag(Y)) < 1e-300
        rep.append({"key": k, "block": [model.blocks[model.keys[k][0]].s, model.blocks[model.keys[k][0]].m],
                    "dim": int(Y.shape[0]), "flags_occurring_in_W": int(occ.sum()),
                    "zero_diag_at_y": int(diag0.sum()), "min_eig": float(w[0]), "max_eig": float(w[-1]),
                    "n_neg_1e-12": int((w < -1e-12).sum()), "n_abs_le_1e-12": int((np.abs(w) <= 1e-12).sum()),
                    "n_pos_1e-12_1e-9": int(((w > 1e-12) & (w <= 1e-9)).sum()),
                    "n_gt_1e-9": int((w > 1e-9).sum()), "n_cuts": int(sum(1 for c in cg.cuts if c[0] == k))})
        print(json.dumps(rep[-1]), flush=True)
    np.savez(a.out, W=W, y=y, b_lp=b_lp, lam=lam, e=cg.e, cdd=cg.cdd,
             **{f"T{k}_{j}": cg.T[k][j] for k in range(len(model.keys)) for j in range(4)})
    with open(a.out.replace(".npz", ".json"), "w") as f:
        json.dump({"b_lp": b_lp, "support": int((y > 0).sum()), "keys": rep, "time_s": time.time() - t0}, f,
                  indent=1)
    print(f"# done {time.time() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
