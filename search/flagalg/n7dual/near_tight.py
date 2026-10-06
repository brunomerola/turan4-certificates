"""Float slack of every admissible class (3,908,438) w.r.t. a certificate (b, Q): candidates for the dual support.
Certificate from a lpcg7 .best.npz (factor rows F_k, Q_k = F_k^T F_k).  Writes the masks with
slack(H) = d(H) - <Q, M(H)> <= b* + delta for several delta, sorted by slack.
Usage: python near_tight.py PREFIX.best.npz --out results/near_tight.npz
"""
import argparse
import json

import numpy as np

import paths
from common import build_model
from kernels import qsum_flat, scan_float, tables


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("best")
    ap.add_argument("--out", required=True)
    ap.add_argument("--keep", type=int, default=200000)
    a = ap.parse_args()
    src = np.load(a.best)
    blocks = [tuple(int(x) for x in b) for b in src["blocks"]]
    model = build_model(int(src["p"]), blocks)
    Qs = [src[f"F{k}"].T @ src[f"F{k}"] for k in range(len(model.keys))]
    qs = qsum_flat(model, Qs)
    cls = np.load(paths.CLASSES_P5).astype(np.int64)
    sl, bad = scan_float(cls, *tables(model), qs, 0.0, 0.0, 0.0, 8)
    assert not bad.any()
    bstar = float(sl.min())
    order = np.argsort(sl)[:a.keep]
    info = {"b_star": bstar, "n": int(cls.size)}
    for dl in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2):
        info[f"n_within_{dl:g}"] = int((sl <= bstar + dl).sum())
    print(json.dumps(info), flush=True)
    np.savez(a.out, masks=cls[order], slack=sl[order], b_star=bstar)


if __name__ == "__main__":
    main()
