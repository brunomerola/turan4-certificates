"""Cross-check of two independent moment-matrix implementations on the same rational dual point:
verify_dual.py (pure Python, own enumeration) against the engine tables (prep_support.py / barrier.Problem, built
from ../n7/common.py + kernels.py).  Moment matrices of the same (block, type) must be permutation-similar: equal
number of occurring flags, equal traces and equal sorted spectra (float).
Usage: python crosscheck_engine.py results/dual_W_final.json results/support_W.npz
"""
import json
import sys
from fractions import Fraction

import numpy as np

import paths  # noqa: F401
from barrier import Problem
from verify_dual import BLOCKS, BlockTables, decode, moment_counts, to_matrix


def main():
    d = json.load(open(sys.argv[1]))
    pb = Problem(sys.argv[2])
    data = np.load(sys.argv[2])
    den = int(d["den"])
    ymap = {int(r["mask"]): int(r["num"]) for r in d["support"]}
    assert set(ymap) == set(int(x) for x in pb.masks)
    y = np.array([ymap[int(x)] / den for x in pb.masks])
    tabs = [BlockTables(s, m) for s, m in BLOCKS]
    Yv = moment_counts([(ymap[int(x)], decode(int(x))) for x in pb.masks], tabs)
    conf = {(t.s, t.m): t.conf for t in tabs}
    ver = []
    for key in sorted(Yv):
        A, flags = to_matrix(Yv[key])
        sc = Fraction(1, den * conf[key[:2]])
        Af = np.array([[float(Fraction(v) * sc) for v in row] for row in A])
        ver.append((key, Af))
    # engine keys in model order (only keys with live flags are in pb.keys)
    kk = [k for k in range(int(data["nkeys"])) if data[f"live{k}"].size]
    eng = [(k, K.Y(y)) for k, K in zip(kk, pb.keys)]
    out = []
    used = set()
    for k, Ye in eng:
        we = np.linalg.eigvalsh(Ye)
        best = None
        for j, (key, Af) in enumerate(ver):
            if j in used or Af.shape != Ye.shape:
                continue
            wv = np.linalg.eigvalsh(Af)
            err = float(np.abs(wv - we).max() / max(np.abs(we).max(), 1e-300))
            if best is None or err < best[0]:
                best = (err, j, key, float(np.trace(Af)), float(np.trace(Ye)))
        assert best is not None, f"no verifier matrix of size {Ye.shape} for engine key {k}"
        used.add(best[1])
        out.append({"engine_key": k, "verifier_block_type": list(best[2]), "dim": Ye.shape[0],
                    "spectrum_rel_err": best[0], "trace_engine": best[4], "trace_verifier": best[3]})
        print(json.dumps(out[-1]), flush=True)
        assert best[0] < 1e-9
    assert len(used) == len(ver), "verifier has (block, type) pairs the engine does not"
    print("CROSSCHECK PASS: all moment matrices permutation-similar (spectra agree)")


if __name__ == "__main__":
    main()
