"""Round a float interior dual point to y = num / 2^K (num >= 1 on the whole support, sum num = 2^K exactly) and
write the JSON read by verify_dual.py.  Optionally mixes with another dual point (e.g. --mix FILE:eps).
Also reports the float scaled minimum eigenvalue (D^-1/2 Y D^-1/2) of every key at the rounded point.
Usage: python round_dual.py results/barrier_W.npz results/support_W.npz --K 72 --out results/dual_W.json
"""
import argparse
import json
from fractions import Fraction

import numpy as np

from barrier import Problem


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ynpz")
    ap.add_argument("support")
    ap.add_argument("--K", type=int, default=72)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    pb = Problem(a.support)
    src = np.load(a.ynpz)
    assert np.array_equal(src["masks"], pb.masks)
    y = np.asarray(src["y"], dtype=np.float64)
    assert (y > 0).all()
    den = 1 << a.K
    num = [max(1, int(Fraction(float(v)) * den)) for v in y]       # floor, at least 1
    diff = den - sum(num)
    i = int(np.argmax(y))
    num[i] += diff
    assert num[i] > 0 and sum(num) == den
    V = Fraction(sum(n * int(e) for n, e in zip(num, np.rint(pb.d * 35).astype(int))), 35 * den)
    yq = np.array([n / den for n in num])
    rel = pb.rel_min_eigs(yq)
    info = {"K": a.K, "support": len(num), "value": str(V), "value_float": float(V),
            "float_obj": float(pb.d @ y), "rel_min_eigs": rel, "min_y": float(y.min())}
    print(json.dumps(info), flush=True)
    json.dump({"den": den, "value": str(V), "source": a.ynpz,
               "support": [{"mask": int(mk), "num": n} for mk, n in zip(pb.masks.tolist(), num)]},
              open(a.out, "w"))


if __name__ == "__main__":
    main()
