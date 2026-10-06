"""Exact 7-vertex distribution of Giraud's construction (complement form, t(5,4) <= 5/16), as a dual point.

Large-n limit of Giraud's construction (see ../n8est/giraud_plugin.py): each vertex is a row or a column (prob. 1/2),
matrix bits i.i.d. fair; complement edges = 4-sets with 0 or 4 rows (all-row / all-column) and 2+2 sets whose 2x2
submatrix has even sum.  Enumerates all 2^7 side vectors and all 2^(|rows| |cols|) bit patterns (330,626
configurations, weights 2^-7 2^-(|rows||cols|)), groups the labelled masks by canonical form, and writes a dual-point
JSON (den = 2^19) for verify_dual.py.  A genuine limit object: every moment matrix is E[f f^T] >= 0, so the
verifier must PASS with value exactly 5/16 (test of the verifier; also an exact interior-type ingredient).
Usage: python giraud_exact.py results/giraud_dual.json
"""
import json
import sys
from collections import Counter
from fractions import Fraction
from itertools import combinations

import numpy as np

import paths  # noqa: F401
from hypergraphs import canon_batch

N = 7
E7 = sorted(combinations(range(N), 4), key=lambda e: e[::-1])


def main():
    out = sys.argv[1]
    cnt = Counter()
    for side in range(1 << N):
        rows = [v for v in range(N) if (side >> v) & 1]
        cols = [v for v in range(N) if not (side >> v) & 1]
        pairs = [(r, c) for r in rows for c in cols]
        pidx = {p: i for i, p in enumerate(pairs)}
        nb = len(pairs)
        w = 1 << (12 - nb)                        # weight * 2^19 = 2^(19 - 7 - nb)
        for bits in range(1 << nb):
            m = 0
            for j, e in enumerate(E7):
                rr = [v for v in e if v in rows]
                if len(rr) in (0, 4):
                    m |= 1 << j
                elif len(rr) == 2:
                    cc = [v for v in e if v not in rr]
                    par = sum((bits >> pidx[(r, c)]) & 1 for r in rr for c in cc) & 1
                    if par == 0:
                        m |= 1 << j
            cnt[m] += w
    lab = np.array(list(cnt.keys()), dtype=np.int64)
    wts = np.array([cnt[int(x)] for x in lab], dtype=object)
    can = canon_batch(lab, N, 4, kind="all")
    cls = Counter()
    for c, wv in zip(can.tolist(), wts):
        cls[int(c)] += int(wv)
    den = 1 << 19
    assert sum(cls.values()) == den
    V = Fraction(sum(v * bin(k).count("1") for k, v in cls.items()), 35 * den)
    print(f"labelled masks {lab.size}, classes {len(cls)}, value {V} = {float(V)}")
    json.dump({"den": den, "value": str(V), "source": "Giraud construction, exact 7-vertex distribution",
               "support": [{"mask": k, "num": v} for k, v in sorted(cls.items())]}, open(out, "w"))


if __name__ == "__main__":
    main()
