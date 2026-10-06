"""Scale-up estimate for N = 7 (r = 3, every p-set contains an edge): number of admissible graphs and SDP block sizes.

Labelled admissible graphs on 7 vertices = sum over the classes C on 6 vertices of (720/|Aut C|) * #{links L on the
6 old vertices : C + L admissible}  (each labelled 7-vertex graph is counted once via its restriction to {0..5}).
The number of isomorphism classes is >= labelled/7! (equality iff no graph has a nontrivial automorphism).
Block sizes for the standard N = 7 types (s, m) = (1, 4), (3, 5), (5, 6): number of sigma-flags, computed exactly.
Usage: python estimate_n7.py [p]
"""
import sys
import time
from math import comb, factorial

import numpy as np

from flags import flag_family
from hypergraphs import EveryPSetHasEdge, bits_matrix, enumerate_classes, perm_weight_matrix

p = int(sys.argv[1]) if len(sys.argv) > 1 else 4
r = 3
t0 = time.time()
pred = EveryPSetHasEdge(p)
cls = enumerate_classes(6, r, pred)
reps = cls[6]
W = perm_weight_matrix(6, r, "all", 0)
aut = (bits_matrix(reps, comb(6, r)) @ W.T == reps[:, None].astype(float)).sum(axis=1)
links = np.arange(1 << comb(6, r - 1), dtype=np.int64) << comb(6, r)
labelled = 0
valid_total = 0
for i, R in enumerate(reps):
    ok = pred(np.int64(R) | links, 7, r)
    v = int(ok.sum())
    valid_total += v
    labelled += (factorial(6) // int(aut[i])) * v
print(f"p={p}: classes on 6 vertices {reps.size}; augmentation candidates passing the predicate {valid_total:,}")
print(f"labelled admissible 3-graphs on 7 vertices: {labelled:,}; iso classes >= {labelled / factorial(7):,.0f}")
for s, m in [(1, 4), (3, 5), (5, 6)]:
    sizes = [flag_family(s, m, int(sig), r, pred).n for sig in cls[s]]
    print(f"  types s={s} (m={m}): {len(sizes)} types, block sizes {sizes}")
print(f"time {time.time() - t0:.1f}s")
