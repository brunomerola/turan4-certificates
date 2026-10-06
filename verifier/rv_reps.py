"""Independent checker: own 6-vertex representatives, labelled counts and Burnside class counts of admissible 7-vertex 4-graphs.

1. All 2^15 4-graphs on {0..5}; canonical form = min over all 720 vertex permutations (brute force); admissible =
   every p-subset of {0..5} contains an edge (vacuous for p = 7).  reps = canonical forms of admissible graphs.
   Written as own 7-vertex masks (lex order of 4-subsets of {0..6}) to reps_p{p}.bin (int64).
2. Labelled count of admissible 7-vertex graphs:  L7 = sum_R (720/|Aut R|) * #{links L in [0,2^20): R+L admissible}.
   For p = 7 it must equal 2^35 - 1.
3. Burnside: #classes = (1/7!) sum_{g in S7} Fix(g); Fix(id) = L7; for g != id, Fix(g) = #admissible unions of
   <g>-orbits of 4-sets (enumerated, numpy).
Usage: python rv_reps.py P
"""
import json
import sys
import time
from itertools import combinations, permutations
from math import factorial

import numpy as np

p = int(sys.argv[1])
LAM = int(sys.argv[2]) if len(sys.argv) > 2 else 1   # every p-set spans >= LAM edges
suf = "" if LAM == 1 else f"_l{LAM}"
t0 = time.time()
E7 = list(combinations(range(7), 4))
I7 = {e: i for i, e in enumerate(E7)}
E6 = list(combinations(range(6), 4))
I6 = {e: i for i, e in enumerate(E6)}


_PC = np.array([bin(v).count("1") for v in range(1 << 16)], dtype=np.int8)


def popc15(a):
    return _PC[a & 0xFFFF]


def popc35(a):
    return _PC[a & 0xFFFF] + _PC[(a >> 16) & 0xFFFF] + _PC[(a >> 32) & 0xFFFF]


def psets(n):
    return [sum(1 << I7[e] for e in combinations(P, 4)) for P in combinations(range(n), p)]


# --- 1. 6-vertex reps
X = np.arange(1 << 15, dtype=np.int64)
canon = X.copy()
for g in permutations(range(6)):
    img = [I6[tuple(sorted(g[v] for v in e))] for e in E6]
    Y = np.zeros_like(X)
    for j in range(15):
        Y |= ((X >> j) & 1) << img[j]
    canon = np.minimum(canon, Y)
adm6 = np.ones(X.size, dtype=bool)
for P in combinations(range(6), p):
    m = sum(1 << I6[e] for e in combinations(P, 4))
    adm6 &= popc15(X & m) >= LAM
reps6 = np.unique(canon[adm6])
orbit = {int(r): int((canon == r).sum()) for r in reps6}
assert sum(orbit.values()) == int(adm6.sum())
to7 = [I7[e] for e in E6]           # 6-vertex own index -> 7-vertex own index


def lift(x):
    return sum(1 << to7[j] for j in range(15) if (x >> j) & 1)


reps7 = np.array([lift(int(r)) for r in reps6], dtype=np.int64)
reps7.tofile(f"reps_p{p}{suf}.bin")

# --- 2. labelled count
link = [I7[t + (6,)] for t in combinations(range(6), 3)]
Lr = np.arange(1 << 20, dtype=np.int64)
LH = np.zeros_like(Lr)
for j in range(20):
    LH |= ((Lr >> j) & 1) << link[j]
PS7 = psets(7)
per_rep = {}
L7 = 0
for r6, r7 in zip(reps6, reps7):
    H = LH | np.int64(r7)
    ok = np.ones(H.size, dtype=bool)
    for m in PS7:
        ok &= popc35(H & m) >= LAM
    c = int(ok.sum())
    per_rep[int(r7)] = c
    L7 += orbit[int(r6)] * c
raw_total = sum(per_rep.values())

# --- 3. Burnside
def cycle_type(g):
    seen, ct = set(), []
    for v in range(7):
        if v in seen:
            continue
        L, w = 0, v
        while w not in seen:
            seen.add(w); w = g[w]; L += 1
        ct.append(L)
    return tuple(sorted(ct, reverse=True))


classes = {}
for g in permutations(range(7)):
    ct = cycle_type(g)
    if ct not in classes:
        classes[ct] = [g, 0]
    classes[ct][1] += 1
fix = {}
for ct, (g, size) in classes.items():
    if ct == (1,) * 7:
        fix[ct] = L7
        continue
    gi = [I7[tuple(sorted(g[v] for v in e))] for e in E7]
    seen, orbs = set(), []
    for i in range(35):
        if i in seen:
            continue
        o, j = 0, i
        while j not in seen:
            seen.add(j); o |= 1 << j; j = gi[j]
        orbs.append(o)
    k = len(orbs)
    cnt = 0
    step = 1 << min(k, 22)
    for base in range(0, 1 << k, step):
        S = np.arange(base, base + step, dtype=np.int64)
        H = np.zeros_like(S)
        for j, o in enumerate(orbs):
            H |= np.where((S >> j) & 1, np.int64(o), np.int64(0))
        ok = np.ones(S.size, dtype=bool)
        for m in PS7:
            ok &= popc35(H & m) >= LAM
        cnt += int(ok.sum())
    fix[ct] = cnt
tot = sum(classes[ct][1] * fix[ct] for ct in classes)
assert tot % 5040 == 0
out = {"p": p, "lam": LAM, "n_reps6": int(reps6.size), "admissible_labelled_6": int(adm6.sum()), "raw_extensions": raw_total,
       "labelled_7": L7, "labelled_7_expected_p7": (1 << 35) - 1 if p == 7 else None,
       "burnside_classes": tot // 5040, "fix": {str(ct): [classes[ct][1], fix[ct]] for ct in classes},
       "per_rep_links": {str(k): v for k, v in per_rep.items()}, "time_s": time.time() - t0}
json.dump(out, open(f"reps_p{p}{suf}.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k not in ("per_rep_links",)}))
