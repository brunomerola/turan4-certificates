"""Independent checker: pure-Python exact value d(H) - sum_b <Q_b, M_b(H)> for a few graphs (cross-check of rv_eval.c).
Flag identity by an isomorphism test (search over bijections fixing the roots), not by canonical forms.
Usage: python rv_pycheck.py RUN_PREFIX MASK_COLEX [MASK_COLEX ...]
"""
import json
import sys
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial

import numpy as np


def colex_rank(e):
    return sum(comb(v, i + 1) for i, v in enumerate(sorted(e)))


def decode(mask, n):
    return {frozenset(e) for e in combinations(range(n), 4) if (mask >> colex_rank(e)) & 1}


def same_flag(A, B, s, m):
    for pu in permutations(range(s, m)):
        g = list(range(s)) + list(pu)
        if {frozenset(g[v] for v in e) for e in A} == B:
            return True
    return False


prefix = sys.argv[1]
cert = json.load(open(prefix + ".cert.json"))
Z = np.load(prefix + ".cert.npz")
M2 = cert["M"] ** 2
keys = cert["keys"]
Q = []
for k, key in enumerate(keys):
    A = Z[f"A{k}"].astype(object)
    Q.append(A.T @ A if A.shape[0] else None)
for hm in map(int, sys.argv[2:]):
    H = decode(hm, 7)
    tot = Fraction(len(H), 35)
    for k, key in enumerate(keys):
        if Q[k] is None or not any(Q[k].ravel()):
            continue
        s, m = key["s"], key["m"]
        sig = decode(key["sigma"], s)
        fl = [decode(f, m) for f in key["flags"]]
        conf = factorial(7) // factorial(7 - s) * comb(7 - s, m - s) * comb(7 - m, m - s)
        acc = 0
        for th in permutations(range(7), s):
            if {e for e in map(frozenset, combinations(range(s), 4)) if frozenset(th[v] for v in e) in H} != sig:
                continue
            rest = [v for v in range(7) if v not in th]
            idx = {}
            for U in combinations(rest, m - s):
                lab = list(th) + list(U)
                F = {e for e in map(frozenset, combinations(range(m), 4)) if frozenset(lab[v] for v in e) in H}
                idx[U] = next((i for i, f in enumerate(fl) if same_flag(F, f, s, m)), None)
            for U1 in idx:
                for U2 in idx:
                    if set(U1) & set(U2) or idx[U1] is None or idx[U2] is None:
                        continue
                    acc += int(Q[k][idx[U1], idx[U2]])
        tot -= Fraction(acc, conf * M2)
    print(hm, tot, float(tot), "== cert bound:", tot == Fraction(cert["bound"]))
