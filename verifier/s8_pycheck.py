"""Independent pure-Python definitional check of the slack at given 7-vertex graphs (independent of s8_prep's tables and
of s8_eval.c; pure Python ints and Fractions, no numpy arithmetic).

usage: s8_pycheck.py PROBLEM CERT_PREFIX MASKS.txt OUT.json
For every mask H in MASKS.txt (one integer per line):
  * admissibility by the definition (no labelled copy of a forbidden graph / every copy meets H);
  * f(H) from the explicit (T, x, y) loop (s8_lib.objective_value);
  * for every block (s, m) of the certificate: types = increasing canonical masks of admissible s-vertex graphs that
    have a flag; flags = sorted canonical (root-fixing) masks; then the sum over EVERY injective root tuple theta of
    [7] and EVERY ordered pair (U1, U2) of disjoint (m-s)-subsets of the rest, with H[theta] equal to a type sigma as
    a labelled graph, of Q_sigma[flag(theta, U1), flag(theta, U2)] / #configurations, with Q = A^T A / M^2 computed
    entrywise on demand in Python ints from the .cert.npz rows;
  * slack = f - sum, an exact Fraction.
Writes {mask: slack (string), slack * DEN} for comparison with s8_eval's per-mask output.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial

import numpy as np

import s8_lib as L

N = 7


def flags_for(prob, s, m):
    r = L.PROBLEMS[prob][0]
    nsig = comb(s, r)
    nfree = comb(m, r) - nsig
    tps = sorted({L.canon_all(c, s, r) for c in range(1 << nsig) if L.admissible(prob, c, s)}) if nsig else [0]
    types, flags = [], []
    for sig in tps:
        fl = sorted({L.canon_fix(sig | (fc << nsig), m, r, s) for fc in range(1 << nfree)
                     if L.admissible(prob, sig | (fc << nsig), m)})
        if fl:
            types.append(sig)
            flags.append({f: i for i, f in enumerate(fl)})
    return types, flags


def main():
    prob, cert, mfile, outp = sys.argv[1:5]
    r, mode, keys, obj = L.PROBLEMS[prob]
    J = json.load(open(cert + ".cert.json"))
    Z = np.load(cert + ".cert.npz")
    A = [[[int(x) for x in row] for row in Z[f"A{k}"].tolist()] for k in range(len(J["keys"]))]
    M2 = int(J["M"]) ** 2
    blocks = [tuple(b) for b in J["blocks"]]
    model = []
    k = 0
    for (s, m) in blocks:
        types, flags = flags_for(prob, s, m)
        kk = list(range(k, k + len(types)))
        for t, kidx in zip(types, kk):
            assert J["keys"][kidx]["sigma"] == t and J["keys"][kidx]["n_flags"] == len(flags[types.index(t)])
        k += len(types)
        model.append((s, m, types, flags, kk))
    assert k == len(J["keys"])
    gcache = {}

    def Gent(key, a, b):
        x = (key, a, b) if a <= b else (key, b, a)
        if x not in gcache:
            gcache[x] = sum(row[a] * row[b] for row in A[key])
        return gcache[x]

    masks = [int(l.split()[0]) for l in open(mfile) if l.strip()]
    out = {}
    for H in masks:
        assert L.admissible(prob, H, N), H
        E = L.edge_set(H, N, r)
        f = L.objective_value(obj, H, N)
        total = Fraction(0)
        for (s, m, types, flags, kk) in model:
            tset = {t: i for i, t in enumerate(types)}
            denom = factorial(N) // factorial(N - s) * comb(N - s, m - s) * comb(N - m, m - s)
            acc = 0
            nconf = 0
            for theta in permutations(range(N), s):
                code = sum(1 << L.cidx(e) for e in L.edges(s, r) if frozenset(theta[i] for i in e) in E)
                rest = [v for v in range(N) if v not in theta]
                Us = list(combinations(rest, m - s))
                if code not in tset:
                    nconf += sum(1 for U1 in Us for U2 in Us if not set(U1) & set(U2))
                    continue
                ti = tset[code]
                fl = {}
                for U in Us:
                    V = tuple(theta) + U
                    fm = sum(1 << L.cidx(e) for e in L.edges(m, r) if frozenset(V[i] for i in e) in E)
                    assert fm & ((1 << comb(s, r)) - 1) == code
                    fl[U] = flags[ti][L.canon_fix(fm, m, r, s)]
                for U1 in Us:
                    for U2 in Us:
                        if set(U1) & set(U2):
                            continue
                        nconf += 1
                        acc += Gent(kk[ti], fl[U1], fl[U2])
            assert nconf == denom
            total += Fraction(acc, denom * M2)
        slack = f - total
        out[str(H)] = {"slack": str(slack), "S": str(slack * L.DEN7[obj] * 5040 * M2)}
        print(H, float(slack), flush=True)
    json.dump(out, open(outp, "w"), indent=1)


if __name__ == "__main__":
    main()
