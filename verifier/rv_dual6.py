"""Independent checker: exact check of the producer's N = 6 dual points (own code; reads only the JSON support lists).

Claim checked: no plain flag-algebra certificate at N = 6 (any blocks (s, m), m > s >= 0, 2m - s <= 6, any types,
any PSD Q on any flag subsets) proves more than V = sum_H y_H d(H), where y >= 0, sum y = 1 is supported on
admissible 6-vertex 4-graphs and every moment matrix Y_sigma = sum_H y_H M_sigma(H) is PSD.
Weak duality: b <= min_H [d(H) - c_Q(H)] <= sum_H y_H [d(H) - c_Q(H)] = V - sum <Q_sigma, Y_sigma> <= V.
M_sigma(H)[a, b] = (number of (theta ordered, U1, U2 ordered disjoint) with H[theta] = sigma exactly, flag(theta,U1)=a,
flag(theta,U2)=b) / conf (all theta counted).  Types: one per isomorphism class (canonical = min mask over S_s);
flags: classes of admissible labelled extensions modulo non-root permutations.  PSD test: exact LDL^T over Fractions
with zero-pivot handling.
Usage: python rv_dual6.py JSON P LAM
"""
import json
import sys
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial

R, N = 4, 6
fn, p, lam = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
d = json.load(open(fn))
sup = d["dual"]["support"]
claimed = Fraction(d["dual"]["value"])


def adm(es, n):
    return all(sum(frozenset(f) in es for f in combinations(P, R)) >= lam for P in combinations(range(n), p))


def mask(es, n):
    E = list(combinations(range(n), R))
    return sum(1 << i for i, e in enumerate(E) if frozenset(e) in es)


def canon_fix(es, s, m):
    return min(mask({frozenset(g[v] for v in e) for e in es}, m)
               for g in (list(range(s)) + list(pu) for pu in permutations(range(s, m))))


Hs = []
for t in sup:
    es = {frozenset(e) for e in t["graph_edges"]}
    assert len(es) == t["edges"] and all(len(e) == 4 and max(e) < N for e in es)
    Hs.append((Fraction(t["y"]), es))
ys = [y for y, _ in Hs]
assert all(y > 0 for y in ys) and sum(ys) == 1
assert all(adm(es, N) for _, es in Hs), "support graph not admissible"
V = sum(y * Fraction(len(es), comb(N, R)) for y, es in Hs)
print(f"support {len(Hs)} graphs, sum y = 1, all admissible; V = sum y d = {V}; claimed {claimed}; equal {V == claimed}")


def ldl_psd(A):
    n = len(A)
    A = [row[:] for row in A]
    for k in range(n):
        piv = A[k][k]
        if piv < 0:
            return False
        if piv == 0:
            if any(A[k][j] != 0 for j in range(k, n)):
                return False
            continue
        for i in range(k + 1, n):
            f = A[i][k] / piv
            if f:
                for j in range(k, n):
                    A[i][j] -= f * A[k][j]
    return True


blocks = [(s, m) for s in range(0, N) for m in range(s + 1, N + 1) if 2 * m - s <= N]
ok_all = True
nmat = 0
for s, m in blocks:
    Es = list(combinations(range(s), R))
    types = sorted({min(mask({frozenset(g[v] for v in e) for e in sig}, s) for g in permutations(range(s)))
                    for sig in ({frozenset(Es[i]) for i in range(len(Es)) if (x >> i) & 1}
                                for x in range(1 << len(Es)))
                    if adm(sig, s)})
    conf = factorial(N) // factorial(N - s) * comb(N - s, m - s) * comb(N - m, m - s)
    for tm in types:
        sig = {frozenset(Es[i]) for i in range(len(Es)) if (tm >> i) & 1}
        Y = {}
        for y, H in Hs:
            for th in permutations(range(N), s):
                if {e for e in map(frozenset, Es) if frozenset(th[v] for v in e) in H} != sig:
                    continue
                rest = [v for v in range(N) if v not in th]
                fl = {}
                for U in combinations(rest, m - s):
                    lab = list(th) + list(U)
                    F = {e for e in map(frozenset, combinations(range(m), R)) if frozenset(lab[v] for v in e) in H}
                    fl[U] = canon_fix(F, s, m)
                for U1 in fl:
                    for U2 in fl:
                        if set(U1) & set(U2):
                            continue
                        Y[(fl[U1], fl[U2])] = Y.get((fl[U1], fl[U2]), 0) + y / conf
        idx = sorted({a for a, _ in Y} | {b for _, b in Y})
        A = [[Y.get((a, b), Fraction(0)) for b in idx] for a in idx]
        sym = all(A[i][j] == A[j][i] for i in range(len(idx)) for j in range(len(idx)))
        psd = ldl_psd(A)
        nmat += 1
        ok_all &= sym and psd
        print(f"  block (s={s}, m={m}) type {tm}: Y {len(idx)}x{len(idx)} nonzero support, symmetric {sym}, PSD {psd}")
print(f"blocks {blocks}; {nmat} moment matrices; ALL PSD: {ok_all}; => no plain N = 6 certificate exceeds {V}")
