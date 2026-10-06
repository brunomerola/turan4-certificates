"""Independent checker: exact re-check of the N = 7 certificates for 5-graphs (certificates/fivegraphs/cert_*_N7_all.json).

Own code; reads only the certificate JSON (blocks: s, m, sigma, flags (colex masks of 5-subsets), den, Qnum).
 * admissible(H): every p-set of [7] spans a number of edges in `allowed` (None = >= 1).
 * ALL 2^21 labelled 5-graphs on [7] are scanned; isomorphism classes found by orbit marking over all 5040
   permutations (numpy), so completeness is by exhaustion; labelled count = sum of orbit sizes is reported.
 * Q_sigma = Qnum / den: symmetric, PSD tested exactly (LDL^T over Fractions).
 * c(H) from the definition: every ordered injective theta with H[theta] == sigma (labelled), every ordered pair of
   disjoint (m-s)-sets, flag identified by canonical form (min colex mask over non-root permutations), divided by
   conf = 7!/(7-s)! C(7-s, m-s) C(7-m, m-s).  b = min_H [e(H)/21 - c(H)], exact.
Usage: python rv_r5.py CERT_JSON
"""
import json
import sys
import time
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial

import numpy as np

R, N = 5, 7
t0 = time.time()
cert = json.load(open(sys.argv[1]))
assert cert["r"] == R and cert["N"] == N
p = cert["p"]
allowed = cert["allowed"]
claimed = Fraction(cert["bound"])


def colex(e):
    return sum(comb(v, i + 1) for i, v in enumerate(sorted(e)))


def decode(mask, n):
    return {frozenset(e) for e in combinations(range(n), R) if (mask >> colex(e)) & 1}


def enc(es, n):
    return sum(1 << colex(e) for e in es)


E7 = list(combinations(range(N), R))          # own order (lex) for the scan
I7 = {frozenset(e): i for i, e in enumerate(E7)}
PS = [sum(1 << I7[frozenset(f)] for f in combinations(P, R)) for P in combinations(range(N), p)]
X = np.arange(1 << 21, dtype=np.int64)
PC = np.array([bin(v).count("1") for v in range(1 << 11)], dtype=np.int16)


def popc(a):
    return PC[a & 2047] + PC[(a >> 11) & 2047]


ok = np.ones(X.size, dtype=bool)
for m_ in PS:
    c = popc(X & m_)
    ok &= (c >= 1) if allowed is None else np.isin(c, allowed)
n_lab = int(ok.sum())
# orbit marking
perms = list(permutations(range(N)))
PM = np.array([[I7[frozenset(g[v] for v in e)] for e in E7] for g in perms], dtype=np.int64)   # 5040 x 21
seen = np.zeros(X.size, dtype=bool)
reps, orbit_sum = [], 0
for x in np.flatnonzero(ok):
    if seen[x]:
        continue
    bits = [(int(x) >> j) & 1 for j in range(21)]
    img = np.zeros(len(perms), dtype=np.int64)
    for j in range(21):
        if bits[j]:
            img |= np.int64(1) << PM[:, j]
    orb = np.unique(img)
    assert ok[orb].all(), "admissibility not invariant?!"
    seen[orb] = True
    reps.append(int(x))
    orbit_sum += orb.size
assert orbit_sum == n_lab
print(f"p={p} allowed={allowed}: labelled admissible {n_lab:,}; classes {len(reps)}; ({time.time()-t0:.0f}s)", flush=True)

# certificate
blocks = []
all_psd = True
for b in cert["blocks"]:
    s, m = b["s"], b["m"]
    assert 2 * m - s <= N
    den = int(b["den"])
    Q = [[Fraction(int(v), 1) for v in row] for row in b["Qnum"]]
    n = len(b["flags"])
    assert len(Q) == n and all(len(r) == n for r in Q)
    sym = all(Q[i][j] == Q[j][i] for i in range(n) for j in range(n))
    A = [r[:] for r in Q]
    psd = True
    for k in range(n):
        piv = A[k][k]
        if piv < 0 or (piv == 0 and any(A[k][j] != 0 for j in range(k, n))):
            psd = False
            break
        if piv == 0:
            continue
        for i in range(k + 1, n):
            f = A[i][k] / piv
            if f:
                for j in range(k, n):
                    A[i][j] -= f * A[k][j]
    all_psd &= sym and psd and den > 0
    if not any(int(v) for row in b["Qnum"] for v in row):
        continue
    sig = decode(b["sigma"], s)
    fl = {}
    for i, f in enumerate(b["flags"]):
        F = decode(f, m)
        assert {e for e in F if e <= frozenset(range(s))} == sig, "flag root part != sigma"
        cf = min(enc({frozenset(g[v] for v in e) for e in F}, m)
                 for g in (list(range(s)) + list(pu) for pu in permutations(range(s, m))))
        assert cf not in fl, "duplicate flag class"
        fl[cf] = i
    conf = factorial(N) // factorial(N - s) * comb(N - s, m - s) * comb(N - m, m - s)
    blocks.append((s, m, sig, fl, [[int(v) for v in row] for row in b["Qnum"]], den, conf))
print(f"Q blocks: {len(cert['blocks'])}, all symmetric PSD (exact): {all_psd}; nonzero blocks {[(x[0], x[1]) for x in blocks]}",
      flush=True)

best, arg, vals = None, None, []
for x in reps:
    H = {frozenset(E7[j]) for j in range(21) if (x >> j) & 1}
    c = Fraction(0)
    for s, m, sig, fl, Qn, den, conf in blocks:
        acc = 0
        Es = [frozenset(e) for e in combinations(range(s), R)]
        Em = [frozenset(e) for e in combinations(range(m), R)]
        for th in permutations(range(N), s):
            if {e for e in Es if frozenset(th[v] for v in e) in H} != sig:
                continue
            rest = [v for v in range(N) if v not in th]
            idx = {}
            for U in combinations(rest, m - s):
                lab = list(th) + list(U)
                F = {e for e in Em if frozenset(lab[v] for v in e) in H}
                cf = min(enc({frozenset(g[v] for v in e) for e in F}, m)
                         for g in (list(range(s)) + list(pu) for pu in permutations(range(s, m))))
                idx[U] = fl.get(cf)
            for U1, a in idx.items():
                for U2, bb in idx.items():
                    if a is not None and bb is not None and not set(U1) & set(U2):
                        acc += Qn[a][bb]
        c += Fraction(acc, den * conf)
    v = Fraction(len(H), comb(N, R)) - c
    vals.append(v)
    if best is None or v < best:
        best, arg = v, x
print(f"exact min over all {len(reps)} classes: {best} = {float(best):.15f}; argmin own mask {arg} "
      f"({bin(arg).count('1')} edges)")
print(f"claimed {claimed} = {float(claimed):.15f}; EQUAL: {best == claimed}; PSD: {all_psd}")
print(f"pi <= 1 - b = {float(1 - best):.15f}; ({time.time()-t0:.0f}s)")
