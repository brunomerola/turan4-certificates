#!/usr/bin/env python3
"""Independent checker: machine check of the LIFTING used for the scope remark (reviewer code).

For a labelled type sigma on [s] and flag sizes s <= m1 < m2 with 2*m2 - s <= 7 (so that the top block (s, m2) is
evaluable on 7-vertex graphs), let D_{a,b}(H) be the exact DENSITY matrix of the configurations (theta, U1, U2):
theta a uniformly random injective s-tuple, U1, U2 disjoint uniformly random subsets of sizes a-s, b-s of the rest.
Let T[F, G] (F of size m1, G of size m2) = fraction of the (m1-s)-subsets W of the free part of G for which
G[[s] + W] is the flag F (an H-independent matrix).  Claimed identities, for EVERY 3-graph H on 7 vertices:
      D_{m1,m2}(H) = T D_{m2,m2}(H)        and        D_{m1,m1}(H) = T D_{m2,m2}(H) T^T.
Hence the mixed-size moment matrix of a type is T_all D_top(H) T_all^T with T_all = [T's; I], and a PSD Q on the
mixed flags is the PSD matrix T_all^T Q T_all on the top block, with the same <., M(H)> for every H (so the same slack
and the same "kills range M(phi)" condition).  m1 <= 2 and m1 = s (1 x 1, the type density) are included.
Checked exactly (Fractions) on 25 random 3-graphs of various densities and on 8 classes of phi; also checks that the
generic joint counts coincide with jw_lib.moment_counts on the 9 blocks (two codings of the moments)."""
from __future__ import annotations

import itertools
import math
import os
import pickle
import random
import sys
from collections import Counter, defaultdict
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import jw_lib as L  # noqa: E402

N = 7


def ftri(s, m):
    return [t for t in L.TRI[m] if t[2] >= s]


def code_of(E, th, U, s, m):
    best = None
    ft = ftri(s, m)
    for oU in itertools.permutations(U):
        w = tuple(th) + oU
        c = 0
        for i, (a, b, cc) in enumerate(ft):
            if E[w[a]][w[b]][w[cc]]:
                c |= 1 << i
        if best is None or c < best:
            best = c
    return best if best is not None else 0


def joint(mask, s, ma, mb):
    E = L.adj3(mask)
    out = defaultdict(Counter)
    for th in itertools.permutations(range(N), s):
        tm = 0
        for i, (a, b, c) in enumerate(L.TRI[s]):
            if E[th[a]][th[b]][th[c]]:
                tm |= 1 << i
        rest = [v for v in range(N) if v not in th]
        A = {U: code_of(E, th, U, s, ma) for U in itertools.combinations(rest, ma - s)}
        B = {U: code_of(E, th, U, s, mb) for U in itertools.combinations(rest, mb - s)}
        for U1, c1 in A.items():
            for U2, c2 in B.items():
                if set(U1).isdisjoint(U2):
                    out[tm][(c1, c2)] += 1
    return out


def nconf(s, ma, mb):
    return math.perm(N, s) * math.comb(N - s, ma - s) * math.comb(N - ma, mb - s)


def T_row(sig, s, m1, m2, G):
    """T[., G]: distribution of the m1-flag of a random (m1-s)-subset of the free part of the m2-flag G."""
    # rebuild the graph on [m2]: type triples + free-part triples from the code (min ordering = identity)
    E = [[[False] * m2 for _ in range(m2)] for _ in range(m2)]

    def add(t):
        for a, b, c in itertools.permutations(t):
            E[a][b][c] = True
    for i, t in enumerate(L.TRI[s]):
        if (sig >> i) & 1:
            add(t)
    for i, t in enumerate(ftri(s, m2)):
        if (G >> i) & 1:
            add(t)
    th = tuple(range(s))
    out = Counter()
    Ws = list(itertools.combinations(range(s, m2), m1 - s))
    for W in Ws:
        out[code_of(E, th, W, s, m1)] += Fraction(1, len(Ws))
    return out


def main():
    rnd = random.Random(7)
    tests = []
    for k in range(25):
        p = [0.1, 0.3, 0.5, 0.7, 0.9][k % 5]
        tests.append(sum(1 << i for i in range(35) if rnd.random() < p))
    phi = pickle.load(open(os.path.join(os.getcwd(), "local", "phi_Fanoc.pkl"), "rb"))["phi"]
    tests += sorted(phi, key=lambda H: -phi[H])[:8]
    triples = [(s, m1, m2) for s in range(0, 6) for m2 in range(max(3, s + 1), N + 1) if 2 * m2 - s <= N
               for m1 in range(s, m2)]
    print("checked (s, m1, m2):", triples)
    nchk = 0
    bad = 0
    # (a) two codings of the 9 blocks agree
    for H in tests[:6] + tests[-3:]:
        mc = L.moment_counts(H)
        for (s, m) in L.BLOCKS:
            jj = joint(H, s, m, m)
            for tm, cnt in jj.items():
                if mc.get((s, m, tm), Counter()) != cnt:
                    bad += 1
            if sum(len(v) for k, v in mc.items() if k[:2] == (s, m)) != sum(len(v) for v in jj.values()):
                bad += 1
    print("two codings of the 9 blocks agree:", bad == 0)
    # (b) lifting identities
    Tcache = {}
    for H in tests:
        for (s, m1, m2) in triples:
            Dt = joint(H, s, m2, m2)
            Dm = joint(H, s, m1, m2)
            D1 = joint(H, s, m1, m1)
            nt, nm, n1 = nconf(s, m2, m2), nconf(s, m1, m2), nconf(s, m1, m1)
            for sig in set(Dt) | set(Dm) | set(D1):
                top = {k: Fraction(v, nt) for k, v in Dt.get(sig, {}).items()}
                mix = {k: Fraction(v, nm) for k, v in Dm.get(sig, {}).items()}
                one = {k: Fraction(v, n1) for k, v in D1.get(sig, {}).items()}
                Gs = {g for kk in top for g in kk}
                for G in Gs:
                    if (sig, s, m1, m2, G) not in Tcache:
                        Tcache[(sig, s, m1, m2, G)] = T_row(sig, s, m1, m2, G)
                TM = defaultdict(Fraction)          # (T D_top)[F, G2]
                for (G1, G2), v in top.items():
                    for F, w in Tcache[(sig, s, m1, m2, G1)].items():
                        TM[(F, G2)] += w * v
                TMT = defaultdict(Fraction)         # (T D_top T^T)[F1, F2]
                for (F1, G2), v in TM.items():
                    for F2, w in Tcache[(sig, s, m1, m2, G2)].items():
                        TMT[(F1, F2)] += v * w
                for kk in set(TM) | set(mix):
                    nchk += 1
                    if TM.get(kk, 0) != mix.get(kk, 0):
                        bad += 1
                for kk in set(TMT) | set(one):
                    nchk += 1
                    if TMT.get(kk, 0) != one.get(kk, 0):
                        bad += 1
    print(f"lifting identities: {nchk} entries checked on {len(tests)} graphs, {bad} mismatches")
    print("LIFT CHECK", "OK" if bad == 0 else "FAILED")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
