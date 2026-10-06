"""Independent checker: exhaustive computer parts of Theorem 2 (local-to-global), own code, own algorithm.

C_n: every locally Giraud 4-graph H on n vertices (every 7-vertex induced subgraph is a Giraud system) is a Giraud
system.  Method (different from the producer's: no oddness system, no linear algebra):
  H locally Giraud on range(n)  =>  R := H[range(n-1)] is locally Giraud, hence Giraud (C_{n-1}; n-1 = 7: by
  definition).  The link L of the last vertex v = n-1 is a 3-graph on range(n-1).  H is locally Giraud iff R is
  and, for EVERY 6-subset W of range(n-1), the 7-vertex graph H[W + v] lies in G_7 (the 2,404 labelled Giraud
  systems on 7 vertices, own enumeration r_cf.py).  H[W + v] consists of R[W] (15 bits) and L restricted to the
  20 triples of W, so each W allows an explicit finite list of patterns for L|W.  All L are obtained by an exact
  hash join of these lists over all C(n-1, 6) sets W (every W is processed, also after all triples are covered).
  Every survivor is then tested for being a Giraud system.
Runs:
  n = 8 : over ALL 2,404 labelled R in G_7 (no isomorphism reduction).  Survivors must be exactly G_8.
  n = 9 : over ALL 32,981 labelled R in G_8.  Survivors must be exactly G_9 (604,426, own enumeration).
  n = 10: over the 33 orbit representatives of G_9 (orbit peeling in r_cf.py).  Survivors tested with two
          independent Giraud tests (st_lib.is_giraud_forced and st_lib.is_giraud_struct) and compared, per R, with
          the set of Giraud extensions of R constructed directly (new vertex added to either side, new row/column
          arbitrary).
Usage: python r_lg.py {8|9|10}
"""
from __future__ import annotations

import json
import os
import sys
import time
from itertools import combinations
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.getcwd()   # outputs (results_*.json, local/) go to the working directory
import st_lib as L  # noqa: E402


def load_set(name):
    return [int(x) for x in open(os.path.join(OUT, "local", name)).read().split()]


def build_tables(m, G7):
    """W-list in join order, and for each W: (care mask, quad indices of R inside W, dict low15 -> list of global
    patterns)."""
    T3 = L.subsets(m, 3)
    t3i = {t: i for i, t in enumerate(T3)}
    loc3 = L.subsets(6, 3)
    loc4 = L.subsets(6, 4)
    by_low = {}
    for g in G7:
        by_low.setdefault(g & ((1 << 15) - 1), []).append(g >> 15)
    Ws = list(combinations(range(m), 6))
    # greedy order: start with (0..5); then always the W with the most already-covered triples (ties: lex)
    trip = {W: frozenset(t3i[tuple(W[j] for j in t)] for t in loc3) for W in Ws}
    order = [Ws[0]]
    covered = set(trip[Ws[0]])
    rest = set(Ws[1:])
    while rest:
        W = max(sorted(rest), key=lambda w: len(trip[w] & covered))
        order.append(W)
        covered |= trip[W]
        rest.remove(W)
    tabs = []
    for W in order:
        g_of_local = [t3i[tuple(W[j] for j in t)] for t in loc3]       # local triple j -> global triple index
        care = 0
        for gi in g_of_local:
            care |= 1 << gi
        qidx = [L.cidx(tuple(W[j] for j in q)) for q in loc4]
        pats = {}
        for low, highs in by_low.items():
            lst = []
            for h in highs:
                v = 0
                for j in range(20):
                    if (h >> j) & 1:
                        v |= 1 << g_of_local[j]
                lst.append(v)
            pats[low] = lst
        tabs.append((care, qidx, pats))
    return tabs, len(T3)


def extensions(R, m, tabs):
    vals = [0]
    C = 0
    for care, qidx, pats in tabs:
        low = 0
        for j, qi in enumerate(qidx):
            if (R >> qi) & 1:
                low |= 1 << j
        lst = pats.get(low)
        if not lst:
            return []
        sh = C & care
        idx = {}
        for w in lst:
            idx.setdefault(w & sh, []).append(w)
        new = set()
        for p in vals:
            for w in idx.get(p & sh, ()):
                new.add(p | w)
        vals = list(new)
        C |= care
        if not vals:
            return []
    n4 = comb(m, 4)
    return [R | (Lk << n4) for Lk in vals]


def direct_giraud_extensions(R, m):
    """All Giraud systems on range(m+1) whose restriction to range(m) is R, built from the representations of R."""
    n = m + 1
    out = set()
    for A in L.is_giraud_forced(R, m):
        A = sorted(A)
        B = [v for v in range(m) if v not in A]
        M = {}
        if B:
            a0, b0 = max(A), max(B)
            for a in A:
                for b in B:
                    if a != a0 and b != b0:
                        M[(a, b)] = 1 - L.bit(R, (a0, a, b0, b))
        for side in (0, 1):
            A2 = A + [m] if side == 0 else A
            B2 = B if side == 0 else B + [m]
            k = len(B) if side == 0 else len(A)          # new row (side A) or new column (side B), arbitrary
            for bits in range(1 << k):
                M2 = dict(M)
                for i, u in enumerate(B if side == 0 else A):
                    key = (m, u) if side == 0 else (u, m)
                    M2[key] = (bits >> i) & 1
                g = L.giraud(n, A2, M2)
                assert g & ((1 << comb(m, 4)) - 1) == R
                out.add(g)
    return out


def main():
    n = int(sys.argv[1])
    m = n - 1
    t0 = time.time()
    G7 = load_set("G7.txt")
    assert len(G7) == 2404
    tabs, nt = build_tables(m, G7)
    res = {"n": n, "W_sets": len(tabs), "triples": nt}
    if n in (8, 9):
        Rs = load_set(f"G{m}.txt")
        Gn = set(load_set(f"G{n}.txt"))
        low = (1 << comb(m, 4)) - 1
        ref = {}
        for g in Gn:
            ref.setdefault(g & low, set()).add(g)
        allsurv = 0
        ok = True
        nonG = []
        maxpart = 0
        for i, R in enumerate(Rs):
            S = set(extensions(R, m, tabs))
            allsurv += len(S)
            eq = S == ref.get(R, set())
            ok &= eq
            if not eq:
                nonG.extend(str(h) for h in S - Gn)
            if i % 2000 == 0:
                print(i, len(S), allsurv, round(time.time() - t0, 1), flush=True)
        res.update({"labelled_R": len(Rs), "survivors_total": allsurv, "G_n_size": len(Gn),
                    "survivors_equal_giraud_extensions_for_every_R": ok, "non_giraud_survivors": nonG[:20],
                    f"C{n}_PASS": bool(ok and allsurv == len(Gn))})
    else:
        reps = json.load(open(os.path.join(OUT, "local", "reps9.json")))
        assert len(reps) == 33 and sum(s for _, s in reps) == 604426
        per = []
        ok = True
        tot = 0
        for i, (R, osz) in enumerate(reps):
            S = extensions(R, m, tabs)
            assert len(S) == len(set(S))
            g1 = [bool(L.is_giraud_forced(h, n)) for h in S]
            g2 = [bool(L.is_giraud_struct(h, n)) for h in S]
            D = direct_giraud_extensions(R, m)
            eq = set(S) == D
            allg = all(g1) and all(g2)
            ok &= allg and eq and g1 == g2
            tot += len(S)
            per.append({"R": str(R), "orbit_size": osz, "locally_giraud_extensions": len(S),
                        "all_giraud_forced_test": all(g1), "all_giraud_struct_test": all(g2),
                        "direct_giraud_extensions": len(D), "equal_sets": eq})
            print(i, len(S), len(D), allg, eq, round(time.time() - t0, 1), flush=True)
        res.update({"reps": per, "locally_giraud_total": tot, "C10_PASS": bool(ok)})
    res["time_s"] = round(time.time() - t0, 1)
    json.dump(res, open(os.path.join(OUT, f"results_lg_n{n}.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "reps"}))


if __name__ == "__main__":
    main()
