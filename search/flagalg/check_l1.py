"""Independent check of a lottery local-packing certificate (lottery.py --cert) on concrete 3-graphs.

Only the certificate JSON is read; rooted flags are canonicalised here by brute force over orderings of the root and
of the extra vertices (pure Python, no engine code).  For each test hypergraph A (every 4-set contains an edge):
  (R) for every K_4^{(3)} K of A and every pair {u1, u2} outside K:  cpart(K; u1, u2) + W[link(u1), link(u2)] <= 1
      exactly (Fractions) -- the packing rows evaluated on configurations that actually occur;
  (Y) y_T = average over 2-sets U outside T of c(A[T + U], T), and S_K = sum_{T subset K} y_T  (must be
      <= 1 + O(1/n));
  (D) for A = shadow of a lottery system G: sum_{T in A} y_T <= sum_{K in G} S_K  (weak duality, exact) and the
      ratio to |G|.
Usage: python check_l1.py CERT.json [n_random]
"""
from __future__ import annotations

import json
import random
import sys
from fractions import Fraction
from itertools import combinations, permutations
from math import comb


def colex(e) -> int:
    return sum(comb(v, i + 1) for i, v in enumerate(sorted(e)))


def rooted_canon(A: set, T: tuple, U: tuple, r: int) -> int:
    m = len(T) + len(U)
    best = None
    pos_edges = list(combinations(range(m), r))
    idx = [colex(e) for e in pos_edges]
    for pt in permutations(T):
        for pu in permutations(U):
            order = pt + pu
            mask = 0
            for e, ix in zip(pos_edges, idx):
                if frozenset(order[i] for i in e) in A:
                    mask |= 1 << ix
            if best is None or mask < best:
                best = mask
    return best


def admissible(A: set, n: int, p: int, r: int) -> bool:
    return all(any(frozenset(e) in A for e in combinations(P, r)) for P in combinations(range(n), p))


def repair(A: set, n: int, p: int, r: int, rng) -> set:
    for P in combinations(range(n), p):
        if not any(frozenset(e) in A for e in combinations(P, r)):
            A.add(frozenset(rng.choice(list(combinations(P, r)))))
    return A


def turan_complement(parts) -> tuple:
    """Complement of Turan's K_4^{(3)}-free construction: triples inside a part, and 2 in V_i + 1 in V_{i-1}."""
    V = []
    for i, sz in enumerate(parts):
        V += [i] * sz
    n = len(V)
    A = set()
    for e in combinations(range(n), 3):
        cls = sorted(V[v] for v in e)
        cnt = {c: cls.count(c) for c in set(cls)}
        if len(cnt) == 1:
            A.add(frozenset(e))
        elif len(cnt) == 2:
            two = [c for c in cnt if cnt[c] == 2][0]
            one = [c for c in cnt if cnt[c] == 1][0]
            if one == (two - 1) % 3:
                A.add(frozenset(e))
    return A, n


def random_lottery_system(n: int, rng) -> list:
    """Random 4-sets added until every 4-set meets one in >= 3 points (k = p = 4, r = 3)."""
    G = []
    shadow = set()
    four = list(combinations(range(n), 4))
    rng.shuffle(four)
    for P in four:
        if not any(frozenset(e) in shadow for e in combinations(P, 3)):
            K = P
            G.append(K)
            for e in combinations(K, 3):
                shadow.add(frozenset(e))
    return G, shadow


class Cert:
    def __init__(self, path):
        c = json.load(open(path))
        self.k, self.r, self.p, self.m = c["k"], c["r"], c["p"], c["m"]
        assert (self.k, self.r, self.m) == (4, 3, 5) and c["tau_sos"], "checker written for k=4, r=3, m=5, tau-SOS"
        self.c = {int(f): Fraction(x) for f, x in zip(c["rooted_flags"], c["c_exact"])}
        den = Fraction(c["W_den"])
        self.W = [[Fraction(int(x)) / den for x in row] for row in c["W_num"]]
        self.link_index = {int(L): i for i, L in enumerate(c["tau1_links"])}
        self.b = Fraction(c["b_exact"])

    def cval(self, A, T, U):
        return self.c[rooted_canon(A, T, U, self.r)]

    def link(self, A, K, u):
        L = 0
        for i, j in combinations(range(4), 2):
            if frozenset((K[i], K[j], u)) in A:
                L |= 1 << colex((i, j))
        return self.link_index[L]


def check(cert: Cert, A: set, n: int, G=None) -> dict:
    assert admissible(A, n, 4, 3)
    V = range(n)
    y = {}
    for T in combinations(V, 3):
        if frozenset(T) in A:
            rest = [v for v in V if v not in T]
            Us = list(combinations(rest, 2))
            y[T] = sum(cert.cval(A, T, U) for U in Us) / len(Us)
    K4 = [K for K in combinations(V, 4) if all(frozenset(e) in A for e in combinations(K, 3))]
    max_row = Fraction(-10)
    max_SK = Fraction(-10)
    for K in K4:
        rest = [v for v in V if v not in K]
        links = {u: cert.link(A, K, u) for u in rest}
        for U in combinations(rest, 2):
            cp = sum(cert.cval(A, T, U) for T in combinations(K, 3))
            row = cp + cert.W[links[U[0]]][links[U[1]]]
            max_row = max(max_row, row)
        SK = sum(y[T] for T in combinations(K, 3))
        max_SK = max(max_SK, SK)
    out = {"n": n, "edges": len(A), "K4": len(K4), "max_row": float(max_row), "row_ok": max_row <= 1,
           "max_SK": float(max_SK), "sum_y": float(sum(y.values())), "sum_y_over_Cn3": float(sum(y.values()) / comb(n, 3))}
    if G is not None:
        SKs = {K: sum(y[T] for T in combinations(K, 3)) for K in G}
        out.update({"G": len(G), "duality_ok": sum(y.values()) <= sum(SKs.values()),
                    "sum_y_over_G": float(sum(y.values()) / len(G))})
    return out


def main():
    cert = Cert(sys.argv[1])
    nr = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    rng = random.Random(20261003)
    print(f"# certificate bound b = {float(cert.b):.10f}; c in [{float(min(cert.c.values())):.4f}, "
          f"{float(max(cert.c.values())):.4f}]", flush=True)
    allok = True
    for parts in [(4, 4, 4), (5, 4, 3), (3, 3, 3), (5, 5, 3)]:
        A, n = turan_complement(parts)
        res = check(cert, A, n)
        allok &= res["row_ok"]
        print("turan", parts, res, flush=True)
    for i in range(nr):
        n = rng.choice([9, 10, 11, 12])
        q = rng.choice([0.2, 0.35, 0.5, 0.65])
        A = {frozenset(e) for e in combinations(range(n), 3) if rng.random() < q}
        A = repair(A, n, 4, 3, rng)
        res = check(cert, A, n)
        allok &= res["row_ok"]
        print("random", q, res, flush=True)
    for i in range(nr):
        n = rng.choice([9, 10, 11, 12])
        G, A = random_lottery_system(n, rng)
        res = check(cert, A, n, G)
        allok &= res["row_ok"] and res["duality_ok"]
        print("lottery", res, flush=True)
    print("ALL ROW CHECKS OK" if allok else "ROW CHECK FAILURE", flush=True)


if __name__ == "__main__":
    main()
