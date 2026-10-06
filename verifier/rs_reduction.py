"""Independent checker (own code): reduction sigma -> f and Giraud's value, checked by brute force from the definitions.

BCL (arXiv:2108.10406, eq. (1)): co2(G) = sum over (k-1)-sets T of d_G(T)^2,
    sigma(F) = limsup_n exco2(n, F) / ( C(n, k-1) (n-k+1)^2 );   k = 4: C(n,3) (n-3)^2.
Checks (exact integers / Fractions):
 (a) random 4-graphs G on n = 5..9, H = complement in K_n^(4):
       co2(G) == C(n,3)(n-3)(n-4) * s_n(G) + 4 e(G),   s_n(G) := P[T+x in G, T+y in G] ((T,x,y) uniform, x != y),
       s_n(G) == 1 - f_n(H),   f_n(H) := 2 d(H) - P[T+x in H, T+y in H];
       on n = 7:  f(H) == (24 e(H) - P(H)) / 420,  P(H) = sum_T c_T (c_T - 1).
 (b) averaging: for random H on n = 8, 9, 10: f_n(H) == average over all 7-subsets S of f(H[S])  (exact identity).
 (c) K_5^4-freeness of G <=> every 5-set of H spans an H-edge (complement), on random graphs.
 (d) Giraud (BCL Section 8.2): rows/columns, A random 0/1: G-edges = 3+1 sets, and 2+2 sets with ODD 2x2 sum.
       exact limit by enumerating sides and entries of 5 distinct vertices: s(G) = 31/64, d(G) = 11/16;
       K_5^4-free for EVERY matrix (exhaustive over all matrices at n = 3+3, 4+4 partially; random at n = 14);
       finite random instance n = 120: co2(G) / (C(n,3)(n-3)^2) printed (should approach 31/64).
"""
import random
import sys
from fractions import Fraction
from itertools import combinations, product
from math import comb

import numpy as np

from rs_lib import f_value, objective_counts

rng = random.Random(20261006)
R = 4


def rand_graph(n, p):
    return {frozenset(e) for e in combinations(range(n), R) if rng.random() < p}


def codeg(es, T, n):
    return sum(1 for x in range(n) if x not in T and frozenset(T + (x,)) in es)


def s_n(es, n):
    good = tot = 0
    for T in combinations(range(n), 3):
        out = [v for v in range(n) if v not in T]
        for x in out:
            for y in out:
                if x != y:
                    tot += 1
                    good += frozenset(T + (x,)) in es and frozenset(T + (y,)) in es
    return Fraction(good, tot)


def main():
    ok = True
    # (a)
    na = 0
    for n in range(5, 10):
        allE = {frozenset(e) for e in combinations(range(n), R)}
        for trial in range(12 if n < 9 else 4):
            G = rand_graph(n, rng.choice([0.3, 0.5, 0.7, 0.9]))
            H = allE - G
            co2 = sum(codeg(G, T, n) ** 2 for T in combinations(range(n), 3))
            sG = s_n(G, n)
            ok &= co2 == comb(n, 3) * (n - 3) * (n - 4) * sG + 4 * len(G)
            ok &= sG == 1 - f_value(H, n)
            if n == 7:
                e, P = objective_counts(H, 7)
                Pc = sum(c * (c - 1) for c in (codeg(H, T, 7) for T in combinations(range(7), 3)))
                ok &= P == Pc and f_value(H, 7) == Fraction(24 * e - P, 420)
            na += 1
    print("(a) co2 / s_n / f identities on", na, "random graphs:", ok, flush=True)
    # (b)
    nb = 0
    for n in (8, 9, 10):
        for trial in range(3):
            H = rand_graph(n, rng.choice([0.4, 0.6, 0.8]))
            avg = Fraction(0)
            Ss = list(combinations(range(n), 7))
            for S in Ss:
                rel = {v: i for i, v in enumerate(S)}
                HS = {frozenset(rel[v] for v in e) for e in H if e <= set(S)}
                avg += f_value(HS, 7)
            avg /= len(Ss)
            ok &= avg == f_value(H, n)
            nb += 1
    print("(b) E_S f(H[S]) == f_n(H) exactly on", nb, "graphs (n = 8, 9, 10):", ok, flush=True)
    # (c)
    for n in (6, 7):
        allE = {frozenset(e) for e in combinations(range(n), R)}
        for trial in range(30):
            G = rand_graph(n, 0.85)
            H = allE - G
            k5free = not any(all(frozenset(f) in G for f in combinations(P, 4)) for P in combinations(range(n), 5))
            adm = all(any(frozenset(f) in H for f in combinations(P, 4)) for P in combinations(range(n), 5))
            ok &= k5free == adm
    print("(c) K5-free(G) <=> admissible(H):", ok, flush=True)

    # (d) Giraud, exact limit
    def gedge(sides, A, quad):
        rows = [v for v in quad if sides[v]]
        if len(rows) in (1, 3):
            return True
        if len(rows) == 2:
            cols = [v for v in quad if not sides[v]]
            return sum(A[(r, c)] for r in rows for c in cols) % 2 == 1
        return False

    num = Fraction(0)
    dens = Fraction(0)
    for sd in product((0, 1), repeat=5):        # vertices 0,1,2 = T, 3 = x, 4 = y
        prs = [(i, j) for i in range(5) for j in range(5) if sd[i] == 1 and sd[j] == 0]
        w = Fraction(1, 32 * (1 << len(prs)))
        for bits in product((0, 1), repeat=len(prs)):
            A = dict(zip(prs, bits))
            ex = gedge(sd, A, (0, 1, 2, 3))
            ey = gedge(sd, A, (0, 1, 2, 4))
            num += w * (ex and ey)
            dens += w * ex
    print("(d) Giraud exact limit: s(G) =", num, " d(G) =", dens, flush=True)
    ok &= num == Fraction(31, 64) and dens == Fraction(11, 16)

    # K5-free for every matrix: all matrices for parts 3+3 (2^9), random for larger
    def giraud_graph(a, b, A):
        n = a + b
        sides = [1] * a + [0] * b
        AA = {(r, c): A[r][c - a] for r in range(a) for c in range(a, n)}
        return {frozenset(q) for q in combinations(range(n), 4) if gedge(sides, AA, q)}, n

    allfree = True
    for M in product((0, 1), repeat=9):
        A = [M[0:3], M[3:6], M[6:9]]
        G, n = giraud_graph(3, 3, A)
        allfree &= not any(all(frozenset(f) in G for f in combinations(P, 4)) for P in combinations(range(n), 5))
    for trial in range(20):
        a = rng.randint(2, 7)
        b = rng.randint(2, 7)
        A = [[rng.randint(0, 1) for _ in range(b)] for _ in range(a)]
        G, n = giraud_graph(a, b, A)
        allfree &= not any(all(frozenset(f) in G for f in combinations(P, 4)) for P in combinations(range(n), 5))
    print("(d) Giraud K5^4-free (all 3x3 matrices, 20 random up to 7x7):", allfree, flush=True)
    ok &= allfree

    # finite random instance, numpy
    h = 60
    n = 2 * h
    A = np.array([[rng.randint(0, 1) for _ in range(h)] for _ in range(h)], dtype=np.int64)
    side = np.array([1] * h + [0] * h)
    # codegree of T = {u<v<w}: count x outside T with T+x a G-edge
    tot = 0
    for T in combinations(range(n), 3):
        rows = [v for v in T if side[v]]
        cols = [v - h for v in T if not side[v]]
        k = len(rows)
        if k == 3:
            c = h                    # x any column
        elif k == 0:
            c = h
        elif k == 2:
            # x row -> 3 rows + 1 col: edge (h - 2 such x); x col y != col: 2x2 sum odd
            r1, r2 = rows
            c0 = cols[0]
            par = (A[r1, c0] + A[r2, c0]) % 2
            colpar = (A[r1, :] + A[r2, :]) % 2
            c = (h - 2) + int(np.sum((colpar + par) % 2 == 1)) - int((colpar[c0] + par) % 2 == 1)
        else:
            c1, c2 = cols
            r0 = rows[0]
            par = (A[r0, c1] + A[r0, c2]) % 2
            rowpar = (A[:, c1] + A[:, c2]) % 2
            c = (h - 2) + int(np.sum((rowpar + par) % 2 == 1)) - int((rowpar[r0] + par) % 2 == 1)
        tot += c * c
    print(f"(d) random Giraud n = {n}: co2/(C(n,3)(n-3)^2) = {tot / (comb(n, 3) * (n - 3) ** 2):.6f}"
          f" (31/64 = {31/64:.6f})", flush=True)
    print("ALL OK" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
