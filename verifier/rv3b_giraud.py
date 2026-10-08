"""Independent checker: own exact checks for Proposition lotup(a), c(6,4,5) <= 7/16 via Kahn 1996 on the Giraud host.

Host G_d on A u B, A = B = F_2^d (N = 2^d, n = 2N).  Edges: 4-subsets of A, of B; cross {a,a',b,b'} with
<a+a', b+b'> = 0.  Hyperedges of X: internal 6-sets (15 internal 4-subsets), gauge-trivial (3,3) 6-sets
(all 9 cross 4-subsets are G-edges -- tested DIRECTLY on the 9 subsets, not via the span argument).

(T) Turan (5,4) by brute force, d = 2, 3, 4 (all 5-sets).
(C) for EVERY cross edge, the number of (3,3) 6-sets containing it whose 9 cross 4-subsets are all edges;
    d = 3, 4, 5 (numpy, exhaustive); compare with (N/2-2)(N/4-2); also |E_cross| vs (N-1)(N/2)^2(N/2-1).
(I) internal: C(N-4,2) (trivial; checked by count at d = 3, 4).
(P) pair weight: max number of cross hyperedges through two distinct cross edges (d = 4 exhaustive over all
    hyperedges; d = 5 via the 5-point unions), compared with the hand bound N.
(V) exact value 15 t(X)/C(2N,4) as a rational function of N; leading coefficients -> 7/16; values d = 3..30.
(E) end-to-end at d = 4: the 6-sets of X form an (n,6,4,5)-lottery system (all 5-sets of 32 points).
"""
import itertools
import math
import sys
from fractions import Fraction as Fr

import numpy as np


def par(v):
    return bin(v).count("1") & 1


def dotmat(N):
    D = np.zeros((N, N), dtype=np.uint8)
    for u in range(N):
        for w in range(N):
            D[u, w] = par(u & w)
    return D


def is_edge(Q, N, D):
    As = [v for v in Q if v < N]
    Bs = [v - N for v in Q if v >= N]
    if len(As) == 4 or len(Bs) == 4:
        return True
    if len(As) == 2:
        return D[As[0] ^ As[1], Bs[0] ^ Bs[1]] == 0
    return False


def turan(d):
    N = 1 << d
    D = dotmat(N)
    for S in itertools.combinations(range(2 * N), 5):
        if not any(is_edge(Q, N, D) for Q in itertools.combinations(S, 4)):
            return False, S
    return True, None


def cross_edges(N, D):
    E = []
    for a, a2 in itertools.combinations(range(N), 2):
        for b, b2 in itertools.combinations(range(N), 2):
            if D[a ^ a2, b ^ b2] == 0:
                E.append((a, a2, b, b2))
    return np.array(E, dtype=np.int64)


def cross_counts(d, chunk=4096):
    """Exhaustive: for each cross edge, count (a3, b3) with a3 not in {a,a2}, b3 not in {b,b2} such that all 9
    cross 4-subsets of {a,a2,a3} x {b,b2,b3} (2+2 choices) are edges."""
    N = 1 << d
    D = dotmat(N)
    E = cross_edges(N, D)
    ar = np.arange(N)
    counts = np.zeros(len(E), dtype=np.int64)
    for s in range(0, len(E), chunk):
        e = E[s:s + chunk]
        a, a2, b, b2 = e[:, 0:1], e[:, 1:2], e[:, 2:3], e[:, 3:4]       # (c,1)
        A3 = ar[None, :]                                                    # (1,N)
        # A-pairs: (a,a2), (a,a3), (a2,a3); B-pairs: (b,b2), (b,b3), (b2,b3)
        Apairs = [np.broadcast_to(a ^ a2, (len(e), N)), a ^ A3, a2 ^ A3]   # each (c,N) over a3
        Bpairs = [np.broadcast_to(b ^ b2, (len(e), N)), b ^ A3, b2 ^ A3]   # each (c,N) over b3
        ok = np.ones((len(e), N, N), dtype=bool)                          # [edge, a3, b3]
        for Ap in Apairs:
            for Bp in Bpairs:
                ok &= (D[Ap[:, :, None], Bp[:, None, :]] == 0)
        # exclude a3 in {a,a2}, b3 in {b,b2}
        ok &= (A3 != a)[:, :, None] & (A3 != a2)[:, :, None]
        ok &= (A3 != b)[:, None, :] & (A3 != b2)[:, None, :]
        counts[s:s + chunk] = ok.reshape(len(e), -1).sum(axis=1)
    return len(E), np.unique(counts)


def hyperedges_cross(N, D):
    H = []
    for A in itertools.combinations(range(N), 3):
        for B in itertools.combinations(range(N), 3):
            if all(D[A[i] ^ A[j], B[k] ^ B[l]] == 0 for i, j in ((0, 1), (0, 2), (1, 2))
                   for k, l in ((0, 1), (0, 2), (1, 2))):
                H.append((A, B))
    return H


def pair_weight_d4():
    from collections import Counter
    N = 16
    D = dotmat(N)
    H = hyperedges_cross(N, D)
    pc = Counter()
    for A, B in H:
        es = [(p, q) for p in itertools.combinations(A, 2) for q in itertools.combinations(B, 2)]
        for e1, e2 in itertools.combinations(es, 2):
            pc[(e1, e2)] += 1
    return len(H), max(pc.values())


def pair_weight_5sets(d):
    """Max over (3,2) point sets {a1,a2,a3,b1,b2} contained in >= 1 cross hyperedge of the number of b3 completing
    it to a gauge-trivial 6-set (two cross edges with a 5-point union lie in exactly that many hyperedges; with a
    6-point union in exactly 1).  By the A<->B symmetry of the host, (2,3) sets give the same maximum."""
    N = 1 << d
    D = dotmat(N).astype(bool)
    best = 0
    ar = np.arange(N)
    for A in itertools.combinations(range(N), 3):
        u = [A[0] ^ A[1], A[0] ^ A[2], A[1] ^ A[2]]
        # B pairs (b1,b2) valid: all <u_i, b1+b2> = 0
        for b1 in range(N):
            for b2 in range(b1 + 1, N):
                w = b1 ^ b2
                if D[u[0], w] or D[u[1], w] or D[u[2], w]:
                    continue
                # b3 completions
                m = np.ones(N, dtype=bool)
                for ui in u:
                    m &= ~D[ui, b1 ^ ar]
                    m &= ~D[ui, b2 ^ ar]
                m[b1] = m[b2] = False
                c = int(m.sum())
                if c > best:
                    best = c
    return best


class Poly:
    """Polynomial in N with Fraction coefficients, c[i] * N^i."""

    def __init__(self, c):
        self.c = [Fr(v) for v in c]

    def __add__(self, o):
        n = max(len(self.c), len(o.c))
        return Poly([(self.c[i] if i < len(self.c) else 0) + (o.c[i] if i < len(o.c) else 0) for i in range(n)])

    def __mul__(self, o):
        if not isinstance(o, Poly):
            return Poly([v * o for v in self.c])
        r = [Fr(0)] * (len(self.c) + len(o.c) - 1)
        for i, a in enumerate(self.c):
            for j, b in enumerate(o.c):
                r[i + j] += a * b
        return Poly(r)

    def ev(self, N):
        return sum((v * Fr(N) ** i for i, v in enumerate(self.c)), Fr(0))


def lin(a, b):  # a*N + b
    return Poly([b, a])


def value_poly():
    # |E_cross| = (N-1)(N/2)^2(N/2-1); t(X) = |E_cross|/9 + 2 C(N,4)/15; c_N = 15 t(X)/C(2N,4)
    Ecross = lin(1, -1) * lin(Fr(1, 2), 0) * lin(Fr(1, 2), 0) * lin(Fr(1, 2), -1)
    CN4 = lin(1, 0) * lin(1, -1) * lin(1, -2) * lin(1, -3) * Fr(1, 24)
    C2N4 = lin(2, 0) * lin(2, -1) * lin(2, -2) * lin(2, -3) * Fr(1, 24)
    tX = Ecross * Fr(1, 9) + CN4 * Fr(2, 15)
    num = tX * 15
    return num, C2N4, Ecross


if __name__ == "__main__":
    which = sys.argv[1:] or ["T", "C", "P", "V", "E"]
    if "T" in which:
        for d in (2, 3, 4):
            ok, S = turan(d)
            print("(T) d=%d n=%d: Turan(5,4) by brute force: %s %s" % (d, 2 << d, ok, S or ""), flush=True)
    if "C" in which:
        for d in (3, 4, 5):
            N = 1 << d
            ne, cs = cross_counts(d)
            print("(C) d=%d: |E_cross| = %d (formula %d) ; counts over ALL cross edges: %s ; formula (N/2-2)(N/4-2) = %d"
                  % (d, ne, (N - 1) * (N // 2) ** 2 * (N // 2 - 1), cs.tolist(), (N // 2 - 2) * (N // 4 - 2)),
                  flush=True)
        for d in (3, 4):
            N = 1 << d
            # internal 4-set in C(N-4,2) internal 6-sets: count directly for one 4-set
            Q = set(range(4))
            cnt = sum(1 for S in itertools.combinations(range(N), 6) if Q <= set(S))
            print("(I) d=%d: internal count %d, C(N-4,2) = %d" % (d, cnt, math.comb(N - 4, 2)))
    if "P" in which:
        nH, mx = pair_weight_d4()
        N = 16
        cnt = (N // 2 - 2) * (N // 4 - 2)
        print("(P) d=4: %d cross hyperedges (= |E_cross|*cnt/9 = %d) ; max # through two cross edges = %d ; "
              "alpha_cross = %s ; internal (N-5)/C(N-4,2) = %s" % (nH, 6720 * cnt // 9, mx, Fr(mx, cnt),
                                                                 Fr(N - 5, math.comb(N - 4, 2))), flush=True)
        for d in (4, 5):
            N = 1 << d
            b = pair_weight_5sets(d)
            cnt = (N // 2 - 2) * (N // 4 - 2)
            print("(P) d=%d: max # cross hyperedges through a (3,2) 5-set = %d (hand bound N-2 = %d) ; "
                  "alpha_cross = max(that,1)/cnt = %s" % (d, b, N - 2, Fr(max(b, 1), cnt)), flush=True)
    if "V" in which:
        num, den, Ec = value_poly()
        print("(V) numerator 15 t(X) coeffs (N^0..N^4): %s" % [str(v) for v in num.c])
        print("    C(2N,4) coeffs: %s" % [str(v) for v in den.c])
        print("    leading ratio = %s" % (num.c[4] / den.c[4]))
        for d in range(3, 31, 3):
            N = 1 << d
            v = num.ev(N) / den.ev(N)
            print("    d=%2d c_N = %.15f  diff from 7/16 = %+.3e" % (d, float(v), float(v - Fr(7, 16))))
        # check the polynomial formula against direct integer formula at d=4,5
        for d in (4, 5):
            N = 1 << d
            tX = Fr((N - 1) * (N // 2) ** 2 * (N // 2 - 1), 9) + Fr(2 * math.comb(N, 4), 15)
            assert 15 * tX / math.comb(2 * N, 4) == num.ev(N) / den.ev(N)
        print("    poly == direct formula at d=4,5: True")
    if "E" in which:
        d = 4
        N = 16
        D = dotmat(N)
        H = hyperedges_cross(N, D)
        # covered 4-sets: all internal 4-sets (internal 6-sets exist since N >= 6) + cross 4-subsets of H
        cov = set()
        for A, B in H:
            for p in itertools.combinations(A, 2):
                for q in itertools.combinations(B, 2):
                    cov.add((p[0], p[1], q[0] + N, q[1] + N))
        bad = 0
        for S in itertools.combinations(range(2 * N), 5):
            As = [v for v in S if v < N]
            Bs = [v for v in S if v >= N]
            if len(As) >= 4 or len(Bs) >= 4:
                continue
            hit = False
            if len(As) == 3:
                for p in itertools.combinations(As, 2):
                    if (p[0], p[1], Bs[0], Bs[1]) in cov:
                        hit = True
                        break
            else:
                for q in itertools.combinations(Bs, 2):
                    if (As[0], As[1], q[0], q[1]) in cov:
                        hit = True
                        break
            bad += not hit
        print("(E) d=4: cross 4-sets covered by gauge-trivial 6-sets: %d (= |E_cross| 6720: %s) ; 5-sets with no "
              "block meeting them in >= 4 points: %d" % (len(cov), len(cov) == 6720, bad), flush=True)
