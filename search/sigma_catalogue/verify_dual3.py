#!/usr/bin/env python3
"""EXACT check of an N = 7 dual point for sigma(J_4) (3-graphs; r = 3 port of the check of Table n7limits).
Producer code (not reviewed).  Self-contained: imports nothing from the flag-algebra engine; only the 35-bit colex
mask convention (3-subsets of {0..6} ordered by their largest element, then the next, ...) is shared.  numpy is used
ONLY for a floating-point Cholesky hint and for int64 limb products; every statement is verified in exact integers.

Claim.  Input y = num_H / den (integers, num_H > 0, sum num_H = den) on 7-vertex 3-graphs H.  Checked:
 (i)   every H is J_4-free (no vertex whose link contains K_4; J_4 = 123 124 125 134 135 145);
 (ii)  for every block (s, m) with m > s, 2m - s <= 7, m >= 3 -- (0,3) (1,3) (2,3) (1,4) (2,4) (3,4) (3,5) (4,5)
       (5,6), all blocks of the plain N = 7 method for 3-graphs -- and every type sigma (one labelled representative
       per isomorphism class: the minimum code over S_s), the moment matrix
           Y_sigma[a, b] = sum_H y_H * #{(theta, U1, U2): H[theta] = sigma exactly, flag(theta, U1) = a,
                                         flag(theta, U2) = b} / conf(s, m)
       (theta an injective s-tuple, (U1, U2) an ordered pair of disjoint (m-s)-sets of the other vertices, flags =
       classes of the labelled extension modulo permutations of the non-root vertices) is positive semidefinite.
       Flags / types that never occur give zero rows / zero matrices.  So Y_sigma is PSD on ALL sigma-flags.
 (iii) V = sum_H y_H f(H), f(H) = F(H)/210, F(H) = 210 - sum over 4-sets Q of C(e_H(Q), 2) (cosig3), exactly.
Weak duality: any plain certificate (b, PSD Q_sigma on any flag sets of these blocks; other labellings of a type only
permute rows and columns; 1 x 1 blocks (m <= 2 or m = s) give terms -q d_tau(H) with q >= 0; mixed flag sizes reduce to
the top block by the lifting T M_top T^T) with f(H) - sum <Q_sigma, M_sigma(H)> >= b for all J_4-free H satisfies
  b <= sum_H y_H (f(H) - sum <Q, M(H)>) = V - sum <Q_sigma, Y_sigma> <= V.
Hence: no plain N = 7 certificate proves sigma(J_4) < 1 - V.  If V < 41/57 the N = 7 optimum is < 41/57, i.e. the
plain N = 7 method cannot prove sigma(J_4) = 16/57 (it cannot even reach any bound below 1 - V).
PSD test (exact), as verify_dual2: zero diagonal => zero row (checked), rows dropped; n <= 40: exact LDL^T; else
B = 2^-e A 2^-e, float Cholesky L of B - delta I (delta = half the float min eigenvalue), Lint = round(2^K L),
C = Lint Lint^T exactly (26-bit limbs, int64), T = B - C / 2^(2K) exactly; PASS iff T is diagonally dominant with
nonnegative diagonal (exact integers) => T PSD => B = C/2^(2K) + T PSD => A PSD.
Usage: python verify_dual3.py DUAL.json [OUT.json] [--ldl-all]
"""
import json
import math
import sys
import time
from fractions import Fraction
from itertools import combinations, permutations

import numpy as np

N, R = 7, 3
E7 = sorted(combinations(range(N), R), key=lambda e: e[::-1])
ALL_BLOCKS = [(s, m) for s in range(N) for m in range(s + 1, N + 1) if 2 * m - s <= N and m >= R]
assert ALL_BLOCKS == [(0, 3), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4), (3, 5), (4, 5), (5, 6)]


def decode(mask):
    return {E7[i] for i in range(len(E7)) if (mask >> i) & 1}


def j4_free(es):
    for x in range(N):
        oth = [v for v in range(N) if v != x]
        link = {(a, b) for a, b in combinations(oth, 2) if tuple(sorted((x, a, b))) in es}
        for Q in combinations(oth, 4):
            if all(p in link for p in combinations(Q, 2)):
                return False
    return True


def F_num(es):
    s = 0
    for Q in combinations(range(N), 4):
        c = sum(t in es for t in combinations(Q, 3))
        s += c * (c - 1) // 2
    return 210 - s


class BlockTables:
    def __init__(self, s, m):
        self.s, self.m = s, m
        Em = list(combinations(range(m), R))
        self.root_edges = [e for e in Em if max(e) < s]
        self.free_edges = [e for e in Em if max(e) >= s]
        ridx = {e: i for i, e in enumerate(self.root_edges)}
        fidx = {e: i for i, e in enumerate(self.free_edges)}
        nr = len(self.root_edges)
        rperm = [[ridx[tuple(sorted(g[v] for v in e))] for e in self.root_edges] for g in permutations(range(s))]
        self.type_canon = [min(sum(1 << pm[i] for i in range(nr) if (c >> i) & 1) for pm in rperm)
                           for c in range(1 << nr)]
        nf = len(self.free_edges)
        fperm = []
        for pu in permutations(range(s, m)):
            g = list(range(s)) + list(pu)
            fperm.append([fidx[tuple(sorted(g[v] for v in e))] for e in self.free_edges])
        self.flag_canon = [min(sum(1 << pm[i] for i in range(nf) if (c >> i) & 1) for pm in fperm)
                           for c in range(1 << nf)]
        self.conf = math.factorial(N) // math.factorial(N - s) * math.comb(N - s, m - s) * math.comb(N - m, m - s)


def moment_counts(support, tabs):
    """Y[(s, m, sigma)][(a, b)] = sum_H num_H * count (integers)."""
    Y = {}
    for num, es in support:
        for tb in tabs:
            s, m = tb.s, tb.m
            cnt_conf = 0
            for th in permutations(range(N), s):
                rc = 0
                for i, e in enumerate(tb.root_edges):
                    if tuple(sorted(th[v] for v in e)) in es:
                        rc |= 1 << i
                rest = [v for v in range(N) if v not in th]
                cnt_conf += math.comb(len(rest), m - s) * math.comb(N - m, m - s)
                if tb.type_canon[rc] != rc:
                    continue
                fc = {}
                for U in combinations(rest, m - s):
                    order = th + U
                    c = 0
                    for i, e in enumerate(tb.free_edges):
                        if tuple(sorted(order[v] for v in e)) in es:
                            c |= 1 << i
                    fc[U] = tb.flag_canon[c]
                M = Y.setdefault((s, m, rc), {})
                for U1 in combinations(rest, m - s):
                    r2 = [v for v in rest if v not in U1]
                    for U2 in combinations(r2, m - s):
                        key = (fc[U1], fc[U2])
                        M[key] = M.get(key, 0) + num
            assert cnt_conf == tb.conf
    return Y


def to_matrix(M):
    flags = sorted({a for a, _ in M} | {b for _, b in M})
    ix = {f: i for i, f in enumerate(flags)}
    n = len(flags)
    A = [[0] * n for _ in range(n)]
    for (a, b), v in M.items():
        A[ix[a]][ix[b]] += v
    return A, flags


def ldl_psd_exact(A):
    n = len(A)
    A = [[Fraction(x) for x in row] for row in A]
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
                for j in range(k + 1, n):
                    A[i][j] -= f * A[k][j]
    return True


def gram_exact(Lint):
    r = len(Lint)
    n = len(Lint[0])
    assert n <= 1024
    big = max(abs(v) for row in Lint for v in row)
    nl = max(1, -(-(big.bit_length() + 1) // 26))
    limbs = []
    rest = [row[:] for row in Lint]
    mask = (1 << 26) - 1
    for _ in range(nl - 1):
        limbs.append(np.array([[v & mask for v in row] for row in rest], dtype=np.int64))
        rest = [[v >> 26 for v in row] for row in rest]
    top = np.array(rest, dtype=np.int64)
    assert int(np.abs(top).max()) < (1 << 26)
    limbs.append(top)
    acc = [[0] * r for _ in range(r)]
    for p in range(nl):
        for q in range(nl):
            Gl = (limbs[p] @ limbs[q].T).tolist()   # int64, |entries| < 1024 * 2^52 = 2^62
            sh = 26 * (p + q)
            for i in range(r):
                ai, gi = acc[i], Gl[i]
                for j in range(r):
                    ai[j] += gi[j] << sh
    return acc


LDL_ALL = False          # --ldl-all: exact LDL^T for every size (slow; used for the singular positive control phi)


def psd_certified(A, K=78):
    n = len(A)
    for i in range(n):
        for j in range(i):
            if A[i][j] != A[j][i]:
                return False, {"reason": "not symmetric"}
    keep = []
    for i in range(n):
        if A[i][i] < 0:
            return False, {"reason": "negative diagonal"}
        if A[i][i] == 0:
            if any(A[i][j] != 0 for j in range(n)):
                return False, {"reason": "zero diagonal with nonzero row"}
        else:
            keep.append(i)
    A = [[A[i][j] for j in keep] for i in keep]
    n = len(A)
    if n == 0:
        return True, {"n": 0}
    if n <= 40 or LDL_ALL:
        return ldl_psd_exact(A), {"n": n, "method": "exact LDL^T (Fractions)"}
    e = [round(math.log2(A[i][i]) / 2) for i in range(n)]
    Bf = np.array([[math.ldexp(float(A[i][j]), -(e[i] + e[j])) for j in range(n)] for i in range(n)])
    lam = float(np.linalg.eigvalsh(Bf)[0])
    if lam <= 0:
        return False, {"n": n, "reason": f"float min eigenvalue of the scaled matrix {lam:.3e} <= 0"}
    delta = lam / 2
    L = np.linalg.cholesky(Bf - delta * np.eye(n))
    Lint = [[int(v) for v in row] for row in np.rint(np.ldexp(L, K)).tolist()]
    C = gram_exact(Lint)
    emax = max(0, 2 * max(e))
    S = 2 * K + emax

    def Tint(i, j):
        sh = S - e[i] - e[j]
        return (A[i][j] << sh) - (C[i][j] << (S - 2 * K))
    worst = None
    for i in range(n):
        row = [Tint(i, j) for j in range(n)]
        di = row[i]
        off = sum(abs(row[j]) for j in range(n) if j != i)
        if di < off:
            return False, {"n": n, "method": "cholesky+residual", "reason": f"row {i} not dominant", "K": K,
                           "lam_float": lam}
        r = Fraction(off, di) if di else Fraction(0)
        worst = r if worst is None or r > worst else worst
    return True, {"n": n, "method": "cholesky+exact residual, diagonal dominance", "K": K, "lam_float": lam,
                  "max_offdiag_over_diag": float(worst)}


def main():
    global LDL_ALL
    t0 = time.time()
    if "--ldl-all" in sys.argv:
        LDL_ALL = True
        sys.argv.remove("--ldl-all")
    d = json.load(open(sys.argv[1]))
    den = int(d["den"])
    support, tot, Vnum = [], 0, 0
    for row in d["support"]:
        num, mask = int(row["num"]), int(row["mask"])
        assert num > 0 and 0 <= mask < (1 << 35)
        es = decode(mask)
        assert j4_free(es), f"support graph {mask} contains J_4"
        support.append((num, es))
        tot += num
        Vnum += num * F_num(es)
    assert tot == den, "sum y != 1"
    assert len({int(r["mask"]) for r in d["support"]}) == len(d["support"]), "duplicate masks"
    V = Fraction(Vnum, 210 * den)
    print(f"# support {len(support)} J_4-free graphs, sum y = 1, V = sum y f = {float(V):.15f}; "
          f"V - 41/57 = {float(V - Fraction(41, 57)):.6e}", flush=True)
    if "value" in d:
        assert Fraction(d["value"]) == V, "claimed value differs"
    tabs = [BlockTables(s, m) for s, m in ALL_BLOCKS]
    Y = moment_counts(support, tabs)
    print(f"# moment matrices for {len(Y)} (block, type) pairs, {time.time() - t0:.1f}s", flush=True)
    ok_all = True
    res = []
    for key in sorted(Y):
        A, flags = to_matrix(Y[key])
        ok, info = psd_certified(A)
        info.update({"block": key[:2], "type_code": key[2], "flags_occurring": len(flags), "PSD": ok})
        res.append(info)
        print(json.dumps(info), flush=True)
        ok_all &= ok
    t = Fraction(41, 57)
    out = {"PASS": ok_all, "problem": "J4", "blocks": ALL_BLOCKS, "support": len(support), "value": str(V),
           "value_float": float(V), "value_ceil12": str(Fraction(math.ceil(V * 10 ** 12), 10 ** 12)),
           "V_below_41_57": V < t, "41_57_minus_V": str(t - V), "41_57_minus_V_float": float(t - V),
           "sigma_lower_end": str(1 - V), "block_checks": res, "time_s": time.time() - t0}
    print(json.dumps({k: v for k, v in out.items() if k != "block_checks"}), flush=True)
    if len(sys.argv) > 2:
        json.dump(out, open(sys.argv[2], "w"), indent=1)
    if ok_all:
        tail = (f"V < 41/57 by {float(t - V):.6e}, so the plain N = 7 optimum is strictly below 41/57."
                if V < t else f"V >= 41/57 (V - 41/57 = {float(V - t):.6e}): nothing below 41/57 is shown.")
        print(f"PASS: no plain N = 7 flag-algebra certificate (blocks {ALL_BLOCKS}) proves sigma(J_4) < 1 - V = "
              f"{float(1 - V):.13f} (V = {float(V):.13f}); " + tail)
    else:
        print("FAIL")
        sys.exit(1)


if __name__ == "__main__":
    main()
