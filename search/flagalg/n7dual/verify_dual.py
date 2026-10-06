"""EXACT check of an N = 7 dual point for t(5, 4) (upper bound on every plain N = 7 flag-algebra certificate).

Self-contained: imports nothing from the flag-algebra engine (no ../flags.py, ../hypergraphs.py, ../n7/*); only
the 35-bit colex mask convention of the class list is shared (4-subsets of {0..6} ordered by their largest element,
then the next ...).  numpy is used ONLY to compute a floating-point Cholesky factor that serves as a hint; every
statement below is then verified in exact integer arithmetic.

Claim.  Input y = num_H / den (integers, num_H > 0, sum num_H = den) on 7-vertex 4-graphs H.  Checked:
 (i)   every H is admissible: every 5-subset of {0..6} contains an edge;
 (ii)  for every block (s, m) with m > s, 2m - s <= 7, m >= 4, i.e. (1,4) (2,4) (3,4) (3,5) (4,5) (5,6), and every
       type sigma (one labelled representative per isomorphism class: the minimum code over S_s), the moment matrix
           Y_sigma[a, b] = sum_H y_H * #{(theta, U1, U2): H[theta] = sigma exactly, flag(theta, U1) = a,
                                         flag(theta, U2) = b} / conf(s, m)
       (theta an injective s-tuple, (U1, U2) an ordered pair of disjoint (m-s)-sets of the other vertices, flags =
       classes of the labelled extension modulo permutations of the non-root vertices) is positive semidefinite.
       Flags/types that never occur give zero rows / zero matrices (PSD); blocks with m <= 3 are 1 x 1 with
       entry >= 0.  So Y_sigma is PSD on ALL admissible sigma-flags.
 (iii) V = sum_H y_H e(H) / 35 is computed exactly;
 (iv)  (reported) how many support graphs are ODD (every 5-set spans 1, 3 or 5 edges): if all are, V also bounds
       every plain N = 7 certificate of the odd-restricted problem (same weak-duality argument).
Weak duality: any plain certificate (b, PSD Q_sigma on any flag sets of these blocks) with
d(H) - sum <Q_sigma, M_sigma(H)> >= b for all admissible H gives b <= V - sum <Q_sigma, Y_sigma> <= V.

PSD test (exact).  For an integer symmetric matrix A (Y_sigma scaled by den * conf):
  * a zero diagonal entry requires the whole row to be zero (checked); such rows are removed;
  * n <= 40: exact LDL^T over Fractions;
  * otherwise: power-of-two diagonal scaling B = 2^-e A 2^-e (e_i = round(log2(A_ii)/2)); float Cholesky
    L of B - delta I (delta = half the float minimum eigenvalue of B); integer Lint = round(2^K L);
    C = Lint Lint^T EXACTLY (int64 limb products, bounds checked); T = B - C / 2^(2K) exactly (integers after a
    common power-of-two shift); PASS iff T is diagonally dominant with nonnegative diagonal (exact integer
    comparisons).  Then T >= 0, hence B = C/2^(2K) + T >= 0 (C is a Gram matrix), hence A >= 0.
Usage: python verify_dual.py DUAL.json
"""
import json
import math
import sys
import time
from fractions import Fraction
from itertools import combinations, permutations

import numpy as np

N, R, P = 7, 4, 5
E7 = sorted(combinations(range(N), R), key=lambda e: e[::-1])          # colex (the class-list convention)
BLOCKS = [(s, m) for s in range(N) for m in range(s + 1, N + 1) if 2 * m - s <= N and m >= R]


def decode(mask):
    return {E7[i] for i in range(len(E7)) if (mask >> i) & 1}


def admissible(es):
    return all(any(e in es for e in combinations(Pset, R)) for Pset in combinations(range(N), P))


class BlockTables:
    """Per block: root edges, free edges (as position tuples in range(m)), canonical maps."""

    def __init__(self, s, m):
        self.s, self.m = s, m
        Em = list(combinations(range(m), R))
        self.root_edges = [e for e in Em if max(e) < s]
        self.free_edges = [e for e in Em if max(e) >= s]
        ridx = {e: i for i, e in enumerate(self.root_edges)}
        fidx = {e: i for i, e in enumerate(self.free_edges)}
        # type canonical form: min code over S_s (labelled root graphs)
        nr = len(self.root_edges)
        rperm = [[ridx[tuple(sorted(g[v] for v in e))] for e in self.root_edges] for g in permutations(range(s))]
        self.type_canon = [min(sum(1 << pm[i] for i in range(nr) if (c >> i) & 1) for pm in rperm)
                           for c in range(1 << nr)]
        # flag canonical form: min code over permutations of the free vertices (roots fixed)
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
    nconf = {}
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
                nU = len(list(combinations(rest, m - s)))
                cnt_conf += nU * math.comb(N - m, m - s)
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
            nconf[(s, m)] = cnt_conf
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
    """Lint Lint^T exactly (Python ints) for an integer matrix Lint (list of lists of Python ints): split into
    26-bit limbs (top limb signed); every limb product < 2^52 and every int64 dot product of <= 1024 terms
    < 2^62 (no overflow)."""
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


def psd_certified(A, K=78):
    """Exact PSD test of an integer symmetric matrix A (list of lists).  Returns (ok, info)."""
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
    if n <= 40:
        return ldl_psd_exact(A), {"n": n, "method": "exact LDL^T (Fractions)"}
    e = [round(math.log2(A[i][i]) / 2) for i in range(n)]
    Bf = np.array([[math.ldexp(float(A[i][j]), -(e[i] + e[j])) for j in range(n)] for i in range(n)])
    lam = float(np.linalg.eigvalsh(Bf)[0])
    if lam <= 0:
        return False, {"n": n, "reason": f"float min eigenvalue of the scaled matrix {lam:.3e} <= 0"}
    delta = lam / 2
    L = np.linalg.cholesky(Bf - delta * np.eye(n))
    Lint = [[int(v) for v in row] for row in np.rint(np.ldexp(L, K)).tolist()]   # exact (integer-valued floats)
    C = gram_exact(Lint)
    # T = B - C / 2^(2K);  B_ij = A_ij / 2^(e_i + e_j).  Common shift S = 2K + max(0, max(e_i + e_j)).
    emax = max(0, 2 * max(e))
    S = 2 * K + emax
    def Tint(i, j):
        sh = S - e[i] - e[j]           # >= 2K >= 0
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
    t0 = time.time()
    d = json.load(open(sys.argv[1]))
    den = int(d["den"])
    support = []
    tot = 0
    Vnum = 0
    for row in d["support"]:
        num, mask = int(row["num"]), int(row["mask"])
        assert num > 0 and 0 <= mask < (1 << 35)
        es = decode(mask)
        assert admissible(es), f"support graph {mask} not admissible"
        support.append((num, es))
        tot += num
        Vnum += num * len(es)
    assert tot == den, "sum y != 1"
    n_odd = sum(1 for _, es in support
                if all(sum(e in es for e in combinations(S5, R)) % 2 == 1 for S5 in combinations(range(N), P)))
    print(f"# support graphs that are ODD (every 5-set spans 1, 3 or 5 edges): {n_odd} of {len(support)}", flush=True)
    assert len({int(r["mask"]) for r in d["support"]}) == len(d["support"]), "duplicate masks"
    V = Fraction(Vnum, 35 * den)
    print(f"# support {len(support)} admissible graphs, sum y = 1, V = sum y d = {float(V):.15f}", flush=True)
    if "value" in d:
        assert Fraction(d["value"]) == V, "claimed value differs"
        print("# claimed value equal: True", flush=True)
    tabs = [BlockTables(s, m) for s, m in BLOCKS]
    print(f"# blocks {BLOCKS}; tables built {time.time() - t0:.1f}s", flush=True)
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
    out = {"PASS": ok_all, "support": len(support), "support_odd": n_odd, "value": str(V), "value_float": float(V),
           "value_ceil12": math.ceil(V * 10 ** 12) / 10 ** 12, "blocks": res, "time_s": time.time() - t0}
    print(json.dumps({k: v for k, v in out.items() if k != "blocks"}), flush=True)
    if len(sys.argv) > 2:
        json.dump(out, open(sys.argv[2], "w"), indent=1)
    if ok_all:
        print(f"PASS: no plain N = 7 flag-algebra certificate proves t(5,4) > {V} (~ {float(V):.12f})")
    else:
        print("FAIL")
        sys.exit(1)


if __name__ == "__main__":
    main()
