"""Independent checker -- extensions of rv1_lib (the own library of rv1_dual.py, used unchanged) for the
predicate "every p-subset spans >= lam edges" and for large moment matrices.  No producer module is imported.

Contents:
  * admissible(h, p, lam), flag_universe_pl(s, m, tau, p, lam)    (own enumeration, predicate-parametrised)
  * canonical forms of 7-vertex 4-graphs (vectorised over all 5040 vertex permutations)
  * accumulation of the integer moment matrices Y_tau = sum_H num_H * count_H,tau  (own counting, rv1_lib)
  * three exact PSD/PD criteria (each one alone is a proof; they are used to cross-check each other):
      - psd_exact      : rv1_lib symmetric elimination over Fractions (exact PSD test, any rank)        [n <= 64]
      - pd_bareiss     : fraction-free (Bareiss) elimination without pivoting; Sylvester's criterion:
                         A is positive definite iff every leading principal minor is > 0            [n <= NB]
      - pd_residual    : own float hint + exact residual.  D = diag(2^-e_i) with 4^e_i <= A_ii < 4^(e_i+1);
                         B ~ D A D (long double), G = Cholesky(B - delta I) in long double (own code),
                         L = round(2^P G) (integers), T = 2^(2P) D A D - L L^T computed EXACTLY (integers; L L^T by
                         int64 matrix products of balanced 21-bit limbs, no float).  If T has nonnegative diagonal
                         and is diagonally dominant, T is PSD (Gershgorin), so 2^(2P) D A D = L L^T + T is PSD and
                         A is PSD; with T_ii > sum_j!=i |T_ij| for all i (strict) T is PD and A is PD.
"""
from fractions import Fraction
from itertools import combinations, permutations
from math import comb

import numpy as np

import rv1_lib as L

N = 7


# ------------------------------------------------------------------------------------------------ predicate
def pset_counts(h, p):
    """Edge counts of the C(7, p) p-subsets of {0..6} (own lex mask h)."""
    out = []
    for P in combinations(range(N), p):
        sp = set(P)
        out.append(sum((h >> i) & 1 for i, e in enumerate(L.L7) if set(e) <= sp))
    return out


def admissible(h, p, lam):
    return all(c >= lam for c in pset_counts(h, p))


def flag_universe_pl(s, m, tau, p, lam):
    """Sorted canonical codes of all tau-flags on m vertices in which every p-subset of [m] spans >= lam edges
    (vacuous for p > m).  tau: own lex code on [s]."""
    root = [L.SUB[m].index(L.SUB[s][j]) for j in range(len(L.SUB[s])) if (tau >> j) & 1]
    nonroot = [j for j, e in enumerate(L.SUB[m]) if not set(e) <= set(range(s))]
    psets = [[j for j, e in enumerate(L.SUB[m]) if set(e) <= set(P)] for P in combinations(range(m), p)] if p <= m else []
    ident = tuple(range(s))
    U = tuple(range(s, m))
    out = set()
    for bits in range(1 << len(nonroot)):
        es = set(L.SUB[m][j] for j in root) | {L.SUB[m][nonroot[i]] for i in range(len(nonroot)) if (bits >> i) & 1}
        if not all(sum(1 for j in F if L.SUB[m][j] in es) >= lam for F in psets):
            continue
        hf = 0
        for e in es:
            hf |= 1 << L.IX[((e[0] * N + e[1]) * N + e[2]) * N + e[3]]
        assert L.code(hf, ident, s) == tau
        out.add(L.canon(hf, ident, U, m))
    return sorted(out)


# ------------------------------------------------------------------------------------------------ isomorphism
_PM = None


def canonical_forms(hs):
    """Canonical form (min over the 5040 vertex permutations of the relabelled own mask) of every graph in hs."""
    global _PM
    if _PM is None:
        _PM = np.array(L.perm_edge_maps(), dtype=np.int64)          # 5040 x 35
    out = []
    for h in hs:
        E = [i for i in range(35) if (h >> i) & 1]
        if not E:
            out.append(0)
            continue
        img = np.left_shift(np.int64(1), _PM[:, E]).sum(axis=1)    # distinct bits: sum == OR
        out.append(int(img.min()))
    return out


# ------------------------------------------------------------------------------------------------ moment matrices
class Acc:
    """Integer moment matrix of one labelled type, rows/cols = flags in order of first occurrence (object ints)."""

    def __init__(self):
        self.idx = {}
        self.cap = 16
        self.M = np.zeros((self.cap, self.cap), dtype=object)

    def _ix(self, f):
        i = self.idx.get(f)
        if i is None:
            i = len(self.idx)
            self.idx[f] = i
            if i >= self.cap:
                nc = 2 * self.cap
                M2 = np.zeros((nc, nc), dtype=object)
                M2[:self.cap, :self.cap] = self.M
                self.M, self.cap = M2, nc
        return i

    def add(self, d, w):
        M = self.M
        for (a, b), c in d.items():
            i = self._ix(a)
            j = self._ix(b)
            M = self.M
            M[i, j] += w * c

    def matrix(self):
        n = len(self.idx)
        flags = sorted(self.idx, key=self.idx.get)
        return flags, self.M[:n, :n]


def build_Y(support, blocks, progress=None):
    """support: list of (own mask, int weight). Returns {(s,m): {tau: Acc}} (all labelled types that occur)."""
    Y = {b: {} for b in blocks}
    for k, (h, w) in enumerate(support):
        for (s, m) in blocks:
            mc = L.moment_counts(h, s, m)
            for tau, d in mc.items():
                acc = Y[(s, m)].get(tau)
                if acc is None:
                    acc = Y[(s, m)][tau] = Acc()
                acc.add(d, w)
        if progress and (k + 1) % 200 == 0:
            progress(k + 1)
    return Y


# ------------------------------------------------------------------------------------------------ exact PSD tests
def pd_bareiss(A):
    """A: square numpy object array of Python ints, symmetric.  Fraction-free Gaussian elimination without
    pivoting (Bareiss); after step k the entry M[k,k] is the leading principal minor of order k+1.  Returns
    (is_PD, first index with minor <= 0 or None).  Sylvester: PD iff all leading principal minors > 0."""
    M = A.copy()
    n = M.shape[0]
    prev = 1
    for k in range(n):
        p = M[k, k]
        if p <= 0:
            return False, k
        if k + 1 < n:
            col = M[k + 1:, k].copy()
            row = M[k, k + 1:].copy()
            sub = M[k + 1:, k + 1:] * p - np.multiply.outer(col, row)
            # exact division (Bareiss); assert exactness on the diagonal as a self-check
            q = sub // prev
            if k < 3 or k % 50 == 0:
                assert all(int(x) == 0 for x in (sub - q * prev).ravel()[:200]), "Bareiss division not exact"
            M[k + 1:, k + 1:] = q
        prev = p
    return True, None


def _limbs(X, nl=3, b=21):
    """Balanced base-2^b digits of an int64 array X (|X| < 2^(b*nl-1)): X = sum_k D_k 2^(b k), |D_k| <= 2^(b-1)."""
    R = X.copy()
    out = []
    half = 1 << (b - 1)
    full = 1 << b
    for _ in range(nl):
        D = ((R + half) % full) - half
        out.append(D)
        R = (R - D) >> b
    assert not R.any(), "limb decomposition incomplete"
    return out


def exact_gram(Lint):
    """Exact L L^T of an int64 matrix with |L| < 2^62, as an object array of Python ints (int64 products only)."""
    n = Lint.shape[0]
    assert n <= 4096
    D = _limbs(Lint, 3, 21)
    for d in D:
        assert int(np.abs(d).max()) <= 1 << 20
    # |D_k D_l^T| entries <= n * 2^40 <= 2^52: exact in int64; group by shift s = k + l
    C = {}
    for k in range(3):
        for l in range(k, 3):
            P = D[k] @ D[l].T
            if l != k:
                P = P + P.T          # D_l D_k^T = (D_k D_l^T)^T
            C[k + l] = C.get(k + l, 0) + P
    G = np.zeros((n, n), dtype=object)
    for s, P in C.items():
        assert int(np.abs(P).max()) < 1 << 62
        G = G + P.astype(object) * (1 << (21 * s))
    return G


def chol_ld(B):
    """Lower Cholesky factor in long double (own right-looking implementation); None if a pivot is <= 0."""
    A = B.copy()
    n = A.shape[0]
    G = np.zeros_like(A)
    for k in range(n):
        d = A[k, k]
        if not d > 0:
            return None
        r = np.sqrt(d)
        G[k, k] = r
        if k + 1 < n:
            col = A[k + 1:, k] / r
            G[k + 1:, k] = col
            A[k + 1:, k + 1:] -= np.multiply.outer(col, col)
    return G


def pd_residual(A, P=61):
    """A: symmetric numpy object array of Python ints with positive diagonal.  Returns dict with PASS (bool),
    strict (bool), worst off-diagonal/diagonal ratio of T, delta, float lam_min of the scaled matrix."""
    n = A.shape[0]
    e = [(int(A[i, i]).bit_length() - 1) // 2 for i in range(n)]
    assert all(4 ** e[i] <= int(A[i, i]) < 4 ** (e[i] + 1) for i in range(n))
    assert 2 * P >= 2 * max(e), ("P too small", max(e))
    # long-double hint B ~ D A D  (|B_ij| < 4 for a PSD A; otherwise report failure)
    q = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(n):
            x = (int(A[i, j]) << 60) >> (e[i] + e[j])
            if abs(x) >= 1 << 63:
                return {"PASS": False, "why": "entry exceeds the PSD bound sqrt(A_ii A_jj)"}
            q[i, j] = x
    B = q.astype(np.longdouble) / np.longdouble(2.0 ** 60)
    lam = float(np.linalg.eigvalsh(B.astype(np.float64)).min())
    if not lam > 0:
        return {"PASS": False, "why": f"float lam_min {lam} <= 0", "lam_min": lam}
    two2P = 1 << (2 * P)
    K = np.zeros((n, n), dtype=object)
    for i in range(n):
        for j in range(n):
            K[i, j] = int(A[i, j]) << (2 * P - e[i] - e[j])
    last = None
    for frac in (0.5, 0.25, 0.1, 0.02):
        delta = lam * frac
        G = chol_ld(B - np.longdouble(delta) * np.eye(n, dtype=np.longdouble))
        if G is None:
            last = {"PASS": False, "why": f"long double Cholesky failed at delta {delta}", "lam_min": lam}
            continue
        Lint = np.rint(G * np.longdouble(2.0 ** P)).astype(np.int64)
        assert int(np.abs(Lint).max()) < 1 << 62
        T = K - exact_gram(Lint)
        diag = [int(T[i, i]) for i in range(n)]
        off = [sum(abs(int(T[i, j])) for j in range(n) if j != i) for i in range(n)]
        ok = all(diag[i] >= off[i] for i in range(n))
        strict = all(diag[i] > off[i] for i in range(n))
        worst = max(Fraction(off[i], diag[i]) if diag[i] > 0 else Fraction(10 ** 9) for i in range(n))
        res = {"PASS": ok, "strict": strict, "worst_ratio": float(worst), "delta": delta, "lam_min": lam,
               "P": P, "min_T_diag_over_2^2P": float(Fraction(min(diag), two2P))}
        if ok:
            return res
        last = res
    return last


# ------------------------------------------------------------------------------------------------ exact failure proofs
def scaled_float(A):
    """(e, B float64) with B = D A D, D = diag(2^-e_i), 4^e_i <= A_ii (A_ii > 0 required)."""
    n = A.shape[0]
    e = [(int(A[i, i]).bit_length() - 1) // 2 for i in range(n)]
    B = np.array([[float(Fraction(int(A[i, j]), 1 << (e[i] + e[j]))) for j in range(n)] for i in range(n)])
    return e, B


def neg_witness(A, bits=40):
    """Exact proof that A (object ints, positive diagonal) is NOT PSD: integer y with y^T A y < 0, built from the float
    eigenvector of the smallest eigenvalue of D A D.  Returns (found, float lam_min, exact y^T A y)."""
    n = A.shape[0]
    e, B = scaled_float(A)
    w, V = np.linalg.eigh(B)
    x = V[:, 0]
    E = max(e)
    y = [int(round(x[i] * 2 ** bits)) << (E - e[i]) for i in range(n)]
    q = 0
    for i in range(n):
        if y[i]:
            q += y[i] * sum(int(A[i, j]) * y[j] for j in range(n) if y[j])
    return q < 0, float(w[0]), q


def perturbed_controls(A):
    """Two exact controls for pd_residual on a PD matrix A: A_minus = A - (2 lam_min) * (scaled rank one along the
    smallest eigenvector), which is indefinite (proved exactly by neg_witness), and A_half = A - (lam_min / 2) * (same),
    which should stay PD.  Returns (A_minus, A_half, lam_min)."""
    n = A.shape[0]
    e, B = scaled_float(A)
    w, V = np.linalg.eigh(B)
    lam, x = float(w[0]), V[:, 0]
    out = []
    for mu in (2.0 * lam, 0.5 * lam):
        Ap = A.copy()
        for i in range(n):
            for j in range(n):
                Ap[i, j] = int(A[i, j]) - int(Fraction(mu * x[i] * x[j]) * (1 << (e[i] + e[j])))
        # keep it exactly symmetric
        for i in range(n):
            for j in range(i + 1, n):
                Ap[j, i] = Ap[i, j]
        out.append(Ap)
    return out[0], out[1], lam
