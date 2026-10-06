"""Exact rational certificates for flag-algebra SDP bounds.

Two rounding methods, both ending with an exact recomputation of
    b_exact = min_H ( obj(H) - sum_sigma <Q_sigma, M_sigma(H)> )
in rational arithmetic, where every M_sigma(H) entry is the exact rational count/denominator and every Q_sigma is a
rational matrix whose positive semidefiniteness is verified exactly (symmetric Gaussian elimination over Q).

1. 'eig' (always applicable, never sharp): Q = V diag(w) V^T (float); drop w_i <= drop_tol * max(w); round w_i >= 0 and
   v_i to dyadic rationals with `bits` bits; Q_exact = sum_i w_i~ v_i~ v_i~^T, PSD by construction (sum of
   rank-one PSD terms with nonnegative weights); b_exact is a rigorous bound slightly below the float bound.
2. 'sharp' (Flagmatic-style, for problems with a nice exact answer b*): rationalise the kernel of each Q_sigma
   (rational reconstruction of its RREF), write Q_sigma = R^T X R with rows of R spanning the orthogonal complement
   of the kernel, and solve the linear equalities  obj(H) - <Q, M(H)> = b*  for the sharp graphs exactly, as the
   least-norm correction of the rounded float X.  Success iff every X_sigma is positive definite (exact check) and
   the exact minimum equals b*.
"""
from __future__ import annotations

from fractions import Fraction
from math import lcm

import numpy as np


# ----------------------------------------------------------------------------------------------- exact PSD test
def is_psd_exact(Q) -> tuple:
    """Exact PSD test of a symmetric rational matrix (list of lists of Fraction/int).

    Symmetric Gaussian elimination: pick a positive diagonal pivot, take the Schur complement; a zero diagonal entry
    requires a zero row; a negative diagonal entry means not PSD.  Returns (is_psd, rank, min_pivot)."""
    n = len(Q)
    A = [[Fraction(Q[i][j]) for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if A[i][j] != A[j][i]:
                return False, None, None
    alive = list(range(n))
    rank = 0
    minpiv = None
    while alive:
        piv = None
        for i in alive:
            if A[i][i] < 0:
                return False, rank, A[i][i]
            if A[i][i] > 0 and (piv is None or A[i][i] > A[piv][piv]):
                piv = i
        if piv is None:  # all remaining diagonals are zero: the remaining block must vanish
            for i in alive:
                for j in alive:
                    if A[i][j] != 0:
                        return False, rank, Fraction(0)
            break
        p = A[piv][piv]
        minpiv = p if minpiv is None or p < minpiv else minpiv
        rank += 1
        alive.remove(piv)
        rowp = A[piv]
        for i in alive:
            f = A[i][piv] / p
            if f:
                Ai = A[i]
                for j in alive:
                    if rowp[j]:
                        Ai[j] -= f * rowp[j]
    return True, rank, minpiv


# ----------------------------------------------------------------------------------------------- method 1: eig
def round_psd_eig(Q: np.ndarray, bits: int = 40, drop_tol: float = 1e-12):
    """Return (Qnum, shift) with Q_exact = Qnum / 2**shift, Qnum an object array of Python ints; exactly PSD."""
    w, V = np.linalg.eigh(Q)
    wmax = max(float(w.max()), 0.0) if w.size else 0.0
    keep = [i for i in range(w.size) if w[i] > drop_tol * max(wmax, 1e-300)]
    scale = 1 << bits
    n = Q.shape[0]
    Qnum = np.zeros((n, n), dtype=object)
    Qnum[:, :] = 0
    for i in keep:
        lam = int(round(float(w[i]) * scale))
        if lam <= 0:
            continue
        v = [int(round(float(x) * scale)) for x in V[:, i]]
        vv = np.array(v, dtype=object)
        Qnum = Qnum + lam * np.outer(vv, vv)
    return Qnum, 3 * bits


# ----------------------------------------------------------------------------------------------- exact evaluation
def block_inner_exact(blk, Qnum, nH: int) -> list:
    """S[h] = sum_{a,b} Qnum[a,b] * cnt_h(a,b)  (Python ints)."""
    S = [0] * nH
    Ql = Qnum.tolist()
    for h, a, b, c in zip(blk.h.tolist(), blk.a.tolist(), blk.b.tolist(), blk.cnt.tolist()):
        S[h] += Ql[a][b] * c
    return S


def exact_slacks(prob, obj_exact: list, Qexact: list) -> list:
    """slack(H) = obj(H) - sum_sigma <Q_sigma, M_sigma(H)>, exactly.  Qexact[i] = (Qnum, den) with Q = Qnum/den."""
    nH = prob.nH
    sl = [Fraction(x) for x in obj_exact]
    for blk, (Qnum, den) in zip(prob.blocks, Qexact):
        S = block_inner_exact(blk, Qnum, nH)
        D = den * blk.denom
        for h in range(nH):
            if S[h]:
                sl[h] -= Fraction(S[h], D)
    return sl


def certify_eig(prob, obj_exact: list, Qs: list, bits: int = 40, drop_tol: float = 1e-12, check_psd: bool = True):
    Qexact, psd = [], []
    for Q in Qs:
        Qnum, sh = round_psd_eig(Q, bits, drop_tol)
        Qexact.append((Qnum, 1 << sh))
        if check_psd:
            ok, rank, _ = is_psd_exact(Qnum.tolist())
            psd.append((ok, rank))
    sl = exact_slacks(prob, obj_exact, Qexact)
    i = min(range(len(sl)), key=lambda k: sl[k])
    return {"b_exact": sl[i], "argmin": i, "psd": psd, "Qexact": Qexact, "slacks": sl}


# ----------------------------------------------------------------------------------------------- exact linear algebra
def rref(rows: list, ncols: int):
    """Exact RREF of a list of rows (Fractions). Returns (R, pivots)."""
    A = [[Fraction(x) for x in row] for row in rows]
    piv = []
    r = 0
    for c in range(ncols):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = 1 / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
        if r == len(A):
            break
    return A[:r], piv


def nullspace_rows(K: list, n: int) -> list:
    """Rational basis (as rows, integer-scaled) of {x : K x = 0}."""
    if not K:
        return [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    Rr, piv = rref(K, n)
    free = [j for j in range(n) if j not in piv]
    out = []
    for f in free:
        x = [Fraction(0)] * n
        x[f] = Fraction(1)
        for row, p in zip(Rr, piv):
            x[p] = -row[f]
        L = lcm(*[v.denominator for v in x])
        out.append([int(v * L) for v in x])
    return out


def float_rref(Kf: np.ndarray, tol: float = 1e-9):
    """Row-reduced form of the row space of Kf (k x n, full row rank), with pivot columns chosen by QR with column
    pivoting (numerically stable; entries stay O(1)).  Returns (rows, pivots)."""
    import scipy.linalg
    k = Kf.shape[0]
    _, _, P = scipy.linalg.qr(Kf, pivoting=True, mode="economic")
    piv = sorted(int(x) for x in P[:k])
    R = np.linalg.solve(Kf[:, piv], Kf)
    return R, piv


def solve_least_norm(Arows: list, rhs: list, u0: list):
    """Exact u = u0 + A^T w minimising ||u - u0|| subject to A u = rhs (A integer rows, rhs Fractions).
    Returns None if inconsistent."""
    m = len(Arows)
    if m == 0:
        return list(u0)
    res = [Fraction(rhs[i]) - sum(Fraction(a) * u for a, u in zip(Arows[i], u0) if a) for i in range(m)]
    Aobj = np.array(Arows, dtype=object)
    G = (Aobj @ Aobj.T).tolist()                               # Gram matrix, exact ints
    aug = [[Fraction(x) for x in G[i]] + [res[i]] for i in range(m)]
    Rr, piv = rref(aug, m + 1)
    if m in piv:
        return None
    w = [Fraction(0)] * m
    for row, p in zip(Rr, piv):
        w[p] = row[m]
    u = list(u0)
    for i in range(m):
        if w[i]:
            for j, a in enumerate(Arows[i]):
                if a:
                    u[j] += a * w[i]
    # verify all equations (including dependent ones)
    for i in range(m):
        if sum(Fraction(a) * x for a, x in zip(Arows[i], u) if a) != rhs[i]:
            return None
    return u


# ----------------------------------------------------------------------------------------------- method 2: sharp
def certify_sharp(prob, obj_exact: list, Qs: list, slacks_float: np.ndarray, target: Fraction,
                  ker_tol: float = 1e-6, sharp_tol: float = 1e-6, maxden: int = 1000, xbits: int = 30,
                  verbose: bool = False):
    """Flagmatic-style sharp rounding.  Returns dict with success flag and b_exact (= target on success)."""
    info = {"target": target}
    Rs, Xs0 = [], []
    for bi, Q in enumerate(Qs):
        n = Q.shape[0]
        w, V = np.linalg.eigh(Q)
        kern = V[:, w < ker_tol]
        if kern.shape[1]:
            Kr, _ = float_rref(kern.T)
            Krat = [[Fraction(float(x)).limit_denominator(maxden) for x in row] for row in Kr]
            err = float(np.abs(Q @ np.array([[float(x) for x in row] for row in Krat]).T).max())
            if err > 1e-4:
                info.setdefault("warnings", []).append(f"block {bi}: kernel rationalisation residual {err:.2e}")
        else:
            Krat = []
        R = nullspace_rows(Krat, n)                              # (n-k) x n integer rows
        Rf = np.array(R, dtype=float).reshape(-1, n)
        if Rf.shape[0]:
            G = Rf @ Rf.T
            Gi = np.linalg.inv(G)
            Xf = Gi @ Rf @ Q @ Rf.T @ Gi
        else:
            Xf = np.zeros((0, 0))
        Rs.append(R)
        Xs0.append(Xf)
        if verbose:
            print(f"  block {bi}: n={n}, kernel dim {len(Krat)}, X size {len(R)}", flush=True)
    sharp = [h for h in range(prob.nH) if slacks_float[h] < sharp_tol]
    info["n_sharp"] = len(sharp)
    # unknowns: upper-triangle entries of every X_sigma
    unk = []
    for bi, R in enumerate(Rs):
        k = len(R)
        for j in range(k):
            for i in range(j + 1):
                unk.append((bi, i, j))
    col = {u: c for c, u in enumerate(unk)}
    L = lcm(*[blk.denom for blk in prob.blocks]) if prob.blocks else 1
    Arows = [[0] * len(unk) for _ in sharp]
    for bi, (blk, R) in enumerate(zip(prob.blocks, Rs)):
        k = len(R)
        if k == 0:
            continue
        Robj = np.array(R, dtype=object)
        fac = L // blk.denom
        sel = np.isin(blk.h, sharp)
        hs, aa, bb, cc = blk.h[sel], blk.a[sel], blk.b[sel], blk.cnt[sel]
        for si, h in enumerate(sharp):
            msk = hs == h
            if not msk.any():
                continue
            n = blk.size
            M = np.zeros((n, n), dtype=object)
            M[:, :] = 0
            for a, b, c in zip(aa[msk].tolist(), bb[msk].tolist(), cc[msk].tolist()):
                M[a, b] += c
            Nm = (Robj @ M @ Robj.T).tolist()                      # <R^T X R, M> = <X, R M R^T>
            row = Arows[si]
            for j in range(k):
                for i in range(j + 1):
                    v = Nm[i][j] if i == j else 2 * Nm[i][j]
                    if v:
                        row[col[(bi, i, j)]] += v * fac
    rhs = [(Fraction(obj_exact[h]) - target) * L for h in sharp]
    scale = 1 << xbits
    u0 = [Fraction(int(round(Xs0[bi][i, j] * scale)), scale) for (bi, i, j) in unk]
    u = solve_least_norm(Arows, rhs, u0)
    if u is None:
        info["success"] = False
        info["reason"] = "sharp equalities inconsistent"
        return info
    Qexact, psd_ok = [], True
    for bi, R in enumerate(Rs):
        k = len(R)
        n = prob.blocks[bi].size
        X = [[Fraction(0)] * k for _ in range(k)]
        for c, (bj, i, j) in enumerate(unk):
            if bj == bi:
                X[i][j] = X[j][i] = u[c]
        ok, rank, minpiv = is_psd_exact(X) if k else (True, 0, None)
        info.setdefault("X_rank", []).append((bi, rank, k))
        if not ok:  # PSD is all that validity needs; rank < k only means a larger kernel than rationalised
            psd_ok = False
            info.setdefault("psd_fail", []).append((bi, ok, rank, k))
        # Q = R^T X R as integer matrix / den
        den = lcm(*[x.denominator for row in X for x in row]) if k else 1
        Xi = np.array([[int(x * den) for x in row] for row in X], dtype=object).reshape(k, k)
        Robj = np.array(R, dtype=object).reshape(k, n)
        Qnum = (Robj.T @ Xi @ Robj) if k else np.zeros((n, n), dtype=object)
        if not k:
            Qnum[:, :] = 0
        Qexact.append((Qnum, den))
    sl = exact_slacks(prob, obj_exact, Qexact)
    i = min(range(len(sl)), key=lambda kk: sl[kk])
    info.update({"success": bool(psd_ok and sl[i] == target), "psd_ok": psd_ok, "b_exact": sl[i], "argmin": i,
                 "Qexact": Qexact, "slacks": sl})
    return info
