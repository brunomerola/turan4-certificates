"""The flag-algebra SDP, solved with Clarabel (interior point, PSD triangle cone).

Primal (variables b, optional LP variables c >= 0, Q_sigma PSD):
    maximise b
    s.t. for every admissible H:  b - sum_j R[H, j] c_j + sum_sigma <Q_sigma, M_sigma(H)> <= d(H)
         G c <= h,  c >= 0,  Q_sigma PSD.
For the Turan problem R is empty and d(H) is the edge density; for the lottery packing problem d = 0 and
R[H, j] = rho_j(H) (see lottery.py).

Clarabel 0.11 conventions (checked on a 3x3 test, see test_flagalg.py): PSDTriangleConeT(n) holds
svec(X) = upper triangle in column-major order (X11, X12, X22, X13, X23, X33, ...), off-diagonals scaled by sqrt(2);
so <Q, M> = svec(Q) . svec(M).  Constraint form: A x + s = b, s in cone.
"""
from __future__ import annotations

import time

import clarabel
import numpy as np
import scipy.sparse as sp

SQ2 = np.sqrt(2.0)


def svec_index(i, j):
    """Index of (i, j) in the column-major upper-triangle vectorisation (i, j may be arrays)."""
    lo = np.minimum(i, j)
    hi = np.maximum(i, j)
    return hi * (hi + 1) // 2 + lo


def unsvec(x: np.ndarray, n: int) -> np.ndarray:
    Q = np.zeros((n, n))
    k = 0
    for j in range(n):
        for i in range(j + 1):
            v = x[k] if i == j else x[k] / SQ2
            Q[i, j] = Q[j, i] = v
            k += 1
    return Q


def solve_flag_sdp(prob, d: np.ndarray, R=None, G=None, h=None, extra_psd=(), tol: float = 1e-10,
                   max_iter: int = 400, verbose: bool = False, settings: dict | None = None, psd_maps=()):
    """extra_psd: sizes of additional PSD blocks W_j that appear only in the rows G (columns of G are
    [c (nc) | svec(W_1) | svec(W_2) ...]).
    psd_maps: alternatively, ((n, S),) with W = unsvec(S @ w), w free variables occupying the G columns after c
    (used for symmetry-reduced W); then extra_psd must be empty."""
    nH = prob.nH
    nc = 0 if R is None else R.shape[1]
    if psd_maps:
        assert not extra_psd
        wtri = [S.shape[1] for (n, S) in psd_maps]   # number of free variables per map
    else:
        wtri = [n * (n + 1) // 2 for n in extra_psd]
    nw = int(sum(wtri))
    sizes = [blk.size for blk in prob.blocks]
    tri = [n * (n + 1) // 2 for n in sizes]
    off = np.cumsum([1 + nc + nw] + tri)     # start column of each Q_sigma block
    nx = int(off[-1])

    rows, cols, vals = [np.arange(nH)], [np.zeros(nH, dtype=np.int64)], [np.ones(nH)]
    if nc:
        Rc = sp.coo_matrix(R)
        rows.append(Rc.row); cols.append(Rc.col + 1); vals.append(-Rc.data)
    for bi, blk in enumerate(prob.blocks):
        coef = blk.cnt / blk.denom
        coef = np.where(blk.a == blk.b, coef, coef / SQ2)
        rows.append(blk.h); cols.append(off[bi] + svec_index(blk.a, blk.b)); vals.append(coef)
    A1 = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nH, nx))
    blocks_A = [A1]
    bvec = [np.asarray(d, dtype=float)]
    cones = [clarabel.NonnegativeConeT(nH)]
    if nc:
        if G is not None and G.shape[0]:
            Gc = sp.coo_matrix(G)
            blocks_A.append(sp.coo_matrix((Gc.data, (Gc.row, Gc.col + 1)), shape=(G.shape[0], nx)))
            bvec.append(np.asarray(h, dtype=float))
            cones.append(clarabel.NonnegativeConeT(G.shape[0]))
        blocks_A.append(sp.coo_matrix((-np.ones(nc), (np.arange(nc), np.arange(nc) + 1)), shape=(nc, nx)))
        bvec.append(np.zeros(nc))
        cones.append(clarabel.NonnegativeConeT(nc))
    woff = 1 + nc
    for n, t in zip(extra_psd, wtri):
        blocks_A.append(sp.coo_matrix((-np.ones(t), (np.arange(t), woff + np.arange(t))), shape=(t, nx)))
        bvec.append(np.zeros(t))
        cones.append(clarabel.PSDTriangleConeT(n))
        woff += t
    for (n, S), t in zip(psd_maps, wtri):
        Sc = sp.coo_matrix(S)
        blocks_A.append(sp.coo_matrix((-Sc.data, (Sc.row, Sc.col + woff)), shape=(n * (n + 1) // 2, nx)))
        bvec.append(np.zeros(n * (n + 1) // 2))
        cones.append(clarabel.PSDTriangleConeT(n))
        woff += t
    for bi, n in enumerate(sizes):
        t = tri[bi]
        blocks_A.append(sp.coo_matrix((-np.ones(t), (np.arange(t), off[bi] + np.arange(t))), shape=(t, nx)))
        bvec.append(np.zeros(t))
        cones.append(clarabel.PSDTriangleConeT(n))
    A = sp.vstack(blocks_A).tocsc()
    b = np.concatenate(bvec)
    q = np.zeros(nx)
    q[0] = -1.0
    P = sp.csc_matrix((nx, nx))

    st = clarabel.DefaultSettings()
    st.verbose = verbose
    st.tol_gap_abs = tol
    st.tol_gap_rel = tol
    st.tol_feas = tol
    st.tol_ktratio = 1e-8
    st.max_iter = max_iter
    for key, val in (settings or {}).items():
        setattr(st, key, val)
    t0 = time.time()
    solver = clarabel.DefaultSolver(P, q, A, b, cones, st)
    sol = solver.solve()
    el = time.time() - t0
    x = np.array(sol.x)
    z = np.array(sol.z)
    Qs = [unsvec(x[off[bi]:off[bi + 1]], n) for bi, n in enumerate(sizes)]
    return {
        "status": str(sol.status),
        "b": float(x[0]),
        "c": x[1:1 + nc].copy(),
        "W": ([unsvec(x[1 + nc + int(sum(wtri[:j])):1 + nc + int(sum(wtri[:j + 1]))], n)
               for j, n in enumerate(extra_psd)] if not psd_maps else
              [unsvec(S @ x[1 + nc + int(sum(wtri[:j])):1 + nc + int(sum(wtri[:j + 1]))], n)
               for j, (n, S) in enumerate(psd_maps)]),
        "Q": Qs,
        "y": z[:nH].copy(),             # dual weights on admissible graphs ("pseudo-densities"), sum ~ 1
        "iterations": int(sol.iterations),
        "solve_time": el,
        "r_prim": float(sol.r_prim),
        "r_dual": float(sol.r_dual),
        "obj_dual": float(sol.obj_val_dual),
        "n_vars": nx,
        "n_rows": int(A.shape[0]),
        "nnz": int(A.nnz),
    }
