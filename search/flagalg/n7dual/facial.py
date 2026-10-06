"""Facial reduction by dead-flag propagation (exact combinatorics, no floats).

For a support set S of admissible 7-vertex graphs and a key k (block (s, m), type sigma), every entry of
M_k(H) is >= 0.  A flag a is DEAD for S if M_k(H)[a, a] = 0 for every H in S.  If Y = sum_H y_H M_k(H) is PSD and
Y[a, a] = 0 then the whole row a of Y vanishes, and since all terms are >= 0, y_H = 0 for every H in S whose M_k(H)
has a nonzero entry in row a.  Removing those H can kill further flags; iterate to a fixed point.  The result is
the largest S' subset of S on which dead flags carry no entry at all, i.e. every PSD-feasible y on S lives on S'
and its moment matrices are zero on the dead rows/columns (so they can be dropped exactly).

Usage (as a module): live_support(T, P, dims, nW) -> (alive_graphs mask, live flag masks per key)
"""
import numpy as np


def sym_entries(T, P):
    """All (graph, a, b) entries of M_k(H) after the automorphism sum (a, b in the full flag index space)."""
    rows, aa, bb = T[0], T[1], T[2]
    G = P.shape[0]
    R = np.repeat(rows[None, :], G, axis=0).ravel()
    A = P[:, aa].ravel()
    B = P[:, bb].ravel()
    return R, A, B


def live_support(T, P, dims, nW, alive=None, verbose=True):
    alive = np.ones(nW, dtype=bool) if alive is None else alive.copy()
    ent = [sym_entries(T[k], P[k]) for k in range(len(dims))]
    it = 0
    while True:
        it += 1
        live = []
        for k, (R, A, B) in enumerate(ent):
            sel = alive[R] & (A == B)
            lv = np.zeros(dims[k], dtype=bool)
            lv[A[sel]] = True
            live.append(lv)
        kill = np.zeros(nW, dtype=bool)
        for k, (R, A, B) in enumerate(ent):
            bad = alive[R] & (~live[k][A] | ~live[k][B])
            kill[R[bad]] = True
        nk = int(kill.sum())
        if verbose:
            print(f"#   facial it {it}: alive {int(alive.sum())}, live flags {[int(l.sum()) for l in live]}, "
                  f"killed {nk}", flush=True)
        if nk == 0:
            return alive, live
        alive &= ~kill
