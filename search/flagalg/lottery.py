"""Lottery flag algebra, first sound formulation: the LOCAL FRACTIONAL PACKING bound.

Problem: l(k, r, p) = lim min |G| / C(n, r) over k-graphs G on [n] such that every p-set meets some edge of G in
>= r points.

Why not a k-graph flag algebra: an optimal G has Theta(n^r) edges, i.e. k-edge density Theta(n^{r-k}) -> 0, so
every fixed-size sample of an optimal G is edgeless with probability -> 1; in the dense limit G is the empty k-graph
and the lottery condition (an existential, non-hereditary condition) cannot be expressed on N-vertex samples.

What is dense: the r-shadow A = d_r(G) (an r-set is in A iff it lies in some edge of G).
 (i)  every p-set S meets an edge E in >= r points  <=>  S contains an r-subset of some E  <=>  S contains an
      A-edge.  So A is an r-graph in which every p-set contains an edge: a hereditary (local) property; this is the
      Turan problem t(p, r).
 (ii) every T in A lies in >= 1 edge of G, and every edge E of G is a K_k^{(r)} of A (all its r-subsets are in A).
Hence for every y >= 0 on A such that sum_{T subset K} y_T <= 1 for every K_k^{(r)} K of A (a fractional packing):
      |G| >= sum_{E in G} sum_{T subset E} y_T >= sum_{T in A} y_T.                                         (*)
Bate's (shadow) bound is y = 1/C(k, r).

Local packing: let y_T depend on the rooted neighbourhood of T in A,
      y_T = E_U [ c( A[T + U], root T ) ],   U a uniformly random (m - r)-subset of V minus T,
with c >= 0 a function on rooted flags (r-graphs on m vertices with a distinguished root edge, roots unordered).
Feasibility, pointwise for every K_k^{(r)} K of A: let U' be a random (m - r)-subset of V minus K; for each T the
law of U differs from U' by total variation O(1/n); for every fixed U', A[K + U'] is an admissible tau-flag F'
(tau = K_k^{(r)} labelled), and if
      sum_{T subset [k], |T| = r} c( F'[T + U'], root T ) <= 1     for every admissible tau-flag F' on k + m - r
                                                                    vertices                           (LP rows)
then sum_{T subset K} y_T <= 1 + O(1/n); rescaling y by 1/(1 + O(1/n)) gives a feasible packing.  Finally
      sum_{T in A} y_T / C(n, r) = E_{T, U}[ 1_A(T) c(A[T + U], T) ] = sum_H p(H; A) rho_c(H)   (exact),
with rho_c(H) = sum_F c_F rho_F(H) and rho_F(H) = Pr over (T, U) inside H of [T in H and (H[T+U], T) ~ F].
So  l(k, r, p) >= min_A sum_H p(H; A) rho_c(H) >= max { b : rho_c(H) - sum <Q_sigma, M_sigma(H)> >= b for all
admissible H, Q_sigma PSD, c >= 0, LP rows },  an SDP jointly in (b, c, Q).  c = 1/C(k,r) recovers Bate with the
same-N Turan bound, so the result is always >= t_FA(N)/C(k, r).

Usage: python lottery.py --k 4 --r 3 --p 4 --N 6 --m 5
"""
from __future__ import annotations

import argparse
import json
import resource
import time
from fractions import Fraction
from itertools import combinations, permutations
from math import comb

import numpy as np
import scipy.sparse as sp

from exact import certify_eig
from flags import build_problem
from hypergraphs import (EveryPSetHasEdge, adjacency_flat, canon_batch, colex_index, edges, induced_mask_index,
                         mask_to_edges)
from sdp import solve_flag_sdp


def rooted_flags(m: int, r: int, pred):
    """Rooted flags on m vertices: masks with bit 0 (= root edge {0..r-1}) set; canonical under S_r x S_{m-r}."""
    nb = comb(m, r)
    free = np.arange(1 << (nb - 1), dtype=np.int64)
    masks = (free << 1) | 1
    ok = pred(masks, m, r)
    can = canon_batch(masks[ok], m, r, kind="setwise", s=r)
    flags, inv = np.unique(can, return_inverse=True)
    lookup = np.full(free.size, -1, dtype=np.int64)
    lookup[np.nonzero(ok)[0]] = inv
    return flags, lookup          # lookup[mask >> 1] -> rooted flag index


def rho_matrix(H: np.ndarray, N: int, r: int, m: int, lookup: np.ndarray, nF: int):
    """Exact counts rho_F(H) * C(N, r) * C(N - r, m - r)  as a sparse (nH, nF) integer matrix, and the denominator."""
    tuples = []
    for T in combinations(range(N), r):
        rest = [v for v in range(N) if v not in T]
        for U in combinations(rest, m - r):
            tuples.append(T + U)
    tuples = np.array(tuples, dtype=np.int64)
    I = induced_mask_index(tuples, N, r)
    pw = np.int64(1) << np.arange(I.shape[1], dtype=np.int64)
    rows, cols = [], []
    bs = max(1, int(1e7 // max(1, I.size)))
    for i in range(0, H.size, bs):
        A = adjacency_flat(H[i:i + bs], N, r)
        M = (A[:, I].astype(np.int64) * pw).sum(axis=2)         # (G, K)
        g, kk = np.nonzero(M & 1)                                 # root T is an edge
        j = lookup[M[g, kk] >> 1]
        assert (j >= 0).all()
        rows.append(g + i)
        cols.append(j)
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    Rm = sp.coo_matrix((np.ones(rows.size, dtype=np.int64), (rows, cols)), shape=(H.size, nF)).tocsr()
    Rm.sum_duplicates()
    return Rm, comb(N, r) * comb(N - r, m - r)


def tau_constraints(k: int, r: int, m: int, pred, lookup: np.ndarray):
    """Rows (as lists of rooted-flag indices, with multiplicity) of the pointwise packing constraints, one per
    admissible tau-flag (tau = K_k^{(r)} on [k]) on k + m - r vertices, up to S_k x S_{m-r}."""
    v = k + m - r
    low = comb(k, r)
    nb = comb(v, r)
    full = (1 << low) - 1
    free = np.arange(1 << (nb - low), dtype=np.int64)
    masks = np.int64(full) | (free << low)
    ok = pred(masks, v, r)
    can = np.unique(canon_batch(masks[ok], v, r, kind="setwise", s=k))
    U = tuple(range(k, v))
    tuples = np.array([T + U for T in combinations(range(k), r)], dtype=np.int64)
    I = induced_mask_index(tuples, v, r)                          # (C(k,r), C(m,r))
    pw = np.int64(1) << np.arange(I.shape[1], dtype=np.int64)
    A = adjacency_flat(can, v, r)
    M = (A[:, I].astype(np.int64) * pw).sum(axis=2)              # (nTau, C(k,r))
    assert (M & 1).all()
    J = lookup[M >> 1]
    assert (J >= 0).all()
    return can, J                                                 # J[t, i] = rooted flag of the i-th r-subset


def tau_sos_rows(k: int, r: int, m: int, pred, lookup: np.ndarray):
    """Averaged packing rows strengthened by a tau-rooted SOS term (tau = K_k^{(r)} labelled on [k]).

    For tau-flags F' on k + 2 vertices (K fixed pointwise, u1 = k, u2 = k + 1, canonical under the swap u1 <-> u2):
        cpart(F') + W[a(F'), b(F')] <= 1,
        cpart(F') = average over U subset {u1, u2}, |U| = m - r, of sum_{T subset [k], |T| = r} c(F'[T + U], T),
    with a, b the tau-flags K + u1, K + u2 (on k + 1 vertices) and W PSD.  For every K_k^{(r)} K of A, averaging over
    a random pair (u1, u2) outside K gives  sum_{T subset K} y_T <= 1 - E[W(a(u1), b(u2))] + O(1/n)
    = 1 - v_K^T W v_K + O(1/n) <= 1 + O(1/n), v_K = law of the link of a random vertex relative to K.
    Requires m - r <= 2.  Returns (ccount (nRows, nF) sparse int, cden, a, b, n_tau1)."""
    assert m - r <= 2
    v = k + 2
    low = comb(k, r)
    full = (1 << low) - 1
    # tau-flags on k + 1 vertices: K + u, link bits = r-sets containing u; labelled, no symmetry
    nb1 = comb(k + 1, r)
    free1 = np.arange(1 << (nb1 - low), dtype=np.int64)
    m1 = np.int64(full) | (free1 << low)
    ok1 = pred(m1, k + 1, r)
    lookup1 = np.full(free1.size, -1, dtype=np.int64)
    lookup1[np.nonzero(ok1)[0]] = np.arange(int(ok1.sum()))
    n1 = int(ok1.sum())
    nb = comb(v, r)
    free = np.arange(1 << (nb - low), dtype=np.int64)
    masks = np.int64(full) | (free << low)
    ok = pred(masks, v, r)
    can = np.unique(canon_batch(masks[ok], v, r, kind="fix", s=k))
    A = adjacency_flat(can, v, r)
    K = tuple(range(k))
    ta = np.array([K + (k,)], dtype=np.int64)
    tb = np.array([K + (k + 1,)], dtype=np.int64)
    pw1 = np.int64(1) << np.arange(nb1, dtype=np.int64)
    Ma = (A[:, induced_mask_index(ta, v, r)].astype(np.int64) * pw1).sum(axis=2)[:, 0]
    Mb = (A[:, induced_mask_index(tb, v, r)].astype(np.int64) * pw1).sum(axis=2)[:, 0]
    a = lookup1[Ma >> low]
    b = lookup1[Mb >> low]
    assert (a >= 0).all() and (b >= 0).all()
    Us = list(combinations((k, k + 1), m - r))
    tuples = np.array([T + U for U in Us for T in combinations(range(k), r)], dtype=np.int64)
    I = induced_mask_index(tuples, v, r)
    pw = np.int64(1) << np.arange(I.shape[1], dtype=np.int64)
    M = (A[:, I].astype(np.int64) * pw).sum(axis=2)
    J = lookup[M >> 1]
    assert (J >= 0).all() and (M & 1).all()
    rr = np.repeat(np.arange(can.size), J.shape[1])
    cc = sp.coo_matrix((np.ones(rr.size, dtype=np.int64), (rr, J.ravel())), shape=(can.size, lookup.max() + 1))
    cc = cc.tocsr()
    cc.sum_duplicates()
    return cc, len(Us), a, b, n1, np.nonzero(ok1)[0]


def tau_sos_rows_sym(k: int, r: int, m: int, pred, lookup: np.ndarray):
    """S_k-symmetry-reduced version of tau_sos_rows.

    The packing rows are permuted among themselves by relabelling K (S_k) and swapping u1, u2, so W may be taken
    S_k-invariant: W[pi a, pi b] = W[a, b] for the induced action of pi on the tau-flags K + u (links).  For invariant
    W, one row per tau-flag class under S_k x S_2 suffices.  W is parametrised by its values on the orbits of S_k on
    unordered pairs {a, b}:  W = sum_o w_o E_o.
    Returns (Cm, cdiv, a, b, n1, tau1_links, orbit (n1 x n1 int), n_orb, link_perms (list of index permutations))."""
    assert m - r <= 2
    v = k + 2
    low = comb(k, r)
    full = (1 << low) - 1
    nb1 = comb(k + 1, r)
    free1 = np.arange(1 << (nb1 - low), dtype=np.int64)
    m1 = np.int64(full) | (free1 << low)
    ok1 = pred(m1, k + 1, r)
    links1 = np.nonzero(ok1)[0]
    lookup1 = np.full(free1.size, -1, dtype=np.int64)
    lookup1[links1] = np.arange(links1.size)
    n1 = int(links1.size)
    # action of S_k on link masks (bits = colex order of the (r-1)-subsets of [k])
    sub = list(combinations(range(k), r - 1))
    sub_idx = {e: i for i, e in enumerate(sorted(sub, key=lambda e: e[::-1]))}
    link_perms = []
    for pi in permutations(range(k)):
        img = np.empty(n1, dtype=np.int64)
        for ai, L in enumerate(links1):
            L2 = 0
            for e, i in sub_idx.items():
                if (int(L) >> i) & 1:
                    L2 |= 1 << sub_idx[tuple(sorted(pi[x] for x in e))]
            img[ai] = lookup1[L2]
        assert (img >= 0).all()
        link_perms.append(img)
    orbit = np.full((n1, n1), -1, dtype=np.int64)
    n_orb = 0
    for a0 in range(n1):
        for b0 in range(a0, n1):
            if orbit[a0, b0] >= 0:
                continue
            for img in link_perms:
                x, y = img[a0], img[b0]
                orbit[x, y] = orbit[y, x] = n_orb
            n_orb += 1
    nb = comb(v, r)
    free = np.arange(1 << (nb - low), dtype=np.int64)
    masks = np.int64(full) | (free << low)
    ok = pred(masks, v, r)
    can = np.unique(canon_batch(masks[ok], v, r, kind="setwise", s=k))
    A = adjacency_flat(can, v, r)
    K = tuple(range(k))
    pw1 = np.int64(1) << np.arange(nb1, dtype=np.int64)
    Ma = (A[:, induced_mask_index(np.array([K + (k,)]), v, r)].astype(np.int64) * pw1).sum(axis=2)[:, 0]
    Mb = (A[:, induced_mask_index(np.array([K + (k + 1,)]), v, r)].astype(np.int64) * pw1).sum(axis=2)[:, 0]
    a = lookup1[Ma >> low]
    b = lookup1[Mb >> low]
    Us = list(combinations((k, k + 1), m - r))
    tuples = np.array([T + U for U in Us for T in combinations(range(k), r)], dtype=np.int64)
    I = induced_mask_index(tuples, v, r)
    pw = np.int64(1) << np.arange(I.shape[1], dtype=np.int64)
    M = (A[:, I].astype(np.int64) * pw).sum(axis=2)
    J = lookup[M >> 1]
    assert (J >= 0).all() and (M & 1).all()
    rr = np.repeat(np.arange(can.size), J.shape[1])
    cc = sp.coo_matrix((np.ones(rr.size, dtype=np.int64), (rr, J.ravel())), shape=(can.size, lookup.max() + 1))
    cc = cc.tocsr()
    cc.sum_duplicates()
    return cc, len(Us), a, b, n1, links1, orbit, n_orb, link_perms


def symmetrize_exact(Wnum, link_perms):
    """Exact average of P_pi W P_pi^T over the group (object array of ints) -> (num, den)."""
    n = Wnum.shape[0]
    acc = np.zeros((n, n), dtype=object)
    acc[:, :] = 0
    for img in link_perms:
        acc = acc + Wnum[np.ix_(np.argsort(img), np.argsort(img))]
    return acc, len(link_perms)


def build_lottery(k: int, r: int, p: int, N: int, m: int, tau_sos: bool = False, verbose: bool = True,
                  tau_sym: bool = False) -> dict:
    """Assemble the lottery local-packing SDP data (independent of the solver)."""
    pred = EveryPSetHasEdge(p)
    prob = build_problem(N, r, pred, type_mode="standard", verbose=verbose)
    flags, lookup = rooted_flags(m, r, pred)
    nF = flags.size
    Rm, rden = rho_matrix(prob.H, N, r, m, lookup, nF)
    ctx = {"k": k, "r": r, "p": p, "N": N, "m": m, "tau_sos": tau_sos, "prob": prob, "flags": flags,
           "lookup": lookup, "nF": nF, "Rm": Rm, "rden": rden, "tau1_links": None}
    ctx["tau_sym"] = tau_sym
    if tau_sym:
        assert tau_sos
        Cm, cdiv, ta, tb, n1, tau1_links, orbit, n_orb, link_perms = tau_sos_rows_sym(k, r, m, pred, lookup)
        nT = Cm.shape[0]
        from sdp import svec_index
        # G columns: [c (nF) | w (n_orb)];  row t: cpart + w_{orbit(a_t, b_t)} <= 1
        Wc = sp.coo_matrix((np.ones(nT), (np.arange(nT), nF + orbit[ta, tb])), shape=(nT, nF + n_orb))
        Gm = (sp.hstack([Cm.astype(float) / cdiv, sp.csr_matrix((nT, n_orb))]) + Wc).tocsr()
        ii, jj = np.triu_indices(n1)
        S = sp.coo_matrix((np.where(ii == jj, 1.0, np.sqrt(2.0)), (svec_index(ii, jj), orbit[ii, jj])),
                          shape=(n1 * (n1 + 1) // 2, n_orb)).tocsr()
        ctx.update({"Cm": Cm, "cdiv": cdiv, "ta": ta, "tb": tb, "n1": n1, "Gm": Gm, "extra": (),
                    "psd_maps": ((n1, S),), "tau1_links": tau1_links, "orbit": orbit, "n_orb": n_orb,
                    "link_perms": link_perms})
    elif not tau_sos:
        taus, J = tau_constraints(k, r, m, pred, lookup)
        nT = taus.size
        gr, gc = np.repeat(np.arange(nT), J.shape[1]), J.ravel()
        Cm = sp.coo_matrix((np.ones(gr.size, dtype=np.int64), (gr, gc)), shape=(nT, nF)).tocsr()
        Cm.sum_duplicates()
        ctx.update({"Cm": Cm, "cdiv": 1, "ta": None, "tb": None, "n1": 0, "Gm": Cm.astype(float), "extra": ()})
    else:
        Cm, cdiv, ta, tb, n1, tau1_links = tau_sos_rows(k, r, m, pred, lookup)
        nT = Cm.shape[0]
        wcoef = np.where(ta == tb, 1.0, 1.0 / np.sqrt(2.0))
        from sdp import svec_index
        Wc = sp.coo_matrix((wcoef, (np.arange(nT), nF + svec_index(ta, tb))),
                           shape=(nT, nF + n1 * (n1 + 1) // 2))
        Gm = (sp.hstack([Cm.astype(float) / cdiv, sp.csr_matrix((nT, n1 * (n1 + 1) // 2))]) + Wc).tocsr()
        ctx.update({"Cm": Cm, "cdiv": cdiv, "ta": ta, "tb": tb, "n1": n1, "Gm": Gm, "extra": (n1,),
                    "tau1_links": tau1_links})
    ctx["nT"] = nT
    if verbose:
        print(f"  rooted flags on m={m}: {nF}; packing rows: {nT} (tau_sos={tau_sos}, W block {ctx['n1']})",
              flush=True)
    return ctx


def certify_lottery(ctx: dict, sol: dict, bits: int = 40, cert_out: str | None = None) -> dict:
    """Exact rounding: c -> nonneg dyadic c_r; W -> exactly PSD W_r (eig); every packing row evaluated exactly;
    alpha = 1/max(1, max row) scales (c_r, W_r, Q_sigma) so that all rows are <= 1 exactly;
    bound = alpha * min_H (rho_{c_r}(H) - <Q_r, M(H)>), Q_r exactly PSD (eig rounding)."""
    prob, Cm, cdiv, nT = ctx["prob"], ctx["Cm"], ctx["cdiv"], ctx["nT"]
    scale = 1 << bits
    cnum = [max(0, int(round(float(x) * scale))) for x in sol["c"]]
    Cl = Cm.tocoo()
    rowsum = [0] * nT
    for i, j, vv in zip(Cl.row.tolist(), Cl.col.tolist(), Cl.data.tolist()):
        rowsum[i] += int(vv) * cnum[j]
    rowval = [Fraction(x, cdiv * scale) for x in rowsum]
    w_psd = None
    Wnum, wsh, wden = None, 0, 1
    if ctx["tau_sos"]:
        from exact import round_psd_eig, is_psd_exact
        Wnum, wsh = round_psd_eig(sol["W"][0], bits)
        wden = 1 << wsh
        if ctx.get("tau_sym"):      # exact group average: still PSD, now invariant, so canonical rows suffice
            Wnum, g = symmetrize_exact(Wnum, ctx["link_perms"])
            wden *= g
        w_psd = is_psd_exact(Wnum.tolist())[:2]
        Wl = Wnum.tolist()
        rowval = [rv + Fraction(int(Wl[int(a_)][int(b_)]), wden) for rv, a_, b_ in zip(rowval, ctx["ta"], ctx["tb"])]
    rowmax = max(rowval) if rowval else Fraction(0)
    alpha = 1 / max(Fraction(1), rowmax)
    assert all(alpha * rv <= 1 for rv in rowval)
    Rl = ctx["Rm"].tocoo()
    objnum = [0] * prob.nH
    for i, j, vv in zip(Rl.row.tolist(), Rl.col.tolist(), Rl.data.tolist()):
        objnum[i] += int(vv) * cnum[j]
    obj_exact = [Fraction(x, ctx["rden"] * scale) for x in objnum]
    ce = certify_eig(prob, obj_exact, sol["Q"], bits=bits)
    ce["b_exact"] = alpha * ce["b_exact"]
    if cert_out:
        cert = {"k": ctx["k"], "r": ctx["r"], "p": ctx["p"], "N": ctx["N"], "m": ctx["m"], "tau_sos": ctx["tau_sos"],
                "b_exact": str(ce["b_exact"]), "alpha": str(alpha), "bits": bits,
                "rooted_flags": [int(x) for x in ctx["flags"]],
                "c_exact": [str(alpha * Fraction(x, scale)) for x in cnum],
                "Q_sigma": [{"s": blk.fam.s, "m": blk.fam.m, "sigma": int(blk.fam.sigma),
                             "flags": [int(x) for x in blk.fam.flags], "den": str(den / alpha),
                             "Qnum": [[str(int(x)) for x in row] for row in Qn.tolist()]}
                            for blk, (Qn, den) in zip(prob.blocks, ce["Qexact"])]}
        if ctx["tau_sos"]:
            cert["tau1_links"] = [int(x) for x in ctx["tau1_links"]]
            cert["W_den"] = str(Fraction(wden) / alpha)
            cert["W_num"] = [[str(int(x)) for x in row] for row in Wnum.tolist()]
        with open(cert_out, "w") as f:
            json.dump(cert, f)
    return {"b_exact": ce["b_exact"], "psd": ce["psd"], "alpha": alpha, "rowmax": rowmax, "w_psd": w_psd}


def run(k: int, r: int, p: int, N: int, m: int, bits: int = 40, tol: float = 1e-10, verbose: bool = True,
        tau_sos: bool = False, cert_out: str | None = None, tau_sym: bool = False, settings=None) -> dict:
    t0 = time.time()
    ctx = build_lottery(k, r, p, N, m, tau_sos, verbose, tau_sym=tau_sym)
    prob = ctx["prob"]
    t_build = time.time() - t0
    t1 = time.time()
    sol = solve_flag_sdp(prob, np.zeros(prob.nH), R=(ctx["Rm"] / ctx["rden"]), G=ctx["Gm"], h=np.ones(ctx["nT"]),
                         extra_psd=ctx["extra"], tol=tol, psd_maps=ctx.get("psd_maps", ()), settings=settings)
    t_sdp = time.time() - t1
    if verbose:
        print(f"  SDP: {sol['status']} b = {sol['b']:.12f} ({sol['iterations']} it, {t_sdp:.2f}s)", flush=True)
    # Bate check with the same N: Turan SDP alone
    num, den = prob.density_num_den()
    t2 = time.time()
    solT = solve_flag_sdp(prob, num / den, tol=tol, settings=settings)
    t_turan = time.time() - t2
    t3 = time.time()
    cl = certify_lottery(ctx, sol, bits, cert_out)
    t_eig = time.time() - t3
    flags = ctx["flags"]
    order = np.argsort(-sol["c"])
    top = [{"c": float(sol["c"][i]), "rooted_flag_edges": mask_to_edges(int(flags[i]), m, r)} for i in order[:6]]
    out = {
        "problem": f"lottery l({k},{r},{p}) local fractional packing (y_T from rooted {m}-vertex flags), N={N}",
        "k": k, "r": r, "p": p, "N": N, "m": m,
        "n_admissible": prob.nH,
        "blocks": [{"s": b.fam.s, "m": b.fam.m, "n_flags": b.size} for b in prob.blocks],
        "n_rooted_flags": int(ctx["nF"]), "n_tau_constraints": int(ctx["nT"]),
        "sdp": {kk: sol[kk] for kk in ("status", "b", "iterations", "solve_time", "r_prim", "r_dual", "obj_dual",
                                       "n_vars", "n_rows", "nnz")},
        "turan_same_N": {"status": solT["status"], "t_lower": solT["b"], "bate": solT["b"] / comb(k, r)},
        "gain_over_bate_same_N": sol["b"] - solT["b"] / comb(k, r),
        "eig_rounding": {"b_exact": str(cl["b_exact"]), "b_exact_float": float(cl["b_exact"]),
                         "psd_checks": [{"psd": ok, "rank": rk} for ok, rk in cl["psd"]],
                         "alpha": float(cl["alpha"]), "max_packing_row_before_scaling": float(cl["rowmax"]),
                         "W_psd_check": cl["w_psd"]},
        "tau_sos": tau_sos, "tau_sym": tau_sym,
        "top_c": top,
        "c_min_max": [float(sol["c"].min()), float(sol["c"].max())],
        "time_s": {"build": t_build, "sdp": t_sdp, "turan_sdp": t_turan, "eig_rounding": t_eig,
                   "total": time.time() - t0},
        "peak_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0,
    }
    if verbose:
        print(f"  Turan same N: {solT['b']:.10f}  Bate {solT['b'] / comb(k, r):.10f}  "
              f"lottery {sol['b']:.10f}  exact {float(cl['b_exact']):.10f}", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--r", type=int, required=True)
    ap.add_argument("--p", type=int, required=True)
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--bits", type=int, default=40)
    ap.add_argument("--tol", type=float, default=1e-10)
    ap.add_argument("--out", default=None)
    ap.add_argument("--tau-sos", action="store_true")
    ap.add_argument("--cert", default=None)
    ap.add_argument("--tau-sym", action="store_true", help="S_k-symmetry-reduced tau-SOS rows (implies --tau-sos)")
    ap.add_argument("--static-reg", type=float, default=None)
    a = ap.parse_args()
    tau_sos = a.tau_sos or a.tau_sym
    st = {"static_regularization_constant": a.static_reg} if a.static_reg else None
    print(f"# lottery k={a.k} r={a.r} p={a.p} N={a.N} m={a.m} tau_sos={tau_sos} tau_sym={a.tau_sym}", flush=True)
    res = run(a.k, a.r, a.p, a.N, a.m, a.bits, a.tol, tau_sos=tau_sos, cert_out=a.cert, tau_sym=a.tau_sym,
              settings=st)
    js = json.dumps(res, indent=1, default=str)
    if a.out:
        with open(a.out, "w") as f:
            f.write(js)
    print(js)


if __name__ == "__main__":
    main()
