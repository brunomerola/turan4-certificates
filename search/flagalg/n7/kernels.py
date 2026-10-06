"""numba kernels for the N = 7, r = 4 flag problems (tables from common.build_model).

Fast path for the three blocks whose flags are a single 4-set (s, m) = (1,4), (2,4), (3,4): their configurations
are exactly the unordered pairs {e, f} of 4-sets with |e & f| = s (e & f = the root set, e - f and f - e the two
free parts), so with x_e = [e in H]
    sum_terms wint * Qsum[x_e, x_f] = wint * (Qsum11 A_s + Qsum01 B_s + Qsum00 C_s),
A_s / B_s / C_s = numbers of edge-edge / edge-nonedge / nonedge-nonedge pairs with |e & f| = s (hardware popcounts
against precomputed 35-bit neighbourhood masks).  cdd(H) = A_1.

Per 4-graph H on 7 vertices (35-bit mask) the kernels compute
    e(H)                      number of edges,
    cdd(H)                    #{(v, {U1, U2}) : v + U1, v + U2 both edges}, U1 | U2 a split of V - v into triples
                              (so  [[d . d]]_1 (H) = 72 cdd / 5040, the stationarity term),
    Z(H) = sum_terms wint * Qsum[a, b]   with  <Q, M(H)> = Z / 5040   (see common.py),
and the float slack
    slack(H) = d(H) - Z/5040 - tau1 * (b0 d - c_dd) - tau2 * (b0 (1 - d) - (d - c_dd)),   d = e/35, c_dd = 72 cdd/5040.
The exact kernels use int64 Qsum and return, per group (e, cdd), the maximum of Z and an H attaining it (the exact
slack is affine in (e, cdd) and decreasing in Z, so the group maxima determine the exact minimum).
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
from llvmlite import ir
from numba import njit, prange, types
from numba.extending import intrinsic

INT64_MIN = np.iinfo(np.int64).min


@intrinsic
def popcnt64(typingctx, x):
    """llvm.ctpop.i64 (hardware popcount)."""
    sig = types.int64(types.int64)

    def codegen(context, builder, signature, args):
        fn = builder.module.declare_intrinsic("llvm.ctpop", [ir.IntType(64)])
        return builder.call(fn, args)
    return sig, codegen


@njit(cache=True)
def _unpack(h, bits):
    for i in range(35):
        bits[i] = np.uint8((h >> i) & 1)


@njit(cache=True)
def _eval_Z(h, bits, nb, nrs, nsig, nU, nfree, npairs, wint, sigpos, freepos, pairs, type_of, flag_of, qoff, qdim,
            aux, Q, fl):
    """Returns (Z, bad, e, cdd).  aux rows 0..2: 35-bit masks of the 4-sets meeting edge i in exactly 1, 2, 3
    vertices; aux[3, bi] = s for a fast-path block (1-bit flags), 0 otherwise."""
    Z = Q[0] * 0
    bad = 0
    for bi in range(nb):
        if aux[3, bi] > 0:
            continue
        for rs in range(nrs[bi]):
            sc = 0
            for k in range(nsig[bi]):
                sc |= np.int64(bits[sigpos[bi, rs, k]]) << k
            t = type_of[bi, sc]
            if t < 0:
                bad = 1
                continue
            ok = True
            for u in range(nU[bi]):
                fc = 0
                for k in range(nfree[bi]):
                    fc |= np.int64(bits[freepos[bi, rs, u, k]]) << k
                f = flag_of[bi, sc, fc]
                fl[u] = f
                if f < 0:
                    ok = False
            if not ok:
                bad = 1
                continue
            base = qoff[bi, t]
            n = qdim[bi, t]
            acc = Q[0] * 0
            for q in range(npairs[bi]):
                acc += Q[base + fl[pairs[bi, q, 0]] * n + fl[pairs[bi, q, 1]]]
            Z += wint[bi] * acc
    # pair statistics by intersection size j = 1, 2, 3
    A1 = 0
    A2 = 0
    A3 = 0
    B1 = 0
    B2 = 0
    B3 = 0
    for i in range(35):
        if (h >> i) & 1:
            A1 += popcnt64(h & aux[0, i])
            A2 += popcnt64(h & aux[1, i])
            A3 += popcnt64(h & aux[2, i])
        else:
            B1 += popcnt64(h & aux[0, i])
            B2 += popcnt64(h & aux[1, i])
            B3 += popcnt64(h & aux[2, i])
    A1 //= 2
    A2 //= 2
    A3 //= 2
    for bi in range(nb):
        j = aux[3, bi]
        if j == 0:
            continue
        base = qoff[bi, 0]
        if j == 1:
            A, B, C = A1, B1, 70 - A1 - B1
        elif j == 2:
            A, B, C = A2, B2, 315 - A2 - B2
        else:
            A, B, C = A3, B3, 210 - A3 - B3
        Z += wint[bi] * (Q[base + 3] * A + Q[base + 1] * B + Q[base] * C)
    return Z, bad, popcnt64(h), A1


@njit(cache=True)
def _cdd(bits, cdd_pos):
    c = 0
    for i in range(cdd_pos.shape[0]):
        c += np.int64(bits[cdd_pos[i, 0]] & bits[cdd_pos[i, 1]])
    return c


@njit(cache=True)
def _popc(bits):
    c = 0
    for i in range(35):
        c += np.int64(bits[i])
    return c


@njit(parallel=True, cache=True)
def scan_float(H, nb, nrs, nsig, nU, nfree, npairs, wint, sigpos, freepos, pairs, type_of, flag_of, qoff, qdim,
               aux, Q, tau1, tau2, b0, nchunks):
    """slack(H) for every H (float64); bad[i] = 1 if H[i] has an inadmissible type/flag (not admissible)."""
    n = H.size
    out = np.empty(n, dtype=np.float64)
    bad = np.zeros(n, dtype=np.uint8)
    step = (n + nchunks - 1) // nchunks
    for c in prange(nchunks):
        bits = np.empty(35, dtype=np.uint8)
        fl = np.empty(20, dtype=np.int64)
        for i in range(c * step, min(n, (c + 1) * step)):
            _unpack(H[i], bits)
            Z, bd, e, cd = _eval_Z(H[i], bits, nb, nrs, nsig, nU, nfree, npairs, wint, sigpos, freepos, pairs,
                                   type_of, flag_of, qoff, qdim, aux, Q, fl)
            d = e / 35.0
            c2 = 72.0 * cd / 5040.0
            out[i] = d - Z / 5040.0 - tau1 * (b0 * d - c2) - tau2 * (b0 * (1.0 - d) - (d - c2))
            bad[i] = bd
    return out, bad


@njit(parallel=True, cache=True)
def scan_exact(H, nb, nrs, nsig, nU, nfree, npairs, wint, sigpos, freepos, pairs, type_of, flag_of, qoff, qdim,
               aux, Qi, nchunks):
    """Exact (int64) group maxima of Z over the given masks.  Returns (maxZ[36, 71], argH[36, 71], count[36, 71],
    nbad)."""
    n = H.size
    step = (n + nchunks - 1) // nchunks
    mz = np.full((nchunks, 36, 71), INT64_MIN, dtype=np.int64)
    ah = np.zeros((nchunks, 36, 71), dtype=np.int64)
    cnt = np.zeros((nchunks, 36, 71), dtype=np.int64)
    nbad = np.zeros(nchunks, dtype=np.int64)
    for c in prange(nchunks):
        bits = np.empty(35, dtype=np.uint8)
        fl = np.empty(20, dtype=np.int64)
        for i in range(c * step, min(n, (c + 1) * step)):
            _unpack(H[i], bits)
            Z, bd, e, cd = _eval_Z(H[i], bits, nb, nrs, nsig, nU, nfree, npairs, wint, sigpos, freepos, pairs,
                                   type_of, flag_of, qoff, qdim, aux, Qi, fl)
            if bd:
                nbad[c] += 1
                continue
            cnt[c, e, cd] += 1
            if Z > mz[c, e, cd]:
                mz[c, e, cd] = Z
                ah[c, e, cd] = H[i]
    return _merge(mz, ah, cnt, nbad)


@njit(cache=True)
def _merge(mz, ah, cnt, nbad):
    MZ = np.full((36, 71), INT64_MIN, dtype=np.int64)
    AH = np.zeros((36, 71), dtype=np.int64)
    CN = np.zeros((36, 71), dtype=np.int64)
    for c in range(mz.shape[0]):
        for e in range(36):
            for d in range(71):
                CN[e, d] += cnt[c, e, d]
                if mz[c, e, d] > MZ[e, d]:
                    MZ[e, d] = mz[c, e, d]
                    AH[e, d] = ah[c, e, d]
    return MZ, AH, CN, nbad.sum()


@njit(parallel=True, cache=True)
def scan_raw_exact(reps, psets_new, nb, nrs, nsig, nU, nfree, npairs, wint, sigpos, freepos, pairs, type_of,
                   flag_of, qoff, qdim, aux, Qi, lchunk):
    """Exact group maxima over ALL one-vertex extensions R | (L << 15), R in reps, L in [0, 2^20), that satisfy the
    p-set condition on the p-sets through vertex 6 (no canonical forms involved)."""
    nper = (1 << 20) // lchunk
    ntask = reps.size * nper
    mz = np.full((ntask, 36, 71), INT64_MIN, dtype=np.int64)
    ah = np.zeros((ntask, 36, 71), dtype=np.int64)
    cnt = np.zeros((ntask, 36, 71), dtype=np.int64)
    nbad = np.zeros(ntask, dtype=np.int64)
    for task in prange(ntask):
        ri = task // nper
        l0 = (task % nper) * lchunk
        Rm = reps[ri]
        bits = np.empty(35, dtype=np.uint8)
        fl = np.empty(20, dtype=np.int64)
        for L in range(l0, l0 + lchunk):
            h = Rm | (L << 15)
            ok = True
            for j in range(psets_new.size):
                if (h & psets_new[j]) == 0:
                    ok = False
                    break
            if not ok:
                continue
            _unpack(h, bits)
            Z, bd, e, cd = _eval_Z(h, bits, nb, nrs, nsig, nU, nfree, npairs, wint, sigpos, freepos, pairs,
                                   type_of, flag_of, qoff, qdim, aux, Qi, fl)
            if bd:
                nbad[task] += 1
                continue
            cnt[task, e, cd] += 1
            if Z > mz[task, e, cd]:
                mz[task, e, cd] = Z
                ah[task, e, cd] = h
    return _merge(mz, ah, cnt, nbad)


@njit(parallel=True, cache=True)
def scan_raw_float_min(reps, psets_new, nb, nrs, nsig, nU, nfree, npairs, wint, sigpos, freepos, pairs, type_of,
                       flag_of, qoff, qdim, aux, Q, tau1, tau2, b0, lchunk):
    """min slack over all raw extensions (float); returns (min per task, argmin per task, n admissible)."""
    nper = (1 << 20) // lchunk
    ntask = reps.size * nper
    mn = np.full(ntask, np.inf)
    am = np.zeros(ntask, dtype=np.int64)
    na = np.zeros(ntask, dtype=np.int64)
    for task in prange(ntask):
        ri = task // nper
        l0 = (task % nper) * lchunk
        Rm = reps[ri]
        bits = np.empty(35, dtype=np.uint8)
        fl = np.empty(20, dtype=np.int64)
        for L in range(l0, l0 + lchunk):
            h = Rm | (L << 15)
            ok = True
            for j in range(psets_new.size):
                if (h & psets_new[j]) == 0:
                    ok = False
                    break
            if not ok:
                continue
            na[task] += 1
            _unpack(h, bits)
            Z, bd, e, cd = _eval_Z(h, bits, nb, nrs, nsig, nU, nfree, npairs, wint, sigpos, freepos, pairs,
                                   type_of, flag_of, qoff, qdim, aux, Q, fl)
            d = e / 35.0
            c2 = 72.0 * cd / 5040.0
            s = d - Z / 5040.0 - tau1 * (b0 * d - c2) - tau2 * (b0 * (1.0 - d) - (d - c2))
            if s < mn[task]:
                mn[task] = s
                am[task] = h
    return mn, am, na


@njit(cache=True)
def terms(H, nb, nrs, nsig, nU, nfree, npairs, wint, sigpos, freepos, pairs, type_of, flag_of, keyid, cdd_pos):
    """All terms of every H: (row, key, a, b, wint) with a <= b, plus e(H) and cdd(H)."""
    n = H.size
    tot = 0
    for bi in range(nb):
        tot += nrs[bi] * npairs[bi]
    rows = np.empty(n * tot, dtype=np.int64)
    keys = np.empty(n * tot, dtype=np.int64)
    aa = np.empty(n * tot, dtype=np.int64)
    bb = np.empty(n * tot, dtype=np.int64)
    ww = np.empty(n * tot, dtype=np.int64)
    ee = np.empty(n, dtype=np.int64)
    cc = np.empty(n, dtype=np.int64)
    bits = np.empty(35, dtype=np.uint8)
    fl = np.empty(20, dtype=np.int64)
    k = 0
    for i in range(n):
        _unpack(H[i], bits)
        ee[i] = _popc(bits)
        cc[i] = _cdd(bits, cdd_pos)
        for bi in range(nb):
            for rs in range(nrs[bi]):
                sc = 0
                for j in range(nsig[bi]):
                    sc |= np.int64(bits[sigpos[bi, rs, j]]) << j
                t = type_of[bi, sc]
                for u in range(nU[bi]):
                    fc = 0
                    for j in range(nfree[bi]):
                        fc |= np.int64(bits[freepos[bi, rs, u, j]]) << j
                    fl[u] = flag_of[bi, sc, fc]
                for q in range(npairs[bi]):
                    a = fl[pairs[bi, q, 0]]
                    b = fl[pairs[bi, q, 1]]
                    rows[k] = i
                    keys[k] = keyid[bi, t]
                    aa[k] = min(a, b)
                    bb[k] = max(a, b)
                    ww[k] = wint[bi]
                    k += 1
    return rows[:k], keys[:k], aa[:k], bb[:k], ww[:k], ee, cc


@njit(parallel=True, cache=True)
def coef_columns(rows, aa, bb, ww, P, V, nrows):
    """out[c, r] = (1/5040) sum_{terms i of row r} ww[i] * sum_g V[c, P[g, aa[i]]] * V[c, P[g, bb[i]]]."""
    k = V.shape[0]
    out = np.zeros((k, nrows))
    G = P.shape[0]
    for c in prange(k):
        for i in range(rows.size):
            a = aa[i]
            b = bb[i]
            acc = 0.0
            for g in range(G):
                acc += V[c, P[g, a]] * V[c, P[g, b]]
            out[c, rows[i]] += ww[i] * acc
    return out / 5040.0


@njit(cache=True)
def moment_raw(rows, aa, bb, ww, y, n):
    """R = sum_terms y[row] * ww / 5040 * (E_ab + E_ba) / 2  (to be summed over the type's automorphisms)."""
    Rm = np.zeros((n, n))
    for i in range(rows.size):
        v = y[rows[i]] * ww[i] / 5040.0 * 0.5
        Rm[aa[i], bb[i]] += v
        Rm[bb[i], aa[i]] += v
    return Rm


def tables(model):
    """Positional kernel arguments (after H) shared by all scans (compact dtypes: positions < 35 and flag indices
    < 1024 fit int8/int16, which keeps the tables in L1/L2)."""
    if getattr(model, "_tab", None) is None:
        m = model
        model._tab = (len(m.blocks), m.nrs, m.nsig, m.nU, m.nfree, m.npairs, m.wint,
                           np.ascontiguousarray(m.sigpos, dtype=np.int8), np.ascontiguousarray(m.freepos, dtype=np.int8),
                           np.ascontiguousarray(m.pairs, dtype=np.int8), np.ascontiguousarray(m.type_of, dtype=np.int16),
                           np.ascontiguousarray(m.flag_of, dtype=np.int16), m.qoff, m.qdim, aux_table(m))
    return model._tab


def aux_table(model):
    """(4, 35) int64: rows 0..2 neighbourhood masks by intersection size 1..3; row 3: fast-path flag per block."""
    E = sorted(combinations(range(7), 4), key=lambda e: e[::-1])
    A = np.zeros((4, 35), dtype=np.int64)
    for i, e in enumerate(E):
        for k, f in enumerate(E):
            j = len(set(e) & set(f))
            if f != e and 1 <= j <= 3:
                A[j - 1, i] |= 1 << k
    for bi, b in enumerate(model.blocks):
        if (b.s, b.m) in ((1, 4), (2, 4), (3, 4)) and len(b.types) == 1 and b.flags[0].tolist() == [0, 1] \
                and b.nsig == 0 and b.nfree == 1:
            A[3, bi] = b.s
    return A


def qsum_flat(model, Qs, dtype=np.float64):
    """Flat Qsum vector from per-key Gram matrices Qs (list aligned with model.keys): Qsum = sum_h Q[pi_h, pi_h]."""
    out = np.zeros(model.qsize, dtype=dtype)
    for (bi, t), Q in zip(model.keys, Qs):
        P = model.blocks[bi].aut[t]
        n = P.shape[1]
        S = np.zeros((n, n), dtype=dtype)
        for g in range(P.shape[0]):
            S += Q[np.ix_(P[g], P[g])]
        off = int(model.qoff[bi, t])
        out[off:off + n * n] = S.ravel()
    return out
