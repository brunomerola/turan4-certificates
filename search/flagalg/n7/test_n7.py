"""Tests of the N = 7, r = 4 tables and kernels against the validated engine (../flags.py block_densities, a
direct configuration enumeration) and of the nauty canonical form.  Run inside tools/flagalg/n7:
    python -m pytest -q test_n7.py
"""
import random
from fractions import Fraction
from math import comb

import numpy as np
import pytest

from common import BLOCK_SM, E7, EIDX, N, admissible7, apply_perm, build_model
from kernels import coef_columns, moment_raw, qsum_flat, scan_exact, scan_float, tables, terms
from flags import block_densities, flag_family  # engine (parent directory, path set by common)
from hypergraphs import EveryPSetHasEdge


@pytest.fixture(scope="module", params=[5, 6, 7])
def model(request):
    return build_model(request.param)


def random_admissible(p, k, seed):
    rng = np.random.default_rng(seed)
    out = []
    while len(out) < k:
        dens = rng.uniform(0.05, 0.6)
        m = 0
        for i in range(35):
            if rng.random() < dens:
                m |= 1 << i
        if admissible7(np.array([m]), p)[0]:
            out.append(m)
    return np.array(out, dtype=np.int64)


def random_grams(model, seed):
    rng = np.random.default_rng(seed)
    Qs = []
    for bi, t in model.keys:
        n = int(model.qdim[bi, t])
        A = rng.normal(size=(n, n))
        Qs.append((A + A.T) / 2)
    return Qs


def engine_inner(model, Qs, H):
    """sum over keys of <Q, M_sigma(H)> via the engine's block_densities (direct configuration enumeration)."""
    pred = EveryPSetHasEdge(model.p)
    tot = np.zeros(H.size)
    for (bi, t), Q in zip(model.keys, Qs):
        b = model.blocks[bi]
        fam = flag_family(b.s, b.m, b.types[t], 4, pred)
        assert np.array_equal(fam.flags, b.flags[t])
        blk = block_densities(H, N, 4, fam)
        assert blk.denom == b.denom
        tot += np.bincount(blk.h, weights=blk.cnt * Q[blk.a, blk.b], minlength=H.size) / blk.denom
    return tot


def test_kernel_matches_engine(model):
    H = random_admissible(model.p, 12, seed=model.p)
    Qs = random_grams(model, seed=1)
    qs = qsum_flat(model, Qs)
    sl, bad = scan_float(H, *tables(model), qs, 0.0, 0.0, 0.0, 4)
    assert not bad.any()
    e = np.array([bin(int(h)).count("1") for h in H])
    kernel_inner = e / 35 - sl
    np.testing.assert_allclose(kernel_inner, engine_inner(model, Qs, H), rtol=1e-10, atol=1e-10)


def test_stationarity_term(model):
    """c_dd(H) = M_(1,4)(H)[edge, edge] (engine), and the float slack's tau terms."""
    H = random_admissible(model.p, 10, seed=7)
    zero = np.zeros(model.qsize)
    s0, _ = scan_float(H, *tables(model), zero, 0.0, 0.0, 0.3, 2)
    s1, _ = scan_float(H, *tables(model), zero, 1.0, 0.0, 0.3, 2)
    s2, _ = scan_float(H, *tables(model), zero, 0.0, 1.0, 0.3, 2)
    pred = EveryPSetHasEdge(model.p)
    fam = flag_family(1, 4, 0, 4, pred)
    blk = block_densities(H, N, 4, fam)
    edge = int(np.searchsorted(fam.flags, 1))
    cdd = np.zeros(H.size)
    sel = (blk.a == edge) & (blk.b == edge)
    np.add.at(cdd, blk.h[sel], blk.cnt[sel] / blk.denom)
    d = np.array([bin(int(h)).count("1") for h in H]) / 35
    np.testing.assert_allclose(s0, d)
    np.testing.assert_allclose(s0 - s1, 0.3 * d - cdd, atol=1e-14)
    np.testing.assert_allclose(s0 - s2, 0.3 * (1 - d) - (d - cdd), atol=1e-14)


def test_exact_kernel_matches_float(model):
    H = random_admissible(model.p, 30, seed=11)
    rng = np.random.default_rng(3)
    Qi = []
    for bi, t in model.keys:
        n = int(model.qdim[bi, t])
        A = rng.integers(-50, 50, size=(n, n))
        Qi.append(A + A.T)
    qi = qsum_flat(model, Qi, dtype=np.int64)
    MZ, AH, CN, nbad = scan_exact(H, *tables(model), qi, 3)
    assert nbad == 0 and CN.sum() == H.size
    sl, _ = scan_float(H, *tables(model), qi.astype(float), 0.0, 0.0, 0.0, 3)
    e = np.array([bin(int(h)).count("1") for h in H])
    Zf = np.rint((e / 35 - sl) * 5040).astype(np.int64)
    for i, h in enumerate(H):
        _, _, _, _, _, ee, cc = terms(H[i:i + 1], *tables(model)[:12], np.zeros((6, 6), dtype=np.int64),
                                      model.cdd_pos)
        assert MZ[ee[0], cc[0]] >= Zf[i]
    assert MZ.max() == Zf.max()


def test_terms_coefficients_and_moments(model):
    """LP-side helpers: per-key terms give the same <v v^T, M(H)> as the scan, and v^T Y v = sum_H y_H coef."""
    H = random_admissible(model.p, 8, seed=5)
    keyid = np.full((len(model.blocks), 6), -1, dtype=np.int64)
    for k, (bi, t) in enumerate(model.keys):
        keyid[bi, t] = k
    rows, keys, aa, bb, ww, ee, cc = terms(H, *tables(model)[:12], keyid, model.cdd_pos)
    rng = np.random.default_rng(2)
    y = rng.random(H.size)
    for k, (bi, t) in enumerate(model.keys):
        n = int(model.qdim[bi, t])
        P = model.blocks[bi].aut[t]
        V = rng.normal(size=(2, n))
        sel = keys == k
        C = coef_columns(rows[sel], aa[sel], bb[sel], ww[sel], P, V, H.size)
        Qs = [np.zeros((int(model.qdim[b2, t2]),) * 2) for b2, t2 in model.keys]
        Qs[k] = np.outer(V[0], V[0])
        sl, _ = scan_float(H, *tables(model), qsum_flat(model, Qs), 0.0, 0.0, 0.0, 2)
        np.testing.assert_allclose(C[0], ee / 35 - sl, rtol=1e-9, atol=1e-12)
        Rm = moment_raw(rows[sel], aa[sel], bb[sel], ww[sel], y, n)
        Y = sum(Rm[np.ix_(P[g], P[g])] for g in range(P.shape[0]))
        np.testing.assert_allclose(V[1] @ Y @ V[1], C[1] @ y, rtol=1e-9, atol=1e-12)


def test_canon7_invariant():
    from enum7 import canon7
    rnd = random.Random(4)
    for _ in range(500):
        m = rnd.getrandbits(35) & (rnd.getrandbits(35) if rnd.random() < 0.5 else (1 << 35) - 1)
        g = list(range(7))
        rnd.shuffle(g)
        assert canon7(m) == canon7(apply_perm(m, 7, g))
        # canonical form is a relabelling of m
        assert bin(canon7(m)).count("1") == bin(m).count("1")


def test_block_structure():
    M = build_model(6)
    assert [(b.s, b.m) for b in M.blocks] == list(BLOCK_SM)
    assert M.gram_dims() == [2, 2, 2, 24, 16, 16, 1023, 1024, 1024, 1024, 1024, 1024]


def test_total_mass(model):
    """With Q = all-ones, <Q, M_sigma(H)> = P(H[theta] = sigma as labelled graphs) for a random injective s-tuple
    theta; summed over the types of a block this is checked against a brute-force count over ordered s-tuples."""
    from itertools import permutations
    from math import factorial
    H = random_admissible(model.p, 6, seed=9)
    ones = [np.ones((int(model.qdim[b, t]),) * 2) for b, t in model.keys]
    sl, bad = scan_float(H, *tables(model), qsum_flat(model, ones), 0.0, 0.0, 0.0, 2)
    e = np.array([bin(int(h)).count("1") for h in H])
    expect = np.zeros(H.size)
    for i, h in enumerate(H):
        h = int(h)
        for b in model.blocks:
            types = set(b.types)
            cnt = 0
            for th in permutations(range(7), b.s):
                sig = 0
                for j, f in enumerate(E7[:comb(b.s, 4)]):
                    if (h >> EIDX[tuple(sorted(th[x] for x in f))]) & 1:
                        sig |= 1 << j
                cnt += sig in types
            expect[i] += cnt / (factorial(7) // factorial(7 - b.s))
    np.testing.assert_allclose(e / 35 - sl, expect, rtol=1e-12)


def test_fixed_point_bound_exact(model):
    """certify7.bound_from_groups (group maxima of Z) equals min_H a_H / (1 + beta_H) computed per graph with
    Fractions, for integer Grams and nonzero tau (the stationarity path of the exact certificate)."""
    from certify7 import bound_from_groups
    H = random_admissible(model.p, 40, seed=21)
    rng = np.random.default_rng(8)
    Qi = []
    for bi, t in model.keys:
        n = int(model.qdim[bi, t])
        A = rng.integers(-3, 4, size=(3, n))
        Qi.append(A.T @ A)
    qi = qsum_flat(model, Qi, dtype=np.int64)
    M2, t1, t2 = 1 << 6, 5, 3
    MZ, AH, CN, nbad = scan_exact(H, *tables(model), qi, 4)
    b, _ = bound_from_groups(MZ, CN, M2, t1, t2)
    _, _, _, _, _, ee, cc = terms(H, *tables(model)[:12], np.zeros((6, 6), dtype=np.int64), model.cdd_pos)
    sl, _ = scan_float(H, *tables(model), qi.astype(float), 0.0, 0.0, 0.0, 2)
    best = None
    for i in range(H.size):
        Z = int(round((ee[i] / 35 - sl[i]) * 5040))
        d, c = Fraction(int(ee[i]), 35), Fraction(int(cc[i]), 70)
        tau1, tau2 = Fraction(t1, M2), Fraction(t2, M2)
        aH = d - Fraction(Z, 5040 * M2) + tau1 * c + tau2 * (d - c)
        val = aH / (1 + tau1 * d + tau2 * (1 - d))
        best = val if best is None or val < best else best
    assert b == best
