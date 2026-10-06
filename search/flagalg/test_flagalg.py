"""Unit tests for the flag-algebra engine.  Run:  python -m pytest -q test_flagalg.py  (inside tools/flagalg)."""
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial

import numpy as np
import pytest

from exact import certify_eig, certify_sharp, is_psd_exact
from flags import build_problem, flag_family
from hypergraphs import EveryPSetHasEdge, NoPredicate, PSetEdgeCounts, colex_index, edges, enumerate_classes
from sdp import solve_flag_sdp


# ------------------------------------------------------------------------------------------ enumeration counts
def test_counts_graphs():                      # OEIS A000088
    R = enumerate_classes(7, 2, NoPredicate())
    assert [R[n].size for n in range(8)] == [1, 1, 2, 4, 11, 34, 156, 1044]


def test_counts_3graphs():                     # OEIS A000665
    R = enumerate_classes(6, 3, NoPredicate())
    assert [R[n].size for n in range(7)] == [1, 1, 1, 2, 5, 34, 2136]


def test_counts_admissible():
    # complements of K_4^(3)-free 3-graphs on 6 vertices: 964 (Baber-Talbot 2012)
    assert enumerate_classes(6, 3, EveryPSetHasEdge(4))[6].size == 964
    # {K_4, induced E_1}-free on 6 vertices: 34 (Baber-Talbot 2012, Razborov 2010)
    assert enumerate_classes(6, 3, PSetEdgeCounts(4, [1, 2, 4]))[6].size == 34
    # triangle-free graphs on 7 vertices (OEIS A006785): 107
    assert enumerate_classes(7, 2, EveryPSetHasEdge(3))[7].size == 107
    # 4-graphs on 6 vertices: every 5-set has an edge (dual graph without isolated vertices) 122; nonempty 155
    assert enumerate_classes(6, 4, EveryPSetHasEdge(5))[6].size == 122
    assert enumerate_classes(6, 4, EveryPSetHasEdge(6))[6].size == 155


def test_flag_count_s2_m4_r3():
    # 4 triples on {0,1,2,3}; the swap 2<->3 fixes 023 and 123 and exchanges 012, 013: 12 orbits, minus empty = 11
    assert flag_family(2, 4, 0, 3, EveryPSetHasEdge(4)).n == 11


# ------------------------------------------------------------------------------------------ densities
def naive_counts(Hmask, N, r, fam):
    """Independent pure-Python count of (theta, U1, U2) configurations per flag pair."""
    E = edges(N, r)
    Hset = {frozenset(E[i]) for i in range(len(E)) if (Hmask >> i) & 1}
    s, m = fam.s, fam.m
    pos = list(combinations(range(m), r))

    def labelled_canon(order):
        best = None
        for pu in permutations(order[s:]):
            o = tuple(order[:s]) + pu
            mask = sum(1 << colex_index(e) for e in pos if frozenset(o[i] for i in e) in Hset)
            best = mask if best is None or mask < best else best
        return best

    index = {int(f): i for i, f in enumerate(fam.flags)}
    cnt = {}
    for th in permutations(range(N), s):
        sig = sum(1 << colex_index(e) for e in combinations(range(s), r) if frozenset(th[i] for i in e) in Hset)
        if sig != fam.sigma:
            continue
        rest = [v for v in range(N) if v not in th]
        for U1 in combinations(rest, m - s):
            rest2 = [v for v in rest if v not in U1]
            for U2 in combinations(rest2, m - s):
                a = index[labelled_canon(th + U1)]
                b = index[labelled_canon(th + U2)]
                cnt[(a, b)] = cnt.get((a, b), 0) + 1
    return cnt


@pytest.mark.parametrize("N,r,p", [(5, 2, 3), (5, 3, 4), (6, 4, 6)])
def test_densities_against_naive(N, r, p):
    prob = build_problem(N, r, EveryPSetHasEdge(p))
    rng = np.random.default_rng(1)
    hs = rng.choice(prob.nH, size=min(6, prob.nH), replace=False)
    for blk in prob.blocks:
        for h in hs:
            sel = blk.h == h
            eng = {(int(a), int(b)): int(c) for a, b, c in zip(blk.a[sel], blk.b[sel], blk.cnt[sel])}
            assert eng == naive_counts(int(prob.H[h]), N, r, blk.fam)


def test_density_identities():
    prob = build_problem(6, 3, EveryPSetHasEdge(4))
    for blk in prob.blocks:
        s, m, N = blk.fam.s, blk.fam.m, prob.N
        # symmetric counts
        d = {(int(h), int(a), int(b)): int(c) for h, a, b, c in zip(blk.h, blk.a, blk.b, blk.cnt)}
        assert all(d.get((h, b, a)) == c for (h, a, b), c in d.items())
        assert blk.denom == factorial(N) // factorial(N - s) * comb(N - s, m - s) * comb(N - m, m - s)
        tot = np.bincount(blk.h, weights=blk.cnt, minlength=prob.nH)
        assert (tot <= blk.denom).all()


# ------------------------------------------------------------------------------------------ SDP conventions
def test_clarabel_psd_triangle_convention():
    import clarabel
    import scipy.sparse as sp
    # x = svec(X) for 3x3; fix X11 = 1, X22 = 1, X33 = 4; minimise x[3]: if x[3] = sqrt2 X13 the optimum is -2 sqrt2
    q = np.zeros(6)
    q[3] = 1.0
    Aeq = sp.lil_matrix((3, 6))
    Aeq[0, 0] = 1
    Aeq[1, 2] = 1
    Aeq[2, 5] = 1
    A = sp.vstack([Aeq, -sp.eye(6)]).tocsc()
    b = np.concatenate([[1, 1, 4], np.zeros(6)])
    st = clarabel.DefaultSettings()
    st.verbose = False
    sol = clarabel.DefaultSolver(sp.csc_matrix((6, 6)), q, A, b, [clarabel.ZeroConeT(3), clarabel.PSDTriangleConeT(3)],
                                 st).solve()
    assert abs(sol.obj_val + 2 * np.sqrt(2)) < 1e-6


def test_is_psd_exact():
    assert is_psd_exact([[1, 1], [1, 1]])[:2] == (True, 1)
    assert is_psd_exact([[1, 2], [2, 1]])[0] is False
    assert is_psd_exact([[0, 1], [1, 0]])[0] is False
    assert is_psd_exact([[Fraction(1, 3), 0], [0, 0]])[:2] == (True, 1)


# ------------------------------------------------------------------------------------------ known bounds
@pytest.mark.parametrize("N,r,p,target", [(3, 2, 3, Fraction(1, 2)), (4, 2, 3, Fraction(1, 2)),
                                          (5, 2, 3, Fraction(1, 2)), (4, 2, 4, Fraction(1, 3)),
                                          (5, 2, 4, Fraction(1, 3))])
def test_sharp_turan_r2(N, r, p, target):
    prob = build_problem(N, r, EveryPSetHasEdge(p))
    num, den = prob.density_num_den()
    sol = solve_flag_sdp(prob, num / den)
    obj = [Fraction(int(x), den) for x in num]
    ce = certify_eig(prob, obj, sol["Q"])
    assert all(ok for ok, _ in ce["psd"])
    assert target - Fraction(1, 10 ** 8) < ce["b_exact"] <= target
    slf = np.array([float(x) for x in ce["slacks"]]) - float(ce["b_exact"])
    cs = certify_sharp(prob, obj, sol["Q"], slf, target)
    assert cs["success"] and cs["b_exact"] == target


def test_no_sos_equals_min_density():
    prob = build_problem(6, 3, EveryPSetHasEdge(4), type_mode="none")
    num, den = prob.density_num_den()
    sol = solve_flag_sdp(prob, num / den)
    assert abs(sol["b"] - 0.3) < 1e-8          # T(6,4,3) = 6 triples, 6/20


def test_k4_minus_e1_sharp():
    """Razborov 2010: pi_ind(K_4^(3), E_1) = 5/9, i.e. min density 4/9 in complement form, from N = 6."""
    prob = build_problem(6, 3, PSetEdgeCounts(4, [1, 2, 4]))
    num, den = prob.density_num_den()
    sol = solve_flag_sdp(prob, num / den)
    obj = [Fraction(int(x), den) for x in num]
    ce = certify_eig(prob, obj, sol["Q"])
    slf = np.array([float(x) for x in ce["slacks"]]) - float(ce["b_exact"])
    cs = certify_sharp(prob, obj, sol["Q"], slf, Fraction(4, 9))
    assert cs["success"] and cs["b_exact"] == Fraction(4, 9)


def test_lottery_r2_never_exceeds_truth():
    """l(k,2,p) = 1/((p-1) C(k,2)) (Fueredi-Szekely-Zubor): the local packing bound must not exceed it."""
    from lottery import run
    for (k, p, N, m, truth) in [(3, 3, 5, 4, Fraction(1, 6)), (4, 3, 5, 4, Fraction(1, 12))]:
        for tau, sym in ((False, False), (True, False), (True, True)):
            res = run(k, 2, p, N, m, verbose=False, tau_sos=tau, tau_sym=sym)
            b = Fraction(res["eig_rounding"]["b_exact"])
            assert b <= truth
            assert float(b) > float(truth) - 1e-6


def test_lottery_m4_equals_bate():
    """(4,3,4) with m = 4 (no tau-SOS): every rooted 4-flag extends symmetrically to a K_4 (link empty, perfect
    matching, C_4 or K_4), so the local packing bound equals Bate's bound at the same N."""
    from lottery import run
    res = run(4, 3, 4, 5, 4, verbose=False)
    assert abs(res["sdp"]["b"] - res["turan_same_N"]["bate"]) < 1e-8
