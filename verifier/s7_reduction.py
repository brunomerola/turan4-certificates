"""Independent checker: exact checks of the reductions used by the certificates (own code, s7_lib only).

r = 3 (cosig3):  for random 3-graphs G on n = 5..9 vertices
  (a) co2(G) = sum_T c_T (c_T - 1) + 3 e(G);
  (b) BCL eq. (1):  co2/(C(n,2)(n-2)^2) = sigma_xy (n-3)/(n-2) + 3e/(C(n,2)(n-2)^2)  <=  sigma_xy + 1/(n-2);
  (c) sigma_xy = E_Q[C(k_Q,2)]/6  and  cosig3 (explicit (T,x,y) loop) = 1 - sigma_xy;
  (d) 7-vertex numerator: 210 * cosig3 = 210 - sum_Q C(k_Q, 2)  (formula used by s7_eval);
  (e) exact averaging: cosig3_n(G) = mean over 7-subsets S of cosig3_7(G[S])   (n = 8, 9).
r = 4 (sigma4, complement form):  for random 4-graphs G, H = complement
  (f) co2(G)/(C(n,3)(n-3)^2) = gamma(G) - (gamma(G) - d(G))/(n-3)  and gamma(G) = 1 - kappa(H);
  (g) 420 * kappa_7(H) = 24 e(H) - sum_T c_T(c_T - 1);  (h) exact averaging over 7-subsets (n = 8, 9);
  (i) G K_5^{4-}-free  <=>  every 5-set of H spans >= 2 edges  <=>  s7_lib.admissible('c4_K5m', H).
Predicates (copy lists of s7_lib) against independent definitional tests on random graphs, n = 5, 6, 7:
  J4: some vertex whose link (graph on the other vertices) contains K_4;   K4m: some 4-set with >= 3 triples;
  C5: some 5 distinct vertices v0..v4 with all five cyclic triples {v_i, v_i+1, v_i+2} edges;
  K5lt: some 5-set whose non-edges are <= 1, or exactly 2 sharing two vertices.
"""
from __future__ import annotations

import json
import random
import sys
from fractions import Fraction
from itertools import combinations, permutations
from math import comb

import s7_lib as L

rng = random.Random(20261006)
res = {}


def induced(mask, n, r, S):
    """mask of G[S] relabelled 0..|S|-1 in increasing order."""
    E = L.edge_set(mask, n, r)
    return sum(1 << L.cidx(tuple(S.index(v) for v in e)) for e in combinations(S, r) if frozenset(e) in E)


# ---------------------------------------------------------------- r = 3
cnt = 0
for n in range(5, 10):
    for trial in range(12):
        p = rng.choice([0.2, 0.5, 0.8, 0.95])
        g = L.rand_mask(n, 3, p, rng)
        c = L.codeg(g, n, 3)
        e = bin(g).count("1")
        co = L.co2(g, n, 3)
        assert co == sum(x * (x - 1) for x in c.values()) + 3 * e
        sxy = L.sigma_xy(g, n, 3)
        s1 = L.sigma_eq1(g, n, 3)
        assert s1 == sxy * Fraction(n - 3, n - 2) + Fraction(3 * e, comb(n, 2) * (n - 2) ** 2)
        assert s1 <= sxy + Fraction(1, n - 2)
        E = L.edge_set(g, n, 3)
        q = Fraction(sum(comb(sum(frozenset(t) in E for t in combinations(Q, 3)), 2) for Q in combinations(range(n), 4)),
                     6 * comb(n, 4))
        assert q == sxy and L.cosig3_value(g, n) == 1 - sxy
        cnt += 1
res["r3_identities_graphs"] = cnt
cnt = 0
for trial in range(40):
    g = L.rand_mask(7, 3, rng.random(), rng)
    E = L.edge_set(g, 7, 3)
    F = 210 - sum(comb(sum(frozenset(t) in E for t in combinations(Q, 3)), 2) for Q in combinations(range(7), 4))
    assert Fraction(F, 210) == L.cosig3_value(g, 7)
    cnt += 1
res["r3_num7_graphs"] = cnt
for n in (8, 9):
    for trial in range(3):
        g = L.rand_mask(n, 3, rng.choice([0.3, 0.6, 0.9]), rng)
        subs = list(combinations(range(n), 7))
        avg = sum(L.cosig3_value(induced(g, n, 3, list(S)), 7) for S in subs) / len(subs)
        assert avg == L.cosig3_value(g, n)
res["r3_averaging"] = "exact for n = 8, 9 (6 graphs)"

# ---------------------------------------------------------------- r = 4
cnt = 0
for n in range(5, 9):
    for trial in range(6):
        g = L.rand_mask(n, 4, rng.choice([0.3, 0.6, 0.9]), rng)
        h = ((1 << comb(n, 4)) - 1) ^ g
        E = L.edge_set(g, n, 4)
        d = Fraction(len(E), comb(n, 4))
        good = tot = 0
        for T in combinations(range(n), 3):
            out = [v for v in range(n) if v not in T]
            for x in out:
                for y in out:
                    if x != y:
                        tot += 1
                        good += (frozenset(T + (x,)) in E) and (frozenset(T + (y,)) in E)
        gam = Fraction(good, tot)
        assert Fraction(L.co2(g, n, 4), comb(n, 3) * (n - 3) ** 2) == gam - (gam - d) / (n - 3)
        assert gam == 1 - L.sigma4_value(h, n)
        cnt += 1
res["r4_identities_graphs"] = cnt
for trial in range(30):
    h = L.rand_mask(7, 4, rng.random(), rng)
    c = L.codeg(h, 7, 4)
    assert Fraction(24 * bin(h).count("1") - sum(x * (x - 1) for x in c.values()), 420) == L.sigma4_value(h, 7)
for n in (8,):
    for trial in range(3):
        h = L.rand_mask(n, 4, rng.choice([0.3, 0.6, 0.9]), rng)
        subs = list(combinations(range(n), 7))
        avg = sum(L.sigma4_value(induced(h, n, 4, list(S)), 7) for S in subs) / len(subs)
        assert avg == L.sigma4_value(h, n)
res["r4_averaging"] = "exact for n = 8 (3 graphs)"
cnt = 0
for n in (5, 6, 7):
    for trial in range(300):
        h = L.rand_mask(n, 4, rng.choice([0.2, 0.5, 0.8, 0.95]), rng)
        g = ((1 << comb(n, 4)) - 1) ^ h
        Eg = L.edge_set(g, n, 4)
        Eh = L.edge_set(h, n, 4)
        k5m_in_g = any(sum(frozenset(q) in Eg for q in combinations(P, 4)) >= 4 for P in combinations(range(n), 5))
        every5 = all(sum(frozenset(q) in Eh for q in combinations(P, 4)) >= 2 for P in combinations(range(n), 5))
        assert (not k5m_in_g) == every5 == L.admissible("c4_K5m", h, n)
        cnt += 1
res["r4_K5m_predicate_graphs"] = cnt


# ---------------------------------------------------------------- 3-graph predicates
def has_J4(E, n):
    for v in range(n):
        O = [u for u in range(n) if u != v]
        for W in combinations(O, 4):
            if all(frozenset((v, a, b)) in E for a, b in combinations(W, 2)):
                return True
    return False


def has_K4m(E, n):
    return any(sum(frozenset(t) in E for t in combinations(Q, 3)) >= 3 for Q in combinations(range(n), 4))


def has_C5(E, n):
    for P in permutations(range(n), 5):
        if all(frozenset((P[i], P[(i + 1) % 5], P[(i + 2) % 5])) in E for i in range(5)):
            return True
    return False


def has_K5lt(E, n):
    for P in combinations(range(n), 5):
        non = [frozenset(t) for t in combinations(P, 3) if frozenset(t) not in E]
        if len(non) <= 1 or (len(non) == 2 and len(non[0] & non[1]) == 2):
            return True
    return False


stats = {}
for prob, tests in (("s3_J4", [has_J4]), ("s3_K4m", [has_K4m]), ("s3_K4mC5", [has_K4m, has_C5]),
                    ("s3_K5lt", [has_K5lt])):
    nok = nbad = 0
    for n in (5, 6, 7):
        for trial in range(250 if n < 7 else 120):
            p = rng.choice([0.1, 0.3, 0.5, 0.7, 0.85, 0.95])
            g = L.rand_mask(n, 3, p, rng)
            E = L.edge_set(g, n, 3)
            free = not any(t(E, n) for t in tests)
            assert free == L.admissible(prob, g, n), (prob, n, g)
            nok += free
            nbad += not free
    stats[prob] = {"admissible": nok, "inadmissible": nbad, "copies7": len(L.problem_copies(prob, 7))}
res["predicates_3graph"] = stats
# the C5 branch of {K4^-, C5}: plant a tight 5-cycle on random 5 vertices of a sparse random graph
nc5 = nboth = 0
for trial in range(400):
    n = rng.choice([5, 6, 7])
    g = L.rand_mask(n, 3, rng.choice([0.0, 0.05, 0.1, 0.2]), rng)
    P = rng.sample(range(n), 5)
    for i in range(5):
        g |= 1 << L.cidx((P[i], P[(i + 1) % 5], P[(i + 2) % 5]))
    E = L.edge_set(g, n, 3)
    if not has_K4m(E, n):
        nc5 += 1
        assert has_C5(E, n) and not L.admissible("s3_K4mC5", g, n) and L.admissible("s3_K4m", g, n)
    # random sparse graphs: compare C5 test with the copy list restricted to C5
    h = L.rand_mask(n, 3, rng.choice([0.1, 0.2, 0.3]), rng)
    Eh = L.edge_set(h, n, 3)
    c5copies = L.copies(L.parse(L.FORB["C5"]), n, 3)
    assert has_C5(Eh, n) == any((h & c) == c for c in c5copies)
    nboth += has_C5(Eh, n)
res["C5_branch"] = {"planted_K4m_free_with_C5": nc5, "random_with_C5": nboth}
res["ALL_OK"] = True
print(json.dumps(res, indent=1))
json.dump(res, open(sys.argv[1], "w"), indent=1)
