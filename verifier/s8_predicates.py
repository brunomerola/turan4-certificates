"""Independent checker: the copy-list predicates of s8_lib against independent definitional tests, new families.

3-graphs (sub mode, G admissible iff no copy of F in E(G)):
  K_5^=  : some 5-set whose non-edges are <= 1, or exactly 2 sharing exactly ONE vertex;
  K_5^{3-}: some 5-set spanning >= 9 triples;      K_6^3: some 6-set spanning all 20 triples;
  C_5    : some 5 distinct vertices v0..v4 with all five cyclic triples {v_i, v_i+1, v_i+2} edges.
4-graphs (complement form, H admissible iff the complement G is F-free):
  K_6^{4-}: every 6-set of H spans >= 2 edges;   K_6^4: every 6-set of H spans >= 1 edge;   K_7^4: H has an edge
  (on n >= 7 vertices: every 7-set spans >= 1 edge).
Random graphs with n = 5..7 at several densities, plus planted copies.
"""
from __future__ import annotations

import json
import random
import sys
from itertools import combinations, permutations
from math import comb

import s8_lib as L

rng = random.Random(20261008)


def es(g, n, r):
    return L.edge_set(g, n, r)


def has_K5eq(E, n):
    for P in combinations(range(n), 5):
        non = [frozenset(t) for t in combinations(P, 3) if frozenset(t) not in E]
        if len(non) <= 1 or (len(non) == 2 and len(non[0] & non[1]) == 1):
            return True
    return False


def has_K5m3(E, n):
    return any(sum(frozenset(t) in E for t in combinations(P, 3)) >= 9 for P in combinations(range(n), 5))


def has_K63(E, n):
    return any(all(frozenset(t) in E for t in combinations(P, 3)) for P in combinations(range(n), 6))


def has_C5(E, n):
    for P in permutations(range(n), 5):
        if all(frozenset((P[i], P[(i + 1) % 5], P[(i + 2) % 5])) in E for i in range(5)):
            return True
    return False


def every_ps(E, n, p, lam):
    return all(sum(frozenset(q) in E for q in combinations(P, 4)) >= lam for P in combinations(range(n), p))


res = {}
tests3 = {"s3_K5eq": has_K5eq, "s3_K5m": has_K5m3, "s3_K6": has_K63, "s3_C5": has_C5}
for prob, t in tests3.items():
    ok = bad = 0
    for n in (5, 6, 7):
        for trial in range(300 if n < 7 else 150):
            p = rng.choice([0.2, 0.5, 0.7, 0.85, 0.95, 0.99])
            g = L.rand_mask(n, 3, p, rng)
            if trial % 5 == 0 and n >= 5:          # plant a copy of the forbidden graph, then delete a random triple
                F = L.parse(L.FORB[L.PROBLEMS[prob][2][0]])
                vs = sorted({v for e in F for v in e})
                if len(vs) <= n:
                    img = rng.sample(range(n), len(vs))
                    mp = dict(zip(vs, img))
                    for e in F:
                        g |= 1 << L.cidx(tuple(mp[v] for v in e))
                    if rng.random() < 0.5:
                        g &= ~(1 << rng.randrange(comb(n, 3)))
            E = es(g, n, 3)
            free = not t(E, n)
            assert free == L.admissible(prob, g, n), (prob, n, g)
            ok += free
            bad += not free
    res[prob] = {"admissible": ok, "inadmissible": bad, "copies7": len(L.problem_copies(prob, 7))}
tests4 = {"c4_K6m": (6, 2), "c4_K6": (6, 1), "c4_K7": (7, 1)}
for prob, (p, lam) in tests4.items():
    ok = bad = 0
    for n in (6, 7):
        for trial in range(400):
            q = rng.choice([0.0, 0.01, 0.03, 0.06, 0.1, 0.3])
            h = L.rand_mask(n, 4, q, rng)
            E = es(h, n, 4)
            adm = every_ps(E, n, p, lam) if n >= p else True
            assert adm == L.admissible(prob, h, n), (prob, n, h)
            ok += adm
            bad += not adm
    res[prob] = {"admissible": ok, "inadmissible": bad, "copies7": len(L.problem_copies(prob, 7))}
res["ALL_OK"] = True
print(json.dumps(res, indent=1))
json.dump(res, open(sys.argv[1], "w"), indent=1)
