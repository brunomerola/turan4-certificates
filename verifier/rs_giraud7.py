"""Independent checker: Giraud's limit on 7 vertices in complement form (own code).

Limit object: each vertex is a row (prob 1/2) or a column; row-column entries i.i.d. fair bits.  G (K_5^4-free) has
the 3+1 sets and the 2+2 sets with ODD 2x2 sum; the complement H (admissible) has the 4+0, 0+4 sets and the 2+2 sets
with EVEN sum.  The exact distribution phi of H on 7 labelled vertices is enumerated (all side patterns, all entries)
and mapped to isomorphism classes by membership in the orbits (all 5040 relabellings) of the producer's ten listed
Giraud class masks.  Checks: support == those ten classes, probabilities, E_phi f == 33/64, every graph admissible.
Since phi is the 7-vertex law of a K_5^4-free construction with s = 31/64, any valid plain bound b satisfies
b <= E_phi[f - c_Q] <= E_phi f = 33/64 (c_Q has nonnegative expectation in the limit), so 33/64 is the best possible.
"""
import json
from fractions import Fraction
from itertools import combinations, permutations, product

from rs_lib import colex_rank, f_value

GIRAUD = {2290237506: "105/2048", 17418883036: "35/4096", 19736566032: "105/1024", 26383237633: "315/1024",
          26894418433: "105/512", 29337331976: "105/1024", 31813832752: "315/4096", 31837153552: "7/64",
          34326183967: "21/1024", 34359738367: "1/64"}
E = list(combinations(range(7), 4))
CR = [colex_rank(e) for e in E]


def orbit(h):
    es = [e for i, e in enumerate(E) if (h >> CR[i]) & 1]
    return {sum(1 << colex_rank(tuple(sorted(g[v] for v in e))) for e in es) for g in permutations(range(7))}


def main():
    where = {}
    for g in GIRAUD:
        for x in orbit(g):
            assert x not in where
            where[x] = g
    phi = {}
    P5 = [sum(1 << colex_rank(f) for f in combinations(P, 4)) for P in combinations(range(7), 5)]
    for sd in product((0, 1), repeat=7):
        prs = [(i, j) for i in range(7) for j in range(7) if sd[i] == 1 and sd[j] == 0]
        w = Fraction(1, 128 * (1 << len(prs)))
        idx = {p: t for t, p in enumerate(prs)}
        for b in range(1 << len(prs)):
            h = 0
            for i, e in enumerate(E):
                rows = [v for v in e if sd[v]]
                if len(rows) in (0, 4):
                    h |= 1 << CR[i]
                elif len(rows) == 2:
                    cols = [v for v in e if not sd[v]]
                    par = sum((b >> idx[(r, c)]) & 1 for r in rows for c in cols) & 1
                    if par == 0:
                        h |= 1 << CR[i]
            assert all(h & m for m in P5), "Giraud complement not admissible?"
            phi[h] = phi.get(h, 0) + w
    cls = {}
    unknown = 0
    for h, w in phi.items():
        g = where.get(h)
        if g is None:
            unknown += 1
            continue
        cls[g] = cls.get(g, 0) + w
    ok = unknown == 0 and sum(cls.values()) == 1 and set(cls) == set(GIRAUD)
    ok &= all(cls[g] == Fraction(GIRAUD[g]) for g in GIRAUD)
    Ef = Fraction(0)
    fv = {}
    for g in GIRAUD:
        es = {frozenset(e) for i, e in enumerate(E) if (g >> CR[i]) & 1}
        fv[g] = f_value(es, 7)
        Ef += cls.get(g, 0) * fv[g]
    ok &= Ef == Fraction(33, 64)
    out = {"labelled_graphs_in_support": len(phi), "outside_the_ten_orbits": unknown,
           "class_probabilities": {str(g): str(cls.get(g)) for g in GIRAUD},
           "f_values": {str(g): str(fv[g]) for g in GIRAUD}, "E_phi_f": str(Ef), "PASS": ok}
    json.dump(out, open("results_giraud7.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
