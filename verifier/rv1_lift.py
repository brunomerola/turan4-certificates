"""Independent checker, CLAIM 1(a): exact check of the lifting identity used in the scope argument (own code).

For a type size s let mt(s) = max{m : 2m - s <= 7}.  For flag sizes s <= ma, mb <= mt(s) (any pair, including mixed
sizes ma != mb and the trivial m = s), the mixed moment matrix of a 7-vertex graph H,
    Mmix[F, F'] = P_{(theta, U1, U2)}[ H[theta] = tau, (H^{U1}_theta) ~ F, (H^{U2}_theta) ~ F' ],
    U1 an (ma-s)-set, U2 an (mb-s)-set, disjoint, (theta, U1, U2) uniform,
must equal  sum_{G, G'} p(F; G) p(F'; G') Mtop[G, G']  with Mtop the moment matrix of the block (s, mt(s)) and
p(F; G) the probability that a uniform (ma-s)-subset W of the non-roots of G gives (G[roots + W]) ~ F.
Checked exactly (Fractions) for every labelled type, on every support graph of the dual point and on random
7-vertex graphs.  If it holds, every (mixed or smaller) block's moment matrix is T Y_top T^T, PSD when Y_top is.
Usage: python rv1_lift.py DUAL.json NRANDOM
"""
import json
import random
import sys
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial

import rv1_lib as L

N = 7


def mixed_counts(h, s, ma, mb):
    out = {}
    tot = 0
    for th in permutations(range(N), s):
        t = L.code(h, th, s)
        rest = tuple(v for v in range(N) if v not in th)
        A = [(set(U), L.canon(h, th, U, ma)) for U in combinations(rest, ma - s)]
        B = [(set(U), L.canon(h, th, U, mb)) for U in combinations(rest, mb - s)]
        d = out.setdefault(t, {})
        for Ua, fa in A:
            for Ub, fb in B:
                if Ua.isdisjoint(Ub):
                    d[(fa, fb)] = d.get((fa, fb), 0) + 1
                    tot += 1
    cf = factorial(N) // factorial(N - s) * comb(N - s, ma - s) * comb(N - ma, mb - s)
    assert tot == cf
    return out, cf


_pcache = {}


def p_restrict(G, s, m, ma):
    """G: canonical code of a flag on m vertices; returns {F: Fraction} for restrictions to ma vertices."""
    key = (G, s, m, ma)
    if key in _pcache:
        return _pcache[key]
    hG = 0
    for j, e in enumerate(L.SUB[m]):
        if (G >> j) & 1:
            hG |= 1 << L.IX[((e[0] * N + e[1]) * N + e[2]) * N + e[3]]
    roots = tuple(range(s))
    Ws = list(combinations(range(s, m), ma - s))
    out = {}
    for W in Ws:
        F = L.canon(hG, roots, W, ma)
        out[F] = out.get(F, 0) + Fraction(1, len(Ws))
    _pcache[key] = out
    return out


def check_graph(h):
    nchk = 0
    for s in range(0, 6):
        mt = (N + s) // 2
        top, cft = mixed_counts(h, s, mt, mt)
        for ma in range(s, mt + 1):
            for mb in range(s, mt + 1):
                if (ma, mb) == (mt, mt):
                    continue
                mix, cfm = mixed_counts(h, s, ma, mb)
                assert set(mix) <= set(top) | set()
                for tau in set(top) | set(mix):
                    lifted = {}
                    for (G, G2), c in top.get(tau, {}).items():
                        pa = p_restrict(G, s, mt, ma)
                        pb = p_restrict(G2, s, mt, mb)
                        for F, x in pa.items():
                            for F2, y in pb.items():
                                lifted[(F, F2)] = lifted.get((F, F2), 0) + Fraction(c, cft) * x * y
                    direct = {k: Fraction(c, cfm) for k, c in mix.get(tau, {}).items()}
                    lifted = {k: v for k, v in lifted.items() if v}
                    if lifted != direct:
                        return False, (s, ma, mb, tau)
                    nchk += 1
    return True, nchk


def main():
    D = json.load(open(sys.argv[1]))
    nr = int(sys.argv[2])
    graphs = [L.colex_to_own(int(x["mask"])) for x in D["support"]]
    rng = random.Random(20261006)
    graphs += [rng.getrandbits(35) for _ in range(nr)]
    ok_all = True
    tot = 0
    for i, h in enumerate(graphs):
        ok, info = check_graph(h)
        if not ok:
            print("MISMATCH", i, h, info, flush=True)
            ok_all = False
        else:
            tot += info
    print(f"graphs checked {len(graphs)} (support {len(D['support'])} + random {nr}); (graph, s, ma, mb, type) "
          f"identities checked {tot}; all equal: {ok_all}")
    print("pairs (s, mt(s)):", [(s, (N + s) // 2) for s in range(6)], "; s = 6, 7 only allow m = s (1x1)")


if __name__ == "__main__":
    main()
