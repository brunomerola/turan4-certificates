"""Independent checker (own code): basic combinatorics for 4-graphs, written from the definitions.

Encodings used by the reviewer (NOT the producer's):
  * "own" 7-vertex mask: bit i <-> i-th 4-subset of {0..6} in LEXICOGRAPHIC order (itertools.combinations order).
  * producer masks are COLEX (bit of {a<b<c<d} = C(a,1)+C(b,2)+C(c,3)+C(d,4)); decoded here with an own colex rank.
  * flag masks on [m] with roots 0..s-1: own "flag order" = the root-only 4-sets first (lex), then the others (lex),
    so the low C(s,4) bits are the labelled type.
No producer module is imported by any of the rs_* checkers.
"""
from fractions import Fraction
from itertools import combinations, permutations
from math import comb

R = 4


def colex_rank(e):
    return sum(comb(v, i + 1) for i, v in enumerate(sorted(e)))


def colex_decode(mask, n):
    """producer colex mask on [n] -> set of frozensets"""
    mask = int(mask)
    return {frozenset(e) for e in combinations(range(n), R) if (mask >> colex_rank(e)) & 1}


def colex_encode(es):
    return sum(1 << colex_rank(e) for e in es)


E7 = list(combinations(range(7), R))
I7 = {frozenset(e): i for i, e in enumerate(E7)}


def own_decode(mask, n=7):
    mask = int(mask)
    return {frozenset(e) for i, e in enumerate(combinations(range(n), R)) if (mask >> i) & 1}


def own_encode(es, n=7):
    idx = {frozenset(e): i for i, e in enumerate(combinations(range(n), R))}
    return sum(1 << idx[frozenset(e)] for e in es)


def flag_edge_order(s, m):
    allE = list(combinations(range(m), R))
    root = [e for e in allE if set(e) <= set(range(s))]
    rest = [e for e in allE if not set(e) <= set(range(s))]
    return root + rest


def enc(es, order):
    idx = {frozenset(e): i for i, e in enumerate(order)}
    return sum(1 << idx[frozenset(e)] for e in es)


def dec(x, order):
    return {frozenset(order[i]) for i in range(len(order)) if (x >> i) & 1}


def canon_flag(es, s, m, order):
    """canonical form of a sigma-flag: min own flag-mask over all permutations of the NON-root vertices"""
    best = None
    for pu in permutations(range(s, m)):
        g = list(range(s)) + list(pu)
        x = enc({frozenset(g[v] for v in e) for e in es}, order)
        if best is None or x < best:
            best = x
    return best


def admissible(es, n, p=5):
    """every p-subset of range(n) contains at least one edge"""
    return all(any(frozenset(f) in es for f in combinations(P, R)) for P in combinations(range(n), p))


def objective_counts(es, n):
    """(e, P): e = #edges, P = #{(T, x, y): T 3-set, x != y outside T, T+x and T+y edges} (explicit loop)"""
    P = 0
    for T in combinations(range(n), 3):
        out = [v for v in range(n) if v not in T]
        for x in out:
            for y in out:
                if x != y and frozenset(T + (x,)) in es and frozenset(T + (y,)) in es:
                    P += 1
    return len(es), P


def f_value(es, n):
    """f(H) = 2 d(H) - P[T+x, T+y in H], (T, x, y) uniform, x != y outside T  (exact Fraction)"""
    e, P = objective_counts(es, n)
    return 2 * Fraction(e, comb(n, R)) - Fraction(P, comb(n, 3) * (n - 3) * (n - 4))
