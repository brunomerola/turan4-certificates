"""Independent checker: exact value of the sharp six-vertex certificates (Remark n6sharp), from the definitions.

Own code, written from Definition cert of the paper; nothing is imported from the producer code (search/) or from
the other checkers. Reads only the certificate JSON: r = 4, N = 6, p, allowed (admissible edge counts of a p-set;
null = at least one), bound, and blocks with s, m, sigma (colex mask on [s]), flags (colex masks on [m], roots
0..s-1) and the rational matrix Q = Qnum / den.

Checks:
  * every Q is symmetric and positive semidefinite: an exact LDL^T decomposition over the rationals is computed and
    the product L D L^T is multiplied back and compared with Q entry by entry (D >= 0), so a PASS does not rely on
    the correctness of the elimination;
  * 2m - s <= 6 for every block; every listed flag is a sigma-flag (its restriction to the roots is sigma) and the
    listed flags of a block are pairwise non-isomorphic (isomorphisms fixing the roots, exhaustive search);
  * all 2^15 labelled 4-graphs on {0..5} are enumerated (own lexicographic edge order); admissible = every p-subset
    spans an allowed number of edges; isomorphism classes by orbit marking over all 720 vertex permutations (the
    orbit sizes must add up to the number of labelled admissible graphs);
  * for every class representative H, c_Q(H) = (1/conf) * sum over ordered injective theta: [s] -> V(H) with
    H[theta] = sigma exactly, over ordered pairs (U1, U2) of disjoint (m-s)-subsets of V(H) - theta([s]), of
    Q(H^U1_theta, H^U2_theta), where a flag not isomorphic to a listed one counts 0, and
    conf = 6!/(6-s)! * C(6-s, m-s) * C(6-m, m-s) (asserted equal to the number of configurations);
    value(H) = e(H)/15 - sum over blocks of c_Q(H), in exact rational arithmetic; b = min over the classes;
  * self-consistency: every representative is re-evaluated after a random relabelling of its vertices.
Usage: python rv_sharp6.py CERT_JSON
"""
import json
import random
import sys
import time
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial

R, N = 4, 6
t0 = time.time()
cert = json.load(open(sys.argv[1]))
assert cert["r"] == R and cert["N"] == N, (cert["r"], cert["N"])
p = cert["p"]
allowed = cert["allowed"]
claimed = Fraction(cert["bound"])


def ok_count(c):
    return c >= 1 if allowed is None else c in allowed


# ---- encodings ------------------------------------------------------------------------------------------------------
def colex_sets(n):
    """The 4-subsets of {0..n-1} in colex order (compare the largest elements first)."""
    return sorted(combinations(range(n), R), key=lambda e: tuple(reversed(e)))


def from_colex(mask, n):
    """Colex bit mask of the certificate -> set of 4-tuples (sorted)."""
    return {e for i, e in enumerate(colex_sets(n)) if (mask >> i) & 1}


LEX6 = list(combinations(range(N), R))             # own edge order on {0..5}: bit i <-> LEX6[i]
POS6 = {e: i for i, e in enumerate(LEX6)}


def graph_of(x):
    return {LEX6[i] for i in range(len(LEX6)) if (x >> i) & 1}


def relabel(es, g):
    """Image of an edge set under the vertex map g (a list or dict)."""
    return {tuple(sorted(g[v] for v in e)) for e in es}


def admissible(es, n):
    return all(ok_count(sum(1 for e in combinations(P, R) if e in es)) for P in combinations(range(n), p))


def flag_key(es, s, m):
    """Canonical form of a sigma-flag on [m] with roots 0..s-1: the least sorted edge list over all orderings of
    the non-roots (roots fixed pointwise)."""
    best = None
    for perm in permutations(range(s, m)):
        g = list(range(s)) + list(perm)
        k = tuple(sorted(relabel(es, g)))
        if best is None or k < best:
            best = k
    return best


# ---- exact PSD test -------------------------------------------------------------------------------------------------
def ldl_psd(Q):
    """Q: list of lists of Fractions. Returns (symmetric, psd). psd is True only if Q = L D L^T was verified exactly
    with L unit lower triangular and D >= 0."""
    n = len(Q)
    if any(Q[i][j] != Q[j][i] for i in range(n) for j in range(n)):
        return False, False
    A = [row[:] for row in Q]
    L = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    D = [Fraction(0)] * n
    for k in range(n):
        d = A[k][k]
        if d < 0:
            return True, False
        if d == 0:
            if any(A[i][k] != 0 for i in range(k + 1, n)):
                return True, False
            continue
        D[k] = d
        for i in range(k + 1, n):
            L[i][k] = A[i][k] / d
        for i in range(k + 1, n):
            if L[i][k]:
                for j in range(k + 1, n):
                    A[i][j] -= L[i][k] * d * L[j][k]
    back = all(sum(L[i][k] * D[k] * L[j][k] for k in range(n)) == Q[i][j] for i in range(n) for j in range(n))
    return True, back and all(x >= 0 for x in D)


# ---- the certificate ------------------------------------------------------------------------------------------------
blocks = []
all_sym = all_psd = flags_ok = True
for b in cert["blocks"]:
    s, m = int(b["s"]), int(b["m"])
    assert 0 <= s <= m <= N and 2 * m - s <= N, (s, m)
    den = int(b["den"])
    assert den > 0
    Qn = [[int(v) for v in row] for row in b["Qnum"]]
    n = len(b["flags"])
    assert len(Qn) == n and all(len(row) == n for row in Qn)
    sym, psd = ldl_psd([[Fraction(v, den) for v in row] for row in Qn])
    all_sym &= sym
    all_psd &= psd
    assert 0 <= int(b["sigma"]) < (1 << comb(s, R)), "sigma mask out of range"
    assert all(0 <= int(f) < (1 << comb(m, R)) for f in b["flags"]), "flag mask out of range"
    sig = from_colex(int(b["sigma"]), s)
    keys = {}
    for i, f in enumerate(b["flags"]):
        F = from_colex(int(f), m)
        if {e for e in F if max(e) < s} != sig:
            flags_ok = False
        k = flag_key(F, s, m)
        if k in keys:
            flags_ok = False
        keys[k] = i
    conf = factorial(N) // factorial(N - s) * comb(N - s, m - s) * comb(N - m, m - s)
    nadm = sum(admissible(from_colex(int(f), m), m) for f in b["flags"])
    blocks.append((s, m, sig, keys, Qn, den, conf))
    print(f"  block (s={s}, m={m}) sigma {sorted(sig)}: {n} flags ({nadm} admissible), den {den}, symmetric {sym}, "
          f"PSD {psd}, nonzero {any(v for row in Qn for v in row)}")
print(f"Q: {len(blocks)} blocks; all symmetric: {all_sym}; all PSD (exact L D L^T, product checked): {all_psd}")
print(f"flags: all sigma-flags, pairwise non-isomorphic: {flags_ok}")


def value(es):
    """e(H)/15 - sum_blocks c_Q(H), exact."""
    tot = Fraction(len(es), comb(N, R))
    for s, m, sig, keys, Qn, den, conf in blocks:
        acc, ncf = 0, 0
        for th in permutations(range(N), s):
            # H[theta] = sigma exactly: for every 4-subset X of [s], theta(X) is an edge iff X is an edge of sigma
            if any((tuple(sorted(th[v] for v in X)) in es) != (X in sig) for X in combinations(range(s), R)):
                ncf += sum(1 for U1 in combinations([v for v in range(N) if v not in th], m - s)
                           for U2 in combinations([v for v in range(N) if v not in th], m - s) if not set(U1) & set(U2))
                continue
            rest = [v for v in range(N) if v not in th]
            idx = {}
            for U in combinations(rest, m - s):
                g = {th[i]: i for i in range(s)}
                g.update({u: s + j for j, u in enumerate(U)})
                F = {tuple(sorted(g[v] for v in e)) for e in es if set(e) <= set(g)}
                idx[U] = keys.get(flag_key(F, s, m))
            for U1, a in idx.items():
                for U2, c in idx.items():
                    if set(U1) & set(U2):
                        continue
                    ncf += 1
                    if a is not None and c is not None:
                        acc += Qn[a][c]
        assert ncf == conf, (s, m, ncf, conf)
        tot -= Fraction(acc, den * conf)
    return tot


# ---- all admissible 4-graphs on six vertices, up to isomorphism ----------------------------------------------------
PERM_MAPS = [[POS6[tuple(sorted(g[v] for v in e))] for e in LEX6] for g in permutations(range(N))]
seen = bytearray(1 << len(LEX6))
reps, nlab, orbit_total = [], 0, 0
for x in range(1 << len(LEX6)):
    es = graph_of(x)
    if not admissible(es, N):
        continue
    nlab += 1
    if seen[x]:
        continue
    bits = [i for i in range(len(LEX6)) if (x >> i) & 1]
    orbit = {sum(1 << pm[i] for i in bits) for pm in PERM_MAPS}
    for y in orbit:
        seen[y] = 1
    orbit_total += len(orbit)
    reps.append(x)
print(f"labelled admissible 4-graphs on 6 vertices: {nlab}; isomorphism classes: {len(reps)} "
      f"(orbit sizes sum to the labelled count: {orbit_total == nlab})")

rng = random.Random(20261006)
vals, relabel_ok = [], True
for x in reps:
    es = graph_of(x)
    v = value(es)
    g = list(range(N))
    rng.shuffle(g)
    relabel_ok &= value(relabel(es, g)) == v
    vals.append(v)
best = min(vals)
att = [x for x, v in zip(reps, vals) if v == best]
print(f"exact min over all {len(reps)} classes: b = {best} = {float(best):.15f}; attained by {len(att)} classes "
      f"(edge counts {sorted(len(graph_of(x)) for x in att)})")
print(f"relabelling self-check (value of a random relabelling of every class equals its value): {relabel_ok}")
print(f"claimed {claimed}; EQUAL: {best == claimed}")
ok = all_sym and all_psd and flags_ok and orbit_total == nlab and relabel_ok and best == claimed
print(f"SHARP N = 6 CERTIFICATE {'OK' if ok else 'NOT OK'} ({time.time() - t0:.0f}s)")
sys.exit(0 if ok else 1)
