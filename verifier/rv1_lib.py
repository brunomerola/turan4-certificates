"""Independent checker -- own library, written from the definitions in the paper.

No producer module is imported.  Conventions of THIS file (deliberately different from the producer's):
  * a 4-graph on a vertex set contained in {0..6} is an int bitmask, bit i = i-th 4-subset of {0..6} in
    LEXICOGRAPHIC order (L7 below); producer / certificate masks use COLEX order and are converted with an own
    colex rank (colex_to_own / own_to_colex);
  * the labelled code of a sequence lab = (lab[0], ..., lab[m-1]) of distinct vertices of H is the bitmask over the
    4-subsets of [m] in lexicographic order: bit j set iff {lab[i] : i in SUB[m][j]} is an edge of H;
  * a labelled type tau on [s] is the code of theta = (theta(0), ..., theta(s-1)); H[theta] = tau means equality of
    codes (exact labelled test, no canonical root labelling, no automorphism tables);
  * a sigma-flag (roots 0..s-1, non-roots s..m-1) is identified up to isomorphism fixing the roots pointwise by the
    canonical code = min over all orderings of the non-root vertices of the labelled code.
Moment counts: for each ORDERED injective theta: [s] -> V(H), each ordered pair (U1, U2) of DISJOINT (m-s)-subsets of
V(H) - theta, the pair (canon(theta, U1), canon(theta, U2)) is counted under the labelled type code(theta).  The
number of configurations is conf(s, m) = 7!/(7-s)! * C(7-s, m-s) * C(7-m, m-s) (asserted).
"""
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial

N = 7
SUB = {m: list(combinations(range(m), 4)) for m in range(N + 1)}
L7 = SUB[N]                                   # own lex order of the 35 4-subsets of {0..6}
IX = [-1] * (N ** 4)                          # ordered quadruple of distinct vertices -> own lex index
for _i, _e in enumerate(L7):
    for _p in permutations(_e):
        IX[((_p[0] * N + _p[1]) * N + _p[2]) * N + _p[3]] = _i
FIVE = [[i for i, e in enumerate(L7) if set(e) <= set(P)] for P in combinations(range(N), 5)]


def colex_rank(e):
    return sum(comb(v, i + 1) for i, v in enumerate(sorted(e)))


def colex_to_own(mask, n=N):
    """Colex-order mask of a 4-graph on [n] (n <= 7) -> own lex mask on [7] (vertex set kept)."""
    h = 0
    for e in combinations(range(n), 4):
        if (mask >> colex_rank(e)) & 1:
            h |= 1 << IX[((e[0] * N + e[1]) * N + e[2]) * N + e[3]]
    return h


def own_to_colex(h):
    x = 0
    for i, e in enumerate(L7):
        if (h >> i) & 1:
            x |= 1 << colex_rank(e)
    return x


def edges_of(h):
    return [L7[i] for i in range(35) if (h >> i) & 1]


def popcount(x):
    return bin(x).count("1")


def five_profile(h):
    """List of edge counts of the 21 five-subsets of {0..6}."""
    return [sum((h >> i) & 1 for i in F) for F in FIVE]


def admissible5(h):
    return all(c >= 1 for c in five_profile(h))


def is_odd(h):
    return all(c % 2 == 1 for c in five_profile(h))


def code(h, lab, m):
    x = 0
    for j, (a, b, c, d) in enumerate(SUB[m]):
        if (h >> IX[((lab[a] * N + lab[b]) * N + lab[c]) * N + lab[d]]) & 1:
            x |= 1 << j
    return x


def canon(h, th, U, m):
    return min(code(h, th + p, m) for p in permutations(U))


def conf(s, m):
    return factorial(N) // factorial(N - s) * comb(N - s, m - s) * comb(N - m, m - s)


def moment_counts(h, s, m):
    """{labelled type code on [s]: {(canon1, canon2): count}} over all configurations (theta, U1, U2) of H."""
    k = m - s
    out = {}
    tot = 0
    for th in permutations(range(N), s):
        t = code(h, th, s)
        rest = tuple(v for v in range(N) if v not in th)
        Us = list(combinations(rest, k))
        fl = [canon(h, th, U, m) for U in Us]
        sets = [set(U) for U in Us]
        d = out.setdefault(t, {})
        for i in range(len(Us)):
            for j in range(len(Us)):
                if sets[i].isdisjoint(sets[j]) and (k > 0 or i == j):
                    key = (fl[i], fl[j])
                    d[key] = d.get(key, 0) + 1
                    tot += 1
    assert tot == conf(s, m), (s, m, tot)
    return out


def flag_universe(s, m, tau, p=5):
    """Sorted canonical codes of all p-admissible tau-flags on m vertices (tau: own lex code on [s])."""
    root = [SUB[m].index(SUB[s][j]) for j in range(len(SUB[s])) if (tau >> j) & 1]
    nonroot = [j for j, e in enumerate(SUB[m]) if not set(e) <= set(range(s))]
    fives = [[j for j, e in enumerate(SUB[m]) if set(e) <= set(P)] for P in combinations(range(m), p)]
    ident = tuple(range(s))
    U = tuple(range(s, m))
    out = set()
    for bits in range(1 << len(nonroot)):
        es = set(SUB[m][j] for j in root) | {SUB[m][nonroot[i]] for i in range(len(nonroot)) if (bits >> i) & 1}
        if not all(any(SUB[m][j] in es for j in F) for F in fives):
            continue
        # embed the flag graph on [m] into [7] (same vertex names) as an own mask
        hf = 0
        for e in es:
            hf |= 1 << IX[((e[0] * N + e[1]) * N + e[2]) * N + e[3]]
        assert code(hf, ident, s) == tau
        out.add(canon(hf, ident, U, m))
    return sorted(out)


def all_labelled_types(s):
    return list(range(1 << len(SUB[s])))


def psd_exact(rows):
    """Exact PSD test of a symmetric matrix (list of lists of ints/Fractions): symmetric Gaussian elimination.
    A is PSD iff a_11 >= 0, a_11 = 0 forces row 1 = 0 (then recurse), and for a_11 > 0 the Schur complement is PSD.
    Returns (is_psd, rank, reason, min_pivot_over_diag (Fraction or None))."""
    n = len(rows)
    for i in range(n):
        for j in range(i):
            if rows[i][j] != rows[j][i]:
                return False, None, "not symmetric", None
    A = [[Fraction(x) for x in r] for r in rows]
    diag0 = [A[i][i] for i in range(n)]
    rank = 0
    worst = None
    for k in range(n):
        p = A[k][k]
        if p < 0:
            return False, rank, f"negative pivot at {k}", None
        if p == 0:
            if any(A[k][j] != 0 for j in range(k + 1, n)):
                return False, rank, f"zero pivot with nonzero row at {k}", None
            continue
        rank += 1
        r = p / diag0[k]
        worst = r if worst is None or r < worst else worst
        rowk = A[k]
        for i in range(k + 1, n):
            a = A[i][k]
            if a == 0:
                continue
            f = a / p
            Ai = A[i]
            for j in range(i, n):
                if rowk[j]:
                    Ai[j] -= f * rowk[j]
            for j in range(i + 1, n):     # keep symmetry explicitly
                A[j][i] = Ai[j]
    return True, rank, "ok", worst


def canonical_graph(h, perm_maps):
    """min over all vertex permutations of the relabelled own mask (perm_maps from perm_edge_maps())."""
    best = None
    bits = [i for i in range(35) if (h >> i) & 1]
    for pm in perm_maps:
        x = 0
        for i in bits:
            x |= 1 << pm[i]
        if best is None or x < best:
            best = x
    return best


def perm_edge_maps():
    out = []
    for g in permutations(range(N)):
        out.append([IX[((g[a] * N + g[b]) * N + g[c]) * N + g[d]] for (a, b, c, d) in L7])
    return out
