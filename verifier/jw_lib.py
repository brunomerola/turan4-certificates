"""Independent checker library (INDEPENDENT re-implementation; written by the reviewer from the definitions).

Imports nothing from the producer code or from the flag-algebra engine.  Python ints /
fractions.Fraction for everything that carries the claim; numpy only for (a) the min-over-S_7 canonical forms
(float64 matmul of 0/1 bits by powers of two <= 2^34: every partial sum is an integer < 2^35, hence exact) and (b)
the mod-p rank estimates, which are used only as upper bounds on null spaces (rank_p <= rank_Q).

Own conventions (deliberately different from the producer's):
  * a 3-graph on [n] is an int whose bit i is the i-th triple of [n] in LEXICOGRAPHIC order (itertools.combinations);
    the witness files use colex bits {a<b<c} -> a + C(b,2) + C(c,3), converted on input by colex_to_lex();
  * a flag of block (s, m) at (theta, U) is coded by the bits of the triples of [m] that meet {s..m-1}, lex order,
    minimised over the orderings of U (theta_i -> i, U -> s..m-1); the type is the lex mask of H[theta] on [s];
  * the 7-vertex law of an iterated blow-up is computed by summing over SET PARTITIONS of the sample (first level at
    which the sample is not inside one part) times injective part labellings, not over functions.
"""
from __future__ import annotations

import itertools
import math
from collections import Counter, defaultdict
from fractions import Fraction

import numpy as np

N = 7
BLOCKS = [(s, m) for s in range(N) for m in range(3, N + 1) if m > s and 2 * m - s <= N]
assert BLOCKS == [(0, 3), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4), (3, 5), (4, 5), (5, 6)], BLOCKS

TRI = {n: list(itertools.combinations(range(n), 3)) for n in range(0, N + 1)}
TIDX = {n: {t: i for i, t in enumerate(TRI[n])} for n in range(0, N + 1)}


def mask_of(tris, n=N) -> int:
    m = 0
    for t in tris:
        m |= 1 << TIDX[n][tuple(sorted(t))]
    return m


def edges_of(mask, n=N):
    return [t for i, t in enumerate(TRI[n]) if (mask >> i) & 1]


# ------------------------------------------------------------------ witness-file encoding (colex) -> own (lex)
_COLEX = {t: t[0] + math.comb(t[1], 2) + math.comb(t[2], 3) for t in TRI[N]}
assert sorted(_COLEX.values()) == list(range(35))


def colex_to_lex(cmask: int) -> int:
    return mask_of([t for t in TRI[N] if (cmask >> _COLEX[t]) & 1])


def lex_to_colex(mask: int) -> int:
    return sum(1 << _COLEX[t] for t in edges_of(mask))


# ------------------------------------------------------------------ canonical forms (min over S_7 of the lex mask)
_PERMS = list(itertools.permutations(range(N)))
_PT = np.array([[TIDX[N][tuple(sorted((p[a], p[b], p[c])))] for (a, b, c) in TRI[N]] for p in _PERMS],
               dtype=np.int64)          # (5040, 35): image bit of each bit


def canon_many(masks):
    masks = np.asarray(list(masks), dtype=np.int64)
    if len(masks) == 0:
        return []
    B = ((masks[:, None] >> np.arange(35)) & 1).astype(np.float64)
    best = np.full(len(masks), np.inf)
    for c0 in range(0, len(_PERMS), 252):
        W = np.exp2(_PT[c0:c0 + 252].T.astype(np.float64))     # (35, chunk)
        best = np.minimum(best, (B @ W).min(axis=1))
    out = [int(x) for x in best]
    for m0, c in zip(masks.tolist()[:3], out[:3]):              # spot check exactness
        assert c == min(sum(1 << int(_PT[k, i]) for i in range(35) if (m0 >> i) & 1) for k in range(len(_PERMS)))
    return out


# ------------------------------------------------------------------ constructions (own definitions)
def fano_lines():
    """PG(2,2): points = nonzero vectors of F_2^3 (1..7 -> parts 0..6), lines = {a, b, a xor b}."""
    L = {frozenset((a - 1, b - 1, (a ^ b) - 1)) for a in range(1, 8) for b in range(a + 1, 8)}
    assert len(L) == 7 and all(len(l) == 3 for l in L)
    return L


def fano_complement():
    L = fano_lines()
    return [t for t in TRI[7] if frozenset(t) not in L]


def s6_design():
    """The 10-triple 3-graph on [6] in which every pair lies in exactly 2 triples and every 4-set spans 0 or 2
    triples (Frankl-Furedi S_6); found by search, first hit (unique up to isomorphism)."""
    T6 = TRI[6]
    for E in itertools.combinations(T6, 10):
        cnt = Counter(p for t in E for p in itertools.combinations(t, 2))
        if len(cnt) != 15 or any(v != 2 for v in cnt.values()):
            continue
        S = set(E)
        if all(sum(t in S for t in itertools.combinations(Q, 3)) in (0, 2) for Q in itertools.combinations(range(6), 4)):
            return list(E)
    raise RuntimeError


def set_partitions(k):
    def rec(i, blocks):
        if i == k:
            yield [list(b) for b in blocks]
            return
        for b in blocks:
            b.append(i)
            yield from rec(i + 1, blocks)
            b.pop()
        blocks.append([i])
        yield from rec(i + 1, blocks)
        blocks.pop()
    yield from rec(0, [])


def law_iterated(pattern, p, kmax=N):
    """Exact law of the k-vertex sample (k <= kmax) of the iterated blow-up of the p-vertex 3-graph `pattern`
    (parts of weight 1/p, transversal triples over pattern edges, recursion inside every part).  For the first level
    at which the k points are not all in one part: P(that level's map is g) = 1/(p^k - p) for each non-constant
    g: [k] -> [p] (self-similarity), the triples inside a fibre are an independent k'-sample of the same law."""
    pat = [[[False] * p for _ in range(p)] for _ in range(p)]
    for t in pattern:
        for a, b, c in itertools.permutations(t):
            pat[a][b][c] = True
    D = {0: {0: Fraction(1)}, 1: {0: Fraction(1)}, 2: {0: Fraction(1)}}
    for k in range(3, kmax + 1):
        acc = defaultdict(Fraction)
        ntot = 0
        for blocks in set_partitions(k):
            j = len(blocks)
            if j < 2:
                continue
            W = {0: Fraction(1)}
            for B in blocks:                     # blocks are increasing lists
                if len(B) >= 3:
                    sub = {}
                    for msk, pr in D[len(B)].items():
                        mm = 0
                        for (a, b, c) in edges_of(msk, len(B)):
                            mm |= 1 << TIDX[k][(B[a], B[b], B[c])]
                        sub[mm] = pr
                    W = {w1 | w2: p1 * p2 for w1, p1 in W.items() for w2, p2 in sub.items()}
            bt = []
            for (x, y, z) in itertools.combinations(range(j), 3):
                mm = 0
                for a in blocks[x]:
                    for b in blocks[y]:
                        for c in blocks[z]:
                            mm |= 1 << TIDX[k][tuple(sorted((a, b, c)))]
                bt.append((x, y, z, mm))
            C = Counter()
            for L in itertools.permutations(range(p), j):
                mm = 0
                for (x, y, z, tm) in bt:
                    if pat[L[x]][L[y]][L[z]]:
                        mm |= tm
                C[mm] += 1
                ntot += 1
            for cm, cnt in C.items():
                for wm, pr in W.items():
                    acc[cm | wm] += cnt * pr
        assert ntot == p ** k - p
        den = p ** k - p
        D[k] = {m: v / den for m, v in acc.items()}
        assert sum(D[k].values()) == 1
    return D


def law_flat(rule, p, k=N):
    """Exact law of the k-vertex sample of the (non-iterated) balanced blow-up: vertices uniform in [p]
    independently, triple {a,b,c} an edge iff rule(part(a), part(b), part(c))."""
    acc = Counter()
    for g in itertools.product(range(p), repeat=k):
        mm = 0
        for i, (a, b, c) in enumerate(TRI[k]):
            if rule(g[a], g[b], g[c]):
                mm |= 1 << i
        acc[mm] += 1
    return {m: Fraction(c, p ** k) for m, c in acc.items()}


def rule_Sn(a, b, c):          # complete balanced 3-partite
    return len({a, b, c}) == 3


def rule_Cn(a, b, c):          # Turan's construction: V1V2V3, ViViV(i+1)
    if len({a, b, c}) == 3:
        return True
    cn = Counter((a, b, c))
    if len(cn) == 2:
        (i, ni), (j, nj) = cn.most_common()
        return ni == 2 and j == (i + 1) % 3
    return False


def rule_Bn(a, b, c):          # two parts, triples meeting both parts
    return len({a, b, c}) == 2


def construction_law7(name):
    """Returns (dict labelled lex mask -> Fraction) for the 7-vertex sample, the problem name and description."""
    if name == "Fanoc":
        return law_iterated(fano_complement(), 7)[7], "J4", "iterated blow-up of the Fano-plane complement"
    if name == "S6it":
        return law_iterated(s6_design(), 6)[7], "K4m", "iterated blow-up of S_6 (2-(6,3,2) design)"
    if name == "edgeit":
        return law_iterated([(0, 1, 2)], 3)[7], "C5m", "iterated blow-up of an edge"
    if name == "Sn":
        return law_flat(rule_Sn, 3), "F5", "balanced complete 3-partite S_n"
    if name == "Cn":
        return law_flat(rule_Cn, 3), "K4", "Turan's construction C_n"
    if name == "Bn":
        return law_flat(rule_Bn, 2), "K5", "B_n: two halves, triples meeting both"
    raise KeyError(name)


def classes_of(law):
    """Group a labelled law into isomorphism classes: dict canonical mask -> Fraction."""
    ms = list(law)
    cs = canon_many(ms)
    out = defaultdict(Fraction)
    for m, c in zip(ms, cs):
        out[c] += law[m]
    return dict(out)


# ------------------------------------------------------------------ objective (own, from the definition)
def adj3(mask, n=N):
    E = [[[False] * n for _ in range(n)] for _ in range(n)]
    for t in edges_of(mask, n):
        for a, b, c in itertools.permutations(t):
            E[a][b][c] = True
    return E


def sigma7(mask) -> Fraction:
    """x != y codegree-squared density: sum over pairs T of c_T (c_T - 1) / (C(7,2) * 5 * 4)."""
    c = Counter(p for t in edges_of(mask) for p in itertools.combinations(t, 2))
    return Fraction(sum(v * (v - 1) for v in c.values()), math.comb(7, 2) * 5 * 4)


def sigma7_4sets(mask) -> Fraction:
    E = set(edges_of(mask))
    s = sum(math.comb(sum(t in E for t in itertools.combinations(Q, 3)), 2) for Q in itertools.combinations(range(7), 4))
    return Fraction(s, 6 * 35)


def fval(mask) -> Fraction:
    v = 1 - sigma7(mask)
    assert v == 1 - sigma7_4sets(mask)
    return v


# ------------------------------------------------------------------ forbidden families (own tests)
def contains_J4(mask):
    E = adj3(mask)
    for v in range(7):
        others = [u for u in range(7) if u != v]
        for Q in itertools.combinations(others, 4):
            if all(E[v][a][b] for a, b in itertools.combinations(Q, 2)):
                return True
    return False


def _contains_pattern(mask, pat, k):
    E = adj3(mask)
    for w in itertools.permutations(range(7), k):
        if all(E[w[a]][w[b]][w[c]] for (a, b, c) in pat):
            return True
    return False


def contains_clique(mask, k):
    S = set(edges_of(mask))
    return any(all(t in S for t in itertools.combinations(Q, 3)) for Q in itertools.combinations(range(7), k))


def admissible(mask, prob):
    if prob == "J4":
        return not contains_J4(mask)
    if prob == "F5":
        return not _contains_pattern(mask, [(0, 1, 2), (0, 1, 3), (2, 3, 4)], 5)
    if prob == "K4":
        return not contains_clique(mask, 4)
    if prob == "K5":
        return not contains_clique(mask, 5)
    if prob == "K4m":
        S = set(edges_of(mask))
        return all(sum(t in S for t in itertools.combinations(Q, 3)) < 3 for Q in itertools.combinations(range(7), 4))
    if prob == "C5m":
        return not _contains_pattern(mask, [(0, 1, 2), (1, 2, 3), (2, 3, 4), (0, 3, 4)], 5)
    raise KeyError(prob)


# ------------------------------------------------------------------ moment counts from the definition
_FTRI = {(s, m): [t for t in TRI[m] if t[2] >= s] for (s, m) in BLOCKS}


def moment_counts(mask, blocks=BLOCKS, mutate=None):
    """dict (s, m, type) -> Counter{(F1, F2): #configurations}; configuration = (theta ordered s-tuple of distinct
    vertices, ordered pair of disjoint (m-s)-subsets U1, U2 of the rest); type = lex mask of H[theta] on [s]."""
    E = adj3(mask)
    out = {}
    for (s, m) in blocks:
        r = m - s
        ft = _FTRI[(s, m)]
        st = TRI[s]
        for th in itertools.permutations(range(N), s):
            tm = 0
            for i, (a, b, c) in enumerate(st):
                if E[th[a]][th[b]][th[c]]:
                    tm |= 1 << i
            rest = [v for v in range(N) if v not in th]
            codes = {}
            for U in itertools.combinations(rest, r):
                best = None
                for oU in itertools.permutations(U):
                    w = th + oU
                    code = 0
                    for i, (a, b, c) in enumerate(ft):
                        if E[w[a]][w[b]][w[c]]:
                            code |= 1 << i
                    if best is None or code < best:
                        best = code
                codes[U] = best
            if mutate == "typeblind":          # MUTATION (sensitivity test only): forget the type
                tm = -1
            cnt = out.setdefault((s, m, tm), Counter())
            Us = list(codes)
            for U1 in Us:
                s1 = set(U1)
                for U2 in Us:
                    if s1.isdisjoint(U2) or mutate == "overlap":   # MUTATION: allow U1, U2 to overlap
                        cnt[(codes[U1], codes[U2])] += 1
    return out


def n_configs(s, m):
    return math.perm(N, s) * math.comb(N - s, m - s) * math.comb(N - m, m - s)


# ------------------------------------------------------------------ exact linear algebra
def psd_bareiss(A):
    """Exact PSD test of a symmetric integer matrix by fraction-free symmetric elimination in the given order.
    Entry (a, b) after eliminating the pivot set P equals det(A_PP) * (Schur complement)_ab, det(A_PP) > 0, so signs
    are exact.  PSD iff every current diagonal entry is >= 0 and a zero diagonal entry has a zero row.  Returns
    (psd, rank, pivots, zero_rows, first_failure)."""
    n = len(A)
    M = [list(r) for r in A]
    for i in range(n):
        for j in range(n):
            assert M[i][j] == M[j][i]
    rem = list(range(n))
    prev = 1
    piv, zero = [], []
    while rem:
        i = rem.pop(0)
        d = M[i][i]
        if d < 0:
            return False, None, piv, zero, ("negative pivot", i)
        if d == 0:
            if any(M[i][j] != 0 for j in rem):
                return False, None, piv, zero, ("zero pivot, nonzero row", i)
            zero.append(i)
            continue
        Mi = M[i]
        for a in rem:
            Ma = M[a]
            mai = Ma[i]
            for b in rem:
                if b < a:
                    continue
                q, rr = divmod(Ma[b] * d - mai * Mi[b], prev)
                assert rr == 0
                Ma[b] = q
                M[b][a] = q
        prev = d
        piv.append(i)
    return True, len(piv), piv, zero, None


def rref_nullspace(A):
    """Exact null space basis (list of integer column vectors, primitive) of a rational matrix A (list of rows)."""
    rows = [[Fraction(x) for x in r] for r in A]
    nr = len(rows)
    nc = len(rows[0]) if nr else 0
    piv = []
    r = 0
    for c in range(nc):
        p = next((i for i in range(r, nr) if rows[i][c] != 0), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        pv = rows[r][c]
        rows[r] = [x / pv for x in rows[r]]
        Rr = rows[r]
        for i in range(nr):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [x - f * y for x, y in zip(rows[i], Rr)]
        piv.append(c)
        r += 1
        if r == nr:
            break
    free = [c for c in range(nc) if c not in set(piv)]
    basis = []
    for fc in free:
        v = [Fraction(0)] * nc
        v[fc] = Fraction(1)
        for i, pc in enumerate(piv):
            v[pc] = -rows[i][fc]
        L = 1
        for x in v:
            L = L * x.denominator // math.gcd(L, x.denominator)
        iv = [int(x * L) for x in v]
        g = 0
        for x in iv:
            g = math.gcd(g, x)
        basis.append([x // g for x in iv])
    return basis, len(piv)


def lcm_den(fracs):
    L = 1
    for x in fracs:
        L = L * x.denominator // math.gcd(L, x.denominator)
    return L


def rank_mod_p(rows, p):
    """Rank mod p of an integer matrix (numpy int64, entries already reduced mod p), Gaussian elimination."""
    A = np.array(rows, dtype=np.int64) % p
    nr, nc = A.shape
    r = 0
    pivrows = []
    idx = np.arange(nr)
    for c in range(nc):
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0:
            continue
        k = r + nz[0]
        A[[r, k]] = A[[k, r]]
        idx[[r, k]] = idx[[k, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        col = A[:, c].copy()
        col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr):
            A[nzr] = (A[nzr] - (col[nzr, None] * A[r][None, :]) % p) % p
        pivrows.append(int(idx[r]))
        r += 1
        if r == nr:
            break
    return r, pivrows
