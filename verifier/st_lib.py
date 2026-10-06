"""Independent checker library (independent of the producer's glib.py / local_global.py / c10_check.py).

Encoding: a 4-graph on range(n) is a Python int; bit i <-> the i-th 4-subset in COLEX order, i.e. the 4-set
x0<x1<x2<x3 has index C(x0,1)+C(x1,2)+C(x2,3)+C(x3,4) (combinatorial number system).  This is the encoding of the
certificate's class masks (checked in r_cf.py on the complete graph and on the co-Fano class).
3-graphs (links) use the same colex rule for 3-subsets.

Giraud system Gamma(A,B,M) (paper, Definition def:giraud): edges = 4-subsets of A, 4-subsets of B and cross sets
{a,a',b,b'} with EVEN minor M_ab+M_ab'+M_a'b+M_a'b'.  B may be empty.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb

import numpy as np


def cidx(t) -> int:
    """colex index of a sorted tuple."""
    return sum(comb(x, i + 1) for i, x in enumerate(t))


@lru_cache(maxsize=None)
def subsets(n: int, r: int) -> tuple:
    """r-subsets of range(n), listed in colex order (position i has colex index i)."""
    out = [None] * comb(n, r)
    for t in combinations(range(n), r):
        out[cidx(t)] = t
    return tuple(out)


@lru_cache(maxsize=None)
def q4index(n: int) -> dict:
    return {t: cidx(t) for t in combinations(range(n), 4)}


def bit(mask: int, t) -> int:
    return (mask >> cidx(tuple(sorted(t)))) & 1


# ---------------------------------------------------------------------------------------------------------------
# Giraud systems
def giraud(n: int, A, M) -> int:
    """Gamma(A, V-A, M) on range(n); M: dict (a,b) -> 0/1 (missing entries = 0)."""
    A = frozenset(A)
    m = 0
    for i, q in enumerate(subsets(n, 4)):
        ins = [v for v in q if v in A]
        k = len(ins)
        if k == 0 or k == 4:
            m |= 1 << i
        elif k == 2:
            outs = [v for v in q if v not in A]
            p = 0
            for a in ins:
                for b in outs:
                    p ^= M.get((a, b), 0)
            if p == 0:
                m |= 1 << i
    return m


def unordered_partitions(n: int):
    """All unordered partitions {A, B} of range(n) with B possibly empty; A is the part containing 0."""
    for am in range(1 << n):
        if am & 1:
            A = tuple(v for v in range(n) if (am >> v) & 1)
            B = tuple(v for v in range(n) if not (am >> v) & 1)
            yield A, B


def _pack_rows(boolrows: np.ndarray) -> list:
    """bool array (R, N) -> list of Python ints with bit i = column i."""
    pk = np.packbits(boolrows.astype(np.uint8), axis=1, bitorder="little")
    return [int.from_bytes(r.tobytes(), "little") for r in pk]


def enum_giraud(n: int, normalised: bool = True) -> dict:
    """dict mask -> set of unordered partitions (A as frozenset, A containing 0) that represent it.
    normalised=True : M ranges over matrices with M[a][b] = 0 when a = max A or b = max B (one per switching class);
    normalised=False: M ranges over ALL 0/1 matrices on A x B (definitional; only feasible for small n)."""
    quads = subsets(n, 4)
    N = len(quads)
    out = {}
    for A, B in unordered_partitions(n):
        As = set(A)
        if normalised and B:
            free = [(a, b) for a in A if a != max(A) for b in B if b != max(B)]
        else:
            free = [(a, b) for a in A for b in B]
        fidx = {p: i for i, p in enumerate(free)}
        k = len(free)
        X = np.arange(1 << k, dtype=np.int64)
        rows = np.zeros((1 << k, N), dtype=bool)
        for i, q in enumerate(quads):
            ins = [v for v in q if v in As]
            if len(ins) in (0, 4):
                rows[:, i] = True
            elif len(ins) == 2:
                outs = [v for v in q if v not in As]
                par = np.zeros(1 << k, dtype=np.int64)
                for a in ins:
                    for b in outs:
                        j = fidx.get((a, b))
                        if j is not None:
                            par ^= (X >> j) & 1
                rows[:, i] = par == 0
        fa = frozenset(A)
        for mk in _pack_rows(rows):
            out.setdefault(mk, set()).add(fa)
    return out


def is_giraud_forced(mask: int, n: int) -> list:
    """Definitional test: list of unordered partitions {A,B} (A containing 0) such that mask = Gamma(A,B,M) for the
    matrix M forced by mask (normalised at a0 = max A, b0 = max B: M_ab = 1 - H({a0,a,b0,b}))."""
    reps = []
    for A, B in unordered_partitions(n):
        M = {}
        if B:
            a0, b0 = max(A), max(B)
            for a in A:
                for b in B:
                    if a != a0 and b != b0:
                        M[(a, b)] = 1 - bit(mask, (a0, a, b0, b))
        if giraud(n, A, M) == mask:
            reps.append(frozenset(A))
    return reps


@lru_cache(maxsize=None)
def five_sets(n: int) -> tuple:
    """for every 5-set (lex order): tuple of the 5 colex indices of its 4-subsets."""
    return tuple(tuple(cidx(q) for q in combinations(S, 4)) for S in combinations(range(n), 5))


def five_counts(mask: int, n: int) -> list:
    return [sum((mask >> i) & 1 for i in idx) for idx in five_sets(n)]


def is_odd(mask: int, n: int) -> bool:
    return all(c % 2 == 1 for c in five_counts(mask, n))


def admissible(mask: int, n: int) -> bool:
    return all(c >= 1 for c in five_counts(mask, n))


def is_giraud_struct(mask: int, n: int) -> list:
    """Structural test (different from the definitional one): {A,B} represents mask iff
      (i) every 4-subset of A and of B is an edge, (ii) no 4-set with 3 points in one part and 1 in the other is an
      edge, (iii) every 5-set with 3 points in one part and 2 in the other spans an odd number of edges.
    (iii) says that c = 1 - H on the cross 4-sets is additive in rows and in columns, hence is the minor function of
    the matrix M_ab = c(a0 a, b0 b) (proof in REPORT.txt, item 2)."""
    reps = []
    fs = list(combinations(range(n), 5))
    fsi = five_sets(n)
    for A, B in unordered_partitions(n):
        As = set(A)
        ok = True
        for i, q in enumerate(subsets(n, 4)):
            k = sum(1 for v in q if v in As)
            e = (mask >> i) & 1
            if (k in (0, 4) and not e) or (k in (1, 3) and e):
                ok = False
                break
        if not ok:
            continue
        for S, idx in zip(fs, fsi):
            k = sum(1 for v in S if v in As)
            if k in (2, 3):
                if sum((mask >> i) & 1 for i in idx) % 2 == 0:
                    ok = False
                    break
        if ok:
            reps.append(frozenset(A))
    return reps


# ---------------------------------------------------------------------------------------------------------------
# objective f = E_S g(e_S), g(e) = e(9-e)/20 and its definitional form 2 d - P[T+x, T+y in H]
G5 = [Fraction(e * (9 - e), 20) for e in range(6)]


def Lambda(mask: int, n: int) -> Fraction:
    return sum((G5[c] for c in five_counts(mask, n)), Fraction(0))


def f5(mask: int, n: int) -> Fraction:
    return Lambda(mask, n) / comb(n, 5)


def f_def(mask: int, n: int) -> Fraction:
    """2 d(H) - P[T+x in H and T+y in H], T uniform 3-set, (x, y) uniform ordered pair of distinct points off T."""
    e = bin(mask).count("1")
    d = Fraction(e, comb(n, 4))
    tot = 0
    for T in combinations(range(n), 3):
        c = sum(bit(mask, T + (x,)) for x in range(n) if x not in T)
        tot += c * (c - 1)
    P = Fraction(tot, comb(n, 3) * (n - 3) * (n - 4))
    return 2 * d - P


def induced(mask: int, n: int, S) -> int:
    S = sorted(S)
    out = 0
    for i, q in enumerate(subsets(len(S), 4)):
        if (mask >> cidx(tuple(S[j] for j in q))) & 1:
            out |= 1 << i
    return out


def relabel(mask: int, n: int, perm) -> int:
    out = 0
    for i, q in enumerate(subsets(n, 4)):
        if (mask >> i) & 1:
            out |= 1 << cidx(tuple(sorted(perm[v] for v in q)))
    return out


# ---------------------------------------------------------------------------------------------------------------
# isomorphism classes by orbit peeling (no canonical form needed)
def perm_quad_table(n: int, chunk: int = 40320) -> np.ndarray:
    """(n!, C(n,4)) uint8 array: row p, column i = colex index of the image of quad i under the p-th permutation of
    itertools.permutations(range(n)).  Built in chunks (n <= 9)."""
    from itertools import islice, permutations
    quads = subsets(n, 4)
    N = len(quads)
    assert N <= 256
    Q = np.array(quads, dtype=np.int64)                      # (N, 4)
    binom = np.array([[comb(x, r + 1) for r in range(4)] for x in range(n)], dtype=np.int64)
    it = permutations(range(n))
    blocks = []
    while True:
        P = np.array(list(islice(it, chunk)), dtype=np.int8)
        if P.size == 0:
            break
        img = P[:, Q]                                        # (c, N, 4) int8
        img.sort(axis=2)
        img = img.astype(np.int64)
        blocks.append(sum(binom[img[:, :, r], r] for r in range(4)).astype(np.uint8))
    return np.concatenate(blocks)


def orbits(family: set, n: int, table: np.ndarray = None, chunk: int = 45360):
    """Partition a relabelling-closed set of labelled 4-graphs into orbits under S_n (orbit peeling).  Returns a
    list of (rep, orbit size); asserts that every image lies in the family (closure under relabelling)."""
    if table is None:
        table = perm_quad_table(n)
    N = table.shape[1]
    remaining = set(family)
    res = []
    while remaining:
        g = min(remaining)
        row = np.array([(g >> i) & 1 for i in range(N)], dtype=bool)
        orb = set()
        for c0 in range(0, table.shape[0], chunk):
            T = table[c0:c0 + chunk].astype(np.intp)
            imgs = np.zeros(T.shape, dtype=bool)
            np.put_along_axis(imgs, T, np.broadcast_to(row, T.shape), axis=1)
            pk = np.unique(np.packbits(imgs.astype(np.uint8), axis=1, bitorder="little"), axis=0)
            orb.update(int.from_bytes(r.tobytes(), "little") for r in pk)
        assert orb <= family, "family not closed under relabelling"
        remaining -= orb
        res.append((g, len(orb)))
    return res
