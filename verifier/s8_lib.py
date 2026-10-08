"""Independent checker library (own code; no producer module is imported).

Conventions (chosen to match the certificate's DATA format, which is defined by the producer's flag ordering; the
format is re-implemented here from its description, not imported):
  * vertices 0..n-1; the r-subsets of range(n) in colex order (sort by reversed tuple); an r-graph on n vertices is
    the integer sum of 2^colex_rank(e) over its edges;
  * a type of size s = the minimum, over all s! relabellings, of the mask of an admissible s-vertex graph; types of a
    block are listed in increasing mask order (only types with at least one flag);
  * a sigma-flag on m vertices = an admissible m-vertex graph whose restriction to 0..s-1 is sigma, canonical form =
    minimum over relabellings fixing 0..s-1 pointwise; flags of a type in increasing canonical-mask order;
  * Gram matrices: G_k = A_k^T A_k (A_k = integer rows of the .cert.npz), Q_k = G_k / M^2.
Families are defined by forbidden graphs (BCL Table 2 edge lists, 1-based).  'sub' problems: G is admissible iff no
labelled copy of a forbidden graph is a subset of E(G).  'co' problems (4-graphs, complement form): H is admissible
iff its complement is F-free, i.e. iff every labelled copy of F meets E(H).
"""
from __future__ import annotations

import random
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import comb

import numpy as np


@lru_cache(maxsize=None)
def edges(n: int, r: int) -> tuple:
    return tuple(sorted(combinations(range(n), r), key=lambda e: e[::-1]))


def cidx(e) -> int:
    e = sorted(e)
    return sum(comb(v, i + 1) for i, v in enumerate(e))


def _check_colex():
    for n in range(1, 9):
        for r in range(1, min(n, 5) + 1):
            for i, e in enumerate(edges(n, r)):
                assert cidx(e) == i


_check_colex()


@lru_cache(maxsize=None)
def emap(n: int, r: int, g: tuple) -> tuple:
    """emap[i] = colex rank of g(e_i) (vertex v -> g[v])."""
    return tuple(cidx(tuple(g[v] for v in e)) for e in edges(n, r))


def apply_perm(mask: int, n: int, r: int, g) -> int:
    em = emap(n, r, tuple(g))
    out = 0
    i = 0
    m = int(mask)
    while m:
        if m & 1:
            out |= 1 << em[i]
        m >>= 1
        i += 1
    return out


def apply_perm_np(masks: np.ndarray, n: int, r: int, g) -> np.ndarray:
    em = emap(n, r, tuple(g))
    masks = np.asarray(masks, dtype=np.int64)
    out = np.zeros_like(masks)
    for i, j in enumerate(em):
        out |= ((masks >> i) & 1) << j
    return out


@lru_cache(maxsize=None)
def perms_fix(m: int, s: int) -> tuple:
    return tuple(tuple(range(s)) + p for p in permutations(range(s, m)))


def canon_all(mask: int, n: int, r: int) -> int:
    return min(apply_perm(mask, n, r, g) for g in permutations(range(n)))


def canon_fix(mask: int, m: int, r: int, s: int) -> int:
    return min(apply_perm(mask, m, r, g) for g in perms_fix(m, s))


def canon_all_np(masks, n, r):
    masks = np.asarray(masks, dtype=np.int64)
    best = masks.copy()
    for g in permutations(range(n)):
        best = np.minimum(best, apply_perm_np(masks, n, r, g))
    return best


def canon_fix_np(masks, m, r, s):
    masks = np.asarray(masks, dtype=np.int64)
    best = masks.copy()
    for g in perms_fix(m, s):
        best = np.minimum(best, apply_perm_np(masks, m, r, g))
    return best


# ------------------------------------------------------------------------------------------------ families
def parse(spec: str):
    return tuple(tuple(sorted(int(ch) - 1 for ch in w)) for w in spec.replace(",", " ").split())


FORB = {  # BCL Table 2 (1-based); read from arXiv:2108.10406v3
    "J4": "123 124 125 134 135 145",
    "K4m": "123 124 134",
    "C5": "123 234 345 145 125",
    # K_5^< : K_5^3 minus two triples meeting in two vertices (BCL Section 3.x: 'the two missing edges intersect in
    # two vertices'); missing 245 and 345 (they share 4, 5)
    "K5lt": "123 124 125 134 135 145 234 235",
    # K_5^{4-}: 4 edges on 5 vertices (4-uniform)
    "K5m4": "1234 1235 1245 1345",
    # families added for the tier 2/3 certificates (BCL arXiv:2108.10406v3, Sections 3.19, 3.21, 3.x; Table 2)
    # K_5^= : K_5^3 minus two triples meeting in exactly ONE vertex (missing 123 and 345, which share 3)
    "K5eq": "124 125 134 135 145 234 235 245",
    # K_5^{3-}: K_5^3 minus one triple (9 edges; missing 123)
    "K5m3": "124 125 134 135 145 234 235 245 345",
    # K_6^3: all 20 triples of [6]
    "K63": " ".join("".join(str(v + 1) for v in e) for e in combinations(range(6), 3)),
    # K_6^{4-}: K_6^4 minus one 4-set (14 edges; missing 3456)
    "K6m4": " ".join("".join(str(v + 1) for v in e) for e in combinations(range(6), 4) if e != (2, 3, 4, 5)),
    # K_6^4 and K_7^4
    "K64": " ".join("".join(str(v + 1) for v in e) for e in combinations(range(6), 4)),
    "K74": " ".join("".join(str(v + 1) for v in e) for e in combinations(range(7), 4)),
}

PROBLEMS = {
    # name: (r, mode, forbidden keys, objective)
    "s3_J4": (3, "sub", ["J4"], "cosig3"),
    "s3_K4m": (3, "sub", ["K4m"], "cosig3"),
    "s3_K4mC5": (3, "sub", ["K4m", "C5"], "cosig3"),
    "s3_K5lt": (3, "sub", ["K5lt"], "cosig3"),
    "c4_K5m": (4, "co", ["K5m4"], "sigma4"),
    "s3_K5eq": (3, "sub", ["K5eq"], "cosig3"),
    "s3_C5": (3, "sub", ["C5"], "cosig3"),
    "s3_K5m": (3, "sub", ["K5m3"], "cosig3"),
    "s3_K6": (3, "sub", ["K63"], "cosig3"),
    "c4_K6m": (4, "co", ["K6m4"], "sigma4"),
    "c4_K6": (4, "co", ["K64"], "sigma4"),
    "c4_K7": (4, "co", ["K74"], "sigma4"),
}


def copies(es, n: int, r: int) -> list:
    """All labelled copies (as n-vertex masks) of the edge set es in K_n^(r) (injective vertex maps)."""
    vs = sorted({v for e in es for v in e})
    out = set()
    for img in permutations(range(n), len(vs)):
        mp = dict(zip(vs, img))
        out.add(sum(1 << cidx(tuple(mp[v] for v in e)) for e in es))
    return sorted(out)


@lru_cache(maxsize=None)
def problem_copies(name: str, n: int) -> tuple:
    r, mode, keys, _ = PROBLEMS[name]
    cs = set()
    for k in keys:
        cs.update(copies(parse(FORB[k]), n, r))
    return tuple(sorted(cs))


def admissible(name: str, mask: int, n: int) -> bool:
    r, mode, keys, _ = PROBLEMS[name]
    mask = int(mask)
    for c in problem_copies(name, n):
        if mode == "sub" and (mask & c) == c:
            return False
        if mode == "co" and (mask & c) == 0:
            return False
    return True


def admissible_np(name: str, masks, n: int) -> np.ndarray:
    r, mode, keys, _ = PROBLEMS[name]
    masks = np.asarray(masks, dtype=np.int64)
    ok = np.ones(masks.shape, dtype=bool)
    for c in problem_copies(name, n):
        c = np.int64(c)
        if mode == "sub":
            ok &= (masks & c) != c
        else:
            ok &= (masks & c) != 0
    return ok


# ------------------------------------------------------------------------------------------------ objectives
def edge_set(mask: int, n: int, r: int) -> set:
    return {frozenset(e) for i, e in enumerate(edges(n, r)) if (int(mask) >> i) & 1}


def codeg(mask: int, n: int, r: int) -> dict:
    E = edge_set(mask, n, r)
    out = {}
    for T in combinations(range(n), r - 1):
        out[T] = sum(1 for v in range(n) if v not in T and frozenset(T + (v,)) in E)
    return out


def co2(mask: int, n: int, r: int) -> int:
    return sum(c * c for c in codeg(mask, n, r).values())


def sigma_eq1(mask: int, n: int, r: int) -> Fraction:
    """BCL eq. (1) normalisation: co2 / (C(n, r-1) (n-r+1)^2)."""
    return Fraction(co2(mask, n, r), comb(n, r - 1) * (n - r + 1) ** 2)


def sigma_xy(mask: int, n: int, r: int) -> Fraction:
    """x != y form: sum_T c(c-1) / (C(n,r-1)(n-r+1)(n-r))."""
    return Fraction(sum(c * (c - 1) for c in codeg(mask, n, r).values()), comb(n, r - 1) * (n - r + 1) * (n - r))


def cosig3_value(mask: int, n: int) -> Fraction:
    """1 - sigma_xy, computed from the explicit (T, x, y) loop (pairs T, ordered x != y outside T)."""
    E = edge_set(mask, n, 3)
    good = tot = 0
    for T in combinations(range(n), 2):
        out = [v for v in range(n) if v not in T]
        for x in out:
            for y in out:
                if x != y:
                    tot += 1
                    good += (frozenset(T + (x,)) in E) and (frozenset(T + (y,)) in E)
    return 1 - Fraction(good, tot)


def sigma4_value(mask: int, n: int) -> Fraction:
    """kappa(H) = 2 d(H) - gamma(H), from the explicit (T, x, y) loop (3-sets T, ordered x != y outside T)."""
    E = edge_set(mask, n, 4)
    d = Fraction(len(E), comb(n, 4))
    good = tot = 0
    for T in combinations(range(n), 3):
        out = [v for v in range(n) if v not in T]
        for x in out:
            for y in out:
                if x != y:
                    tot += 1
                    good += (frozenset(T + (x,)) in E) and (frozenset(T + (y,)) in E)
    return 2 * d - Fraction(good, tot)


def objective_value(obj: str, mask: int, n: int) -> Fraction:
    return cosig3_value(mask, n) if obj == "cosig3" else sigma4_value(mask, n)


# 7-vertex integer numerators: cosig3 = F/210 with F = 210 - sum_Q C(k_Q, 2); sigma4 = F/420 with F = 24 e - P
DEN7 = {"cosig3": 210, "sigma4": 420}


def num7(obj: str, mask: int) -> int:
    v = objective_value(obj, mask, 7) * DEN7[obj]
    assert v.denominator == 1
    return int(v)


def rand_mask(n: int, r: int, p: float, rng: random.Random) -> int:
    return sum(1 << i for i in range(comb(n, r)) if rng.random() < p)
