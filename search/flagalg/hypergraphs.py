"""r-graphs on labelled vertex sets as bitmasks; canonical forms; isomorphism-class enumeration.

Conventions
-----------
* Vertices are 0..n-1. The r-subsets of [n] are ordered in **colex** order (sorted by the reversed tuple), so the
  r-subsets of [n-1] are exactly the first C(n-1, r) of [n], and the r-subsets inside the first s vertices are the
  first C(s, r) bits of any mask.  An r-graph is the integer  sum_{e in E} 2^{colex_index(e)}.
* Canonical form of a mask under a permutation group G = min over g in G of the permuted mask (as an integer).
  Brute force over G, vectorised as a float64 matrix product (exact while C(n, r) <= 52 bits).
* Admissibility is a *hereditary* predicate (closed under induced subgraphs), evaluated on batches of masks.
  The one used here: "every p-subset of the vertices contains an edge" (complement form of the Turan problem).
"""
from __future__ import annotations

from functools import lru_cache
from itertools import combinations, permutations
from math import comb

import numpy as np

MAX_BITS = 52  # float64 exactness limit for the matrix-product canonical form


# ----------------------------------------------------------------------------------------------- edges, indices
@lru_cache(maxsize=None)
def edges(n: int, r: int) -> tuple:
    """All r-subsets of range(n) in colex order (tuple of sorted tuples)."""
    return tuple(sorted(combinations(range(n), r), key=lambda e: e[::-1]))


def colex_index(e) -> int:
    """Colex rank of a sorted r-tuple."""
    return sum(comb(v, i + 1) for i, v in enumerate(e))


@lru_cache(maxsize=None)
def edge_pos_arrays(m: int, r: int) -> tuple:
    """For the r-subsets of range(m) in colex order: r position arrays (j-th smallest element of every edge)."""
    E = edges(m, r)
    return tuple(np.array([e[j] for e in E], dtype=np.int64) for j in range(r))


def popcount(x: int) -> int:
    return bin(int(x)).count("1")


def mask_to_edges(mask: int, n: int, r: int) -> list:
    return [e for i, e in enumerate(edges(n, r)) if (int(mask) >> i) & 1]


def edges_to_mask(es, n: int, r: int) -> int:
    return sum(1 << colex_index(tuple(sorted(e))) for e in es)


def bits_matrix(masks: np.ndarray, nbits: int) -> np.ndarray:
    """(G,) int64 masks -> (G, nbits) float64 0/1 matrix."""
    masks = np.asarray(masks, dtype=np.int64)
    return ((masks[:, None] >> np.arange(nbits, dtype=np.int64)[None, :]) & 1).astype(np.float64)


# ----------------------------------------------------------------------------------------------- permutation groups
@lru_cache(maxsize=None)
def group_perms(n: int, kind: str = "all", s: int = 0) -> tuple:
    """Permutation groups on range(n), as tuples of images.

    kind = 'all'      : the symmetric group S_n;
    kind = 'fix'      : permutations fixing 0..s-1 pointwise (label-preserving maps of sigma-flags);
    kind = 'setwise'  : permutations mapping {0..s-1} onto itself (S_s x S_{n-s}; unordered roots).
    """
    if kind == "all":
        return tuple(permutations(range(n)))
    if kind == "fix":
        return tuple(tuple(range(s)) + p for p in permutations(range(s, n)))
    if kind == "setwise":
        return tuple(a + b for a in permutations(range(s)) for b in permutations(range(s, n)))
    raise ValueError(kind)


@lru_cache(maxsize=None)
def perm_weight_matrix(n: int, r: int, kind: str = "all", s: int = 0) -> np.ndarray:
    """W[g, e] = 2^{index of image of edge e under g}; permuted mask = bits @ W.T (exact in float64)."""
    E = edges(n, r)
    if len(E) > MAX_BITS:
        raise ValueError(f"C({n},{r}) = {len(E)} > {MAX_BITS}: float canonical form not exact")
    perms = group_perms(n, kind, s)
    W = np.empty((len(perms), len(E)), dtype=np.float64)
    for gi, g in enumerate(perms):
        for ei, e in enumerate(E):
            W[gi, ei] = float(1 << colex_index(tuple(sorted(g[v] for v in e))))
    return W


def canon_batch(masks, n: int, r: int, kind: str = "all", s: int = 0, batch: int = 4096) -> np.ndarray:
    """Canonical forms (min over the group of permuted masks) of a batch of masks."""
    masks = np.asarray(masks, dtype=np.int64)
    nb = comb(n, r)
    if nb == 0 or masks.size == 0:
        return masks.copy()
    W = perm_weight_matrix(n, r, kind, s)
    out = np.empty_like(masks)
    for i in range(0, masks.size, batch):
        B = bits_matrix(masks[i:i + batch], nb)
        out[i:i + batch] = (B @ W.T).min(axis=1).astype(np.int64)
    return out


# ----------------------------------------------------------------------------------------------- admissibility
@lru_cache(maxsize=None)
def pset_masks(n: int, r: int, p: int) -> np.ndarray:
    """(E, #p-sets) 0/1 float matrix: column = r-subsets contained in that p-subset of range(n)."""
    E = edges(n, r)
    ps = list(combinations(range(n), p))
    M = np.zeros((len(E), len(ps)), dtype=np.float64)
    for j, P in enumerate(ps):
        for e in combinations(P, r):
            M[colex_index(e), j] = 1.0
    return M


class EveryPSetHasEdge:
    """Hereditary predicate: every p-subset of the vertex set contains at least one edge.

    For r-graphs this is the complement form of the Turan problem for K_p^{(r)}: H is admissible iff its complement
    is K_p^{(r)}-free.
    """

    def __init__(self, p: int):
        self.p = p
        self.name = f"every {p}-set contains an edge"

    def __call__(self, masks, n: int, r: int) -> np.ndarray:
        masks = np.asarray(masks, dtype=np.int64)
        if n < self.p or masks.size == 0:
            return np.ones(masks.shape, dtype=bool)
        B = bits_matrix(masks, comb(n, r))
        return (B @ pset_masks(n, r, self.p) > 0.5).all(axis=1)


class PSetEdgeCounts:
    """Hereditary predicate: every p-subset spans a number of edges in `allowed`.

    Example (Razborov 2010, complement form): r = 3, p = 4, allowed = {1, 2, 4} <=> the complement is K_4^{(3)}-free
    (no 4-set with 0 edges here) and has no induced E_1 = 4 vertices spanning exactly one edge (no 4-set with 3 edges
    here).
    """

    def __init__(self, p: int, allowed):
        self.p = p
        self.allowed = sorted(set(allowed))
        self.name = f"every {p}-set spans a number of edges in {self.allowed}"

    def __call__(self, masks, n: int, r: int) -> np.ndarray:
        masks = np.asarray(masks, dtype=np.int64)
        if n < self.p or masks.size == 0:
            return np.ones(masks.shape, dtype=bool)
        B = bits_matrix(masks, comb(n, r))
        cnt = np.rint(B @ pset_masks(n, r, self.p)).astype(np.int64)
        return np.isin(cnt, self.allowed).all(axis=1)


class NoPredicate:
    name = "all r-graphs"

    def __call__(self, masks, n, r):
        return np.ones(np.asarray(masks).shape, dtype=bool)


# ----------------------------------------------------------------------------------------------- enumeration
def enumerate_classes(N: int, r: int, pred, chunk: int = 1 << 16, verbose: bool = False) -> dict:
    """Canonical representatives of the isomorphism classes of admissible r-graphs on n vertices, n = 0..N.

    Vertex augmentation: every admissible graph on n vertices, minus its last vertex, is admissible (hereditary
    predicate), hence isomorphic to a representative R on n-1 vertices; so every class on n vertices contains some
    R + (new vertex with an arbitrary link).  Candidates are filtered by the predicate and deduplicated by canonical
    form.  Returns {n: sorted int64 array of canonical masks}.
    """
    reps = {0: np.array([0], dtype=np.int64)}
    for n in range(1, N + 1):
        old = comb(n - 1, r)
        new = comb(n - 1, r - 1)
        links = np.arange(1 << new, dtype=np.int64) << old
        found = []
        prev = reps[n - 1]
        # process in chunks of candidates
        per = max(1, chunk // max(1, links.size))
        for i in range(0, prev.size, per):
            cand = (prev[i:i + per, None] | links[None, :]).ravel()
            cand = cand[pred(cand, n, r)]
            if cand.size:
                found.append(np.unique(canon_batch(cand, n, r)))
        reps[n] = np.unique(np.concatenate(found)) if found else np.zeros(0, dtype=np.int64)
        if verbose:
            print(f"  n={n}: {reps[n].size} classes", flush=True)
    return reps


def adjacency_flat(masks: np.ndarray, n: int, r: int) -> np.ndarray:
    """(G, n**r) bool: entry [g, sum_j v_j n^(r-1-j)] = 1 iff {v_0..v_{r-1}} (distinct) is an edge of graph g."""
    masks = np.asarray(masks, dtype=np.int64)
    A = np.zeros((masks.size, n ** r), dtype=bool)
    E = edges(n, r)
    for ei, e in enumerate(E):
        bit = ((masks >> ei) & 1).astype(bool)
        for pe in permutations(e):
            idx = 0
            for v in pe:
                idx = idx * n + v
            A[:, idx] = bit
    return A


def induced_mask_index(tuples: np.ndarray, n: int, r: int) -> np.ndarray:
    """For (T, m) vertex tuples of a graph on n vertices: (T, C(m, r)) flat adjacency indices of the induced
    labelled subgraph's possible edges (colex order on positions)."""
    tuples = np.asarray(tuples, dtype=np.int64)
    m = tuples.shape[1]
    P = edge_pos_arrays(m, r)
    idx = np.zeros((tuples.shape[0], len(P[0])), dtype=np.int64)
    for j in range(r):
        idx = idx * n + tuples[:, P[j]]
    return idx
