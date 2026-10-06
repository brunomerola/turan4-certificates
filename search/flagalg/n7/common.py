"""Tables for the N = 7, r = 4 flag-algebra problems t(p, 4), p in {5, 6, 7} (complement form: minimise the edge
density of 4-graphs in which every p-set contains an edge).

Conventions (same as ../hypergraphs.py): vertices 0..6, 4-subsets in colex order, a 4-graph on 7 vertices is a
35-bit integer.  The 4-subsets of {0..5} are the first 15 bits, so a one-vertex extension of a 6-vertex graph R by
vertex 6 with link L (a set of 3-subsets of {0..5}, 20 bits in colex order) is  R | (L << 15).

Blocks.  A block is a pair (s, m) with m > s, 2m - s <= 7 and m >= 4 (smaller m give 1x1 Gram matrices).  At N = 7:
    (1,4) (2,4) (3,4) (3,5) (4,5) (5,6).
For a type sigma (canonical admissible 4-graph on roots 0..s-1) the sigma-flags on m vertices are those of
../flags.py (flag_family).  The flag-pair density matrix M_sigma(H) is averaged over all configurations
(theta, U1, U2): theta an injective s-tuple, (U1, U2) an ordered pair of disjoint (m-s)-sets; their number is
denom(s, m) and DEN = 5040 = lcm of all denominators.

Symmetric evaluation (the idea of Jeong-Park-Im-Lee-Yang 2026, arXiv 2609.27495, github.com/taeyool/tetrahedron-turan,
Apache-2.0, "canonical root labelling followed by the full type automorphism action"; independent re-implementation,
no code copied):  for an unordered root set S whose induced graph is isomorphic to sigma, the orderings theta of S
with H[theta] = sigma are theta0 o h, h in Aut(sigma), for one canonical theta0.  Ordering theta0 o h relabels the
flags by the permutation pi_h of the flag list, so
    <Q, M_sigma(H)> = sum over S, unordered {U1, U2}:  (2 / denom) * sum_{h in Aut(sigma)} Q[pi_h a, pi_h b],
where a, b are the flags of (theta0, U1), (theta0, U2).  We store Qsum = sum_h Q[pi_h, pi_h] and integer term weights
wint = 2 * DEN / denom, so that  <Q, M(H)> = (1/DEN) * sum_terms wint * Qsum[a, b].  This is exact for every Q (no
invariance of Q is assumed).
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, permutations
from math import comb, factorial

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
if PARENT not in sys.path:
    sys.path.insert(0, PARENT)

from flags import flag_family  # noqa: E402
from hypergraphs import EveryPSetHasEdge, canon_batch, colex_index, edges, enumerate_classes  # noqa: E402

N = 7
R = 4
NE = comb(N, R)            # 35 possible edges
DEN = 5040                 # common denominator of all block densities at N = 7
BLOCK_SM = ((1, 4), (2, 4), (3, 4), (3, 5), (4, 5), (5, 6))
E7 = edges(N, R)
EIDX = {e: i for i, e in enumerate(E7)}
DATA = os.path.join(HERE, "data")
RESULTS = os.path.join(HERE, "results")
LOGS = os.path.join(PARENT, "logs")


def apply_perm(mask: int, n: int, g, r: int = R) -> int:
    """Image of a labelled r-graph on range(n) under the vertex map v -> g[v]."""
    out = 0
    for i, e in enumerate(edges(n, r)):
        if (mask >> i) & 1:
            out |= 1 << colex_index(tuple(sorted(g[v] for v in e)))
    return out


def n_configurations(s: int, m: int, n: int = N) -> int:
    return factorial(n) // factorial(n - s) * comb(n - s, m - s) * comb(n - m, m - s)


@lru_cache(maxsize=None)
def pset_masks7(p: int) -> np.ndarray:
    """35-bit masks of the 4-subsets inside each p-subset of {0..6} that contains vertex 6 (the p-sets that a
    one-vertex extension of an admissible 6-vertex graph must also hit)."""
    out = []
    for P in combinations(range(N), p):
        if N - 1 not in P:
            continue
        out.append(sum(1 << EIDX[e] for e in combinations(P, R)))
    return np.array(out, dtype=np.int64)


@lru_cache(maxsize=None)
def all_pset_masks7(p: int) -> np.ndarray:
    return np.array([sum(1 << EIDX[e] for e in combinations(P, R)) for P in combinations(range(N), p)],
                    dtype=np.int64)


def admissible7(masks, p: int) -> np.ndarray:
    masks = np.asarray(masks, dtype=np.int64)
    ok = np.ones(masks.shape, dtype=bool)
    for pm in all_pset_masks7(p):
        ok &= (masks & pm) != 0
    return ok


@lru_cache(maxsize=None)
def reps6(p: int) -> np.ndarray:
    """Canonical representatives (engine's brute-force canonical form) of admissible 4-graphs on 6 vertices."""
    return enumerate_classes(6, R, EveryPSetHasEdge(p))[6]


# ----------------------------------------------------------------------------------------------- blocks
@dataclass
class N7Block:
    s: int
    m: int
    types: list            # canonical sigma masks (types with at least one flag)
    flags: list            # per type: sorted int64 array of flag masks (m-vertex, roots 0..s-1)
    aut: list              # per type: (|Aut sigma|, n_t) int32 array, row h = pi_h on flag indices
    aut_roots: list        # per type: list of root permutations h (tuples)
    denom: int             # number of configurations (theta, U1, U2)
    wint: int              # 2 * DEN / denom
    nsig: int              # C(s, 4)
    nfree: int             # C(m, 4) - C(s, 4)
    type_of: np.ndarray    # (2^nsig,) labelled sigma code -> type index, -1 if not admissible
    flag_of: np.ndarray    # (2^nsig, 2^nfree) -> flag index in canonical root labelling, -1 if not admissible
    rootsets: np.ndarray   # (nrs, s) sorted root sets
    sigpos: np.ndarray     # (nrs, nsig) edge indices (in H) of the 4-subsets of the root set, colex
    freepos: np.ndarray    # (nrs, nU, nfree) edge indices of the free 4-subsets of root set + U
    pairs: np.ndarray      # (npairs, 2) disjoint unordered pairs of U indices

    @property
    def dims(self):
        return [len(f) for f in self.flags]


def _flag_index(flags: np.ndarray, mask: int, s: int, m: int) -> int:
    c = int(canon_batch(np.array([mask], dtype=np.int64), m, R, kind="fix", s=s)[0])
    i = int(np.searchsorted(flags, c))
    return i if i < len(flags) and int(flags[i]) == c else -1


def build_block(s: int, m: int, p: int) -> N7Block:
    pred = EveryPSetHasEdge(p)
    nsig = comb(s, R)
    nfree = comb(m, R) - nsig
    classes = enumerate_classes(s, R, pred)[s]
    types, flags, auts, aut_roots = [], [], [], []
    for sig in classes:
        fam = flag_family(s, m, int(sig), R, pred)
        if fam.n == 0:
            continue
        fl = np.asarray(fam.flags, dtype=np.int64)
        hs = [h for h in permutations(range(s)) if apply_perm(int(sig), s, h) == int(sig)]
        P = np.empty((len(hs), fl.size), dtype=np.int32)
        for hi, h in enumerate(hs):
            g = tuple(h) + tuple(range(s, m))
            for a, f in enumerate(fl):
                j = _flag_index(fl, apply_perm(int(f), m, g), s, m)
                assert j >= 0
                P[hi, a] = j
            assert sorted(P[hi].tolist()) == list(range(fl.size))
        types.append(int(sig))
        flags.append(fl)
        auts.append(P)
        aut_roots.append(hs)
    tindex = {t: i for i, t in enumerate(types)}
    type_of = np.full(1 << nsig, -1, dtype=np.int32)
    flag_of = np.full((1 << nsig, 1 << nfree), -1, dtype=np.int32)
    perms_s = list(permutations(range(s)))
    for sc in range(1 << nsig):
        c = int(canon_batch(np.array([sc], dtype=np.int64), s, R, kind="all")[0]) if nsig else 0
        if c not in tindex:
            continue
        t = tindex[c]
        type_of[sc] = t
        g0 = next(h for h in perms_s if apply_perm(sc, s, h) == c)
        g = tuple(g0) + tuple(range(s, m))
        for fc in range(1 << nfree):
            flag_of[sc, fc] = _flag_index(flags[t], apply_perm(sc | (fc << nsig), m, g), s, m)
    # geometry
    rootsets = list(combinations(range(N), s))
    Esm = edges(m, R)
    sig_sub = Esm[:nsig]
    free_sub = Esm[nsig:]
    nU = None
    sigpos, freepos = [], []
    pairs = None
    for rs in rootsets:
        O = [v for v in range(N) if v not in rs]
        Us = list(combinations(O, m - s))
        if nU is None:
            nU = len(Us)
            idx = list(combinations(range(len(O)), m - s))
            pairs = [(i, j) for i in range(nU) for j in range(i + 1, nU) if not set(idx[i]) & set(idx[j])]
        sigpos.append([EIDX[tuple(sorted(rs[x] for x in e))] for e in sig_sub])
        fp = []
        for U in Us:
            V = tuple(rs) + tuple(U)
            fp.append([EIDX[tuple(sorted(V[x] for x in e))] for e in free_sub])
        freepos.append(fp)
    denom = n_configurations(s, m)
    assert (2 * DEN) % denom == 0
    assert len(rootsets) * factorial(s) * 2 * len(pairs) == denom
    return N7Block(s, m, types, flags, auts, aut_roots, denom, 2 * DEN // denom, nsig, nfree, type_of, flag_of,
                   np.array(rootsets, dtype=np.int64).reshape(len(rootsets), s),
                   np.array(sigpos, dtype=np.int64).reshape(len(rootsets), nsig),
                   np.array(freepos, dtype=np.int64).reshape(len(rootsets), nU, nfree),
                   np.array(pairs, dtype=np.int64).reshape(-1, 2))


@dataclass
class N7Model:
    """All blocks for one p, plus the flat padded tables consumed by the numba kernels (kernels.py)."""
    p: int
    blocks: list
    keys: list                 # (block index, type index) for every Gram matrix, in order
    qoff: np.ndarray           # (nb, maxT) offset of Gram (bi, t) in the flat Q vector, -1 if absent
    qdim: np.ndarray           # (nb, maxT)
    qsize: int
    nrs: np.ndarray
    nsig: np.ndarray
    nU: np.ndarray
    nfree: np.ndarray
    npairs: np.ndarray
    wint: np.ndarray
    sigpos: np.ndarray         # (nb, 35, 5)
    freepos: np.ndarray        # (nb, 35, 20, 10)
    pairs: np.ndarray          # (nb, 15, 2)
    type_of: np.ndarray        # (nb, 32)
    flag_of: np.ndarray        # (nb, 32, 1024)
    psets: np.ndarray          # 35-bit masks of all p-subsets (for admissibility checks)
    psets_new: np.ndarray      # p-subsets containing vertex 6
    cdd_pos: np.ndarray        # (70, 2) edge index pairs {v}+U1, {v}+U2 for the stationarity term

    def gram_dims(self):
        return [int(self.qdim[b, t]) for b, t in self.keys]


def build_model(p: int, sm=BLOCK_SM) -> N7Model:
    blocks = [build_block(s, m, p) for s, m in sm]
    nb = len(blocks)
    maxT = max(len(b.types) for b in blocks)
    qoff = np.full((nb, maxT), -1, dtype=np.int64)
    qdim = np.zeros((nb, maxT), dtype=np.int64)
    keys = []
    off = 0
    for bi, b in enumerate(blocks):
        for t, fl in enumerate(b.flags):
            qoff[bi, t] = off
            qdim[bi, t] = len(fl)
            off += len(fl) ** 2
            keys.append((bi, t))
    sigpos = np.zeros((nb, 35, 5), dtype=np.int64)
    freepos = np.zeros((nb, 35, 20, 10), dtype=np.int64)
    pairs = np.zeros((nb, 15, 2), dtype=np.int64)
    type_of = np.full((nb, 32), -1, dtype=np.int32)
    flag_of = np.full((nb, 32, 1024), -1, dtype=np.int32)
    nrs, nsig, nU, nfree, npairs, wint = (np.zeros(nb, dtype=np.int64) for _ in range(6))
    for bi, b in enumerate(blocks):
        nrs[bi] = b.rootsets.shape[0]
        nsig[bi] = b.nsig
        nU[bi] = b.freepos.shape[1]
        nfree[bi] = b.nfree
        npairs[bi] = b.pairs.shape[0]
        wint[bi] = b.wint
        sigpos[bi, :nrs[bi], :nsig[bi]] = b.sigpos
        freepos[bi, :nrs[bi], :nU[bi], :nfree[bi]] = b.freepos
        pairs[bi, :npairs[bi]] = b.pairs
        type_of[bi, :b.type_of.size] = b.type_of
        flag_of[bi, :b.flag_of.shape[0], :b.flag_of.shape[1]] = b.flag_of
    cdd = []
    for v in range(N):
        O = [u for u in range(N) if u != v]
        for U in combinations(O, 3):
            Uc = tuple(u for u in O if u not in U)
            if U < Uc:
                cdd.append((EIDX[tuple(sorted((v,) + U))], EIDX[tuple(sorted((v,) + Uc))]))
    return N7Model(p, blocks, keys, qoff, qdim, off, nrs, nsig, nU, nfree, npairs, wint, sigpos, freepos, pairs,
                   type_of, flag_of, all_pset_masks7(p), pset_masks7(p), np.array(cdd, dtype=np.int64))


def summary(model: N7Model) -> list:
    out = []
    for bi, b in enumerate(model.blocks):
        out.append({"s": b.s, "m": b.m, "denom": b.denom, "wint": b.wint, "types": b.types, "dims": b.dims,
                    "aut_orders": [len(a) for a in b.aut]})
    return out


if __name__ == "__main__":
    import json
    import time
    for p in (5, 6, 7):
        t0 = time.time()
        M = build_model(p)
        print(f"p={p}: build {time.time() - t0:.1f}s, Gram dims {M.gram_dims()}, flat Q size {M.qsize}")
        print(json.dumps(summary(M)))
