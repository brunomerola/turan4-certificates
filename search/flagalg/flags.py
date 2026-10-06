"""Types, flags and flag-product densities (Razborov's flag algebras) for hereditary r-graph problems.

Definitions (standard, Razborov 2007):
* A type sigma is an admissible r-graph on the labelled vertex set {0..s-1}.  One labelled representative per
  isomorphism class (the canonical representative), as usual.
* A sigma-flag on m vertices is an admissible r-graph F on {0..m-1} whose restriction to {0..s-1} *equals* sigma
  (labels = vertices 0..s-1), up to isomorphisms fixing 0..s-1 pointwise.
* For an admissible H on N >= 2m - s vertices, p(F_a, F_b; H) is the probability that, for a uniformly random
  injection theta: [s] -> V(H) and a uniformly random ordered pair (U1, U2) of disjoint (m-s)-subsets of
  V(H) minus im(theta), both (H[theta + U1], theta) ~ F_a and (H[theta + U2], theta) ~ F_b (as labelled flags;
  this includes the event H[im theta] = sigma under the labelling theta).
  M_sigma(H) = (p(F_a, F_b; H))_{a,b} is symmetric.
* Basic inequality: for every PSD Q and every large admissible G,
      sum_H p(H; G) <Q, M_sigma(H)> = E_theta[ v_theta^T Q v_theta ] + O(1/n) >= -O(1/n),
  where v_theta = (p(F_a; G, theta))_a.  Hence for every admissible G,  d(G) >= min_H (d(H) - sum <Q, M(H)>) - o(1).

All densities are kept as exact integer counts with an integer denominator per (N, s, m).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations, permutations
from math import comb, factorial

import numpy as np

from hypergraphs import (adjacency_flat, canon_batch, edges, enumerate_classes, induced_mask_index, popcount)


@dataclass
class FlagFamily:
    """All sigma-flags on m vertices for one type sigma."""
    s: int
    m: int
    sigma: int                      # labelled mask of sigma on {0..s-1}
    flags: np.ndarray               # canonical masks (sorted)
    lookup: np.ndarray              # free-bits -> flag index (-1 if not admissible)

    @property
    def n(self) -> int:
        return int(self.flags.size)


def flag_family(s: int, m: int, sigma: int, r: int, pred) -> FlagFamily:
    """Enumerate the sigma-flags on m vertices: masks sigma | free << C(s, r) (colex: edges inside [s] come first),
    admissible, canonical under permutations fixing 0..s-1 pointwise."""
    low = comb(s, r)
    nfree = comb(m, r) - low
    free = np.arange(1 << nfree, dtype=np.int64)
    masks = np.int64(sigma) | (free << low)
    ok = pred(masks, m, r)
    lookup = np.full(free.size, -1, dtype=np.int64)
    can = canon_batch(masks[ok], m, r, kind="fix", s=s)
    flags, inv = np.unique(can, return_inverse=True)
    lookup[np.nonzero(ok)[0]] = inv
    return FlagFamily(s, m, int(sigma), flags, lookup)


def configurations(N: int, s: int, m: int):
    """All (theta, U1, U2): theta an ordered s-tuple of distinct vertices, U1, U2 disjoint (m-s)-subsets of the rest
    (ordered pair).  Returns (T1, T2) as (K, m) int arrays: T1 = theta + U1, T2 = theta + U2."""
    T1, T2 = [], []
    for theta in permutations(range(N), s):
        rest = [v for v in range(N) if v not in theta]
        for U1 in combinations(rest, m - s):
            rest2 = [v for v in rest if v not in U1]
            for U2 in combinations(rest2, m - s):
                T1.append(theta + U1)
                T2.append(theta + U2)
    return np.array(T1, dtype=np.int64).reshape(-1, m), np.array(T2, dtype=np.int64).reshape(-1, m)


def n_configurations(N: int, s: int, m: int) -> int:
    return (factorial(N) // factorial(N - s)) * comb(N - s, m - s) * comb(N - m, m - s)


@dataclass
class Block:
    """One SDP block: type sigma (s vertices), flags on m vertices, and the sparse matrices M_sigma(H)."""
    fam: FlagFamily
    denom: int                       # common denominator of all entries
    h: np.ndarray                    # COO: graph index
    a: np.ndarray                    # COO: flag index a
    b: np.ndarray                    # COO: flag index b
    cnt: np.ndarray                  # COO: integer count (entry = cnt / denom); symmetric (both (a,b) and (b,a))

    @property
    def size(self) -> int:
        return self.fam.n


def block_densities(H: np.ndarray, N: int, r: int, fam: FlagFamily, batch: int = 0) -> Block:
    """Exact counts for M_sigma(H) for all H (canonical masks on N vertices)."""
    s, m = fam.s, fam.m
    assert 2 * m - s <= N, (s, m, N)
    T1, T2 = configurations(N, s, m)
    if batch <= 0:  # keep the (G, K, C(m, r)) int64 temporaries around 100 MB
        batch = max(1, int(1.5e7 // max(1, T1.shape[0] * comb(m, r))))
    I1 = induced_mask_index(T1, N, r)
    I2 = induced_mask_index(T2, N, r)
    nb = comb(m, r)
    low = comb(s, r)
    lowmask = (1 << low) - 1
    pw = (np.int64(1) << np.arange(nb, dtype=np.int64))
    n = fam.n
    hs, codes_all = [], []
    for i in range(0, H.size, batch):
        A = adjacency_flat(H[i:i + batch], N, r)                          # (G, N^r)
        M1 = (A[:, I1].astype(np.int64) * pw).sum(axis=2)                # (G, K)
        M2 = (A[:, I2].astype(np.int64) * pw).sum(axis=2)
        valid = (M1 & lowmask) == fam.sigma                              # theta induces sigma (same for M2)
        a = fam.lookup[M1 >> low]
        b = fam.lookup[M2 >> low]
        g, k = np.nonzero(valid)
        aa, bb = a[g, k], b[g, k]
        if (aa < 0).any() or (bb < 0).any():
            raise RuntimeError("induced flag not admissible: predicate is not hereditary?")
        codes = (g + i) * (n * n) + aa * n + bb
        codes_all.append(codes)
    codes = np.concatenate(codes_all) if codes_all else np.zeros(0, dtype=np.int64)
    u, c = np.unique(codes, return_counts=True)
    h = u // (n * n)
    ab = u % (n * n)
    return Block(fam, n_configurations(N, s, m), h, ab // n, ab % n, c.astype(np.int64))


def standard_type_sizes(N: int, mode: str = "standard"):
    """(s, m) pairs. 'standard': 2m - s = N, m > s (Flagmatic default).  'all': every 2m - s <= N, m > s.
    'none': no SOS blocks (the bound is then min_H d(H), a sanity check)."""
    out = []
    if mode == "none":
        return out
    for s in range(0, N):
        for m in range(s + 1, N + 1):
            if (mode == "standard" and 2 * m - s == N) or (mode == "all" and 2 * m - s <= N):
                out.append((s, m))
    return out


@dataclass
class FlagProblem:
    N: int
    r: int
    pred: object
    H: np.ndarray                    # canonical masks of admissible r-graphs on N vertices
    edges_count: np.ndarray          # e(H)
    blocks: list = field(default_factory=list)
    classes: dict = field(default_factory=dict)

    @property
    def nH(self) -> int:
        return int(self.H.size)

    def density_num_den(self):
        """Edge density d(H) = e(H) / C(N, r) as (numerators, denominator)."""
        return self.edges_count, comb(self.N, self.r)


def build_problem(N: int, r: int, pred, type_mode: str = "standard", verbose: bool = False) -> FlagProblem:
    classes = enumerate_classes(N, r, pred, verbose=verbose)
    H = classes[N]
    ec = np.array([popcount(x) for x in H], dtype=np.int64)
    prob = FlagProblem(N, r, pred, H, ec, classes=classes)
    for s, m in standard_type_sizes(N, type_mode):
        for sigma in classes[s]:
            fam = flag_family(s, m, int(sigma), r, pred)
            if fam.n == 0:
                continue
            blk = block_densities(H, N, r, fam)
            prob.blocks.append(blk)
            if verbose:
                print(f"  type s={s} sigma={int(sigma)} m={m}: {fam.n} flags, {blk.cnt.size} nonzeros", flush=True)
    return prob
