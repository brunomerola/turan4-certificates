"""Generalised predicate for the N = 7, r = 4 machinery of ../n7: "every p-set spans at least LAM edges".

Problem: t_LAM(p, 4) = lim min{ e(G)/C(n,4) : every p-set of V(G) spans >= LAM edges }.
For p = 5, LAM = 2 this is the complement form of the Turan density of H^4_4 = K_5^{4-} (the 4-graph on 5 vertices
with 4 edges):  G is admissible iff its edge-complement has at most 3 edges in every 5-set, i.e. is K_5^{4-}-free,
so  pi(K_5^{4-}) = 1 - t_2(5, 4).  LAM = 1 is the original K_p^(4) problem (../n7, unchanged).

How: the n7 modules take their predicate from the name EveryPSetHasEdge (in ../hypergraphs.py and n7/common.py);
configure(p, lam) rebinds that name in both modules to a threshold predicate (PSetEdgeCounts with allowed counts
lam..C(p,4), a hereditary property), BEFORE any n7 table is built.  Everything predicate-dependent is then rebuilt
from it: 6-vertex representatives (common.reps6), types and flags (common.build_block -> flags.flag_family, hence
type_of / flag_of, which also make the kernels flag inadmissible configurations as 'bad'), and the N = 6 seed.
The only hard-coded "contains an edge" tests in n7 are the one-vertex-extension filters of the raw scans
(kernels.scan_raw_*: (h & pset) == 0); this module provides threshold versions (popcount(h & pset) < lam).
No n7 file is modified.  Imports numba/numpy only (no highspy, no ortools).
"""
from __future__ import annotations

import os
import sys
from math import comb

import numpy as np
from numba import njit, prange

HERE = os.path.dirname(os.path.abspath(__file__))
FLAGALG = os.path.dirname(HERE)
N7DIR = os.path.join(FLAGALG, "n7")
for _d in (N7DIR, FLAGALG):
    if _d not in sys.path:
        sys.path.insert(0, _d)

import hypergraphs  # noqa: E402
from hypergraphs import PSetEdgeCounts  # noqa: E402

DATA = os.path.join(HERE, "data")
RESULTS = os.path.join(HERE, "results")
LOGS = os.path.join(FLAGALG, "logs")
CONFIG = {"p": None, "lam": None}
INT64_MIN = np.iinfo(np.int64).min


def make_pred_class(lam: int):
    class AtLeastLam(PSetEdgeCounts):
        def __init__(self, p: int):
            super().__init__(p, range(lam, comb(p, 4) + 1))
            self.name = f"every {p}-set spans at least {lam} edges"
    return AtLeastLam


def configure(p: int, lam: int):
    """Rebind the predicate name in ../hypergraphs.py and n7/common.py (must run before any table is built)."""
    import common
    if CONFIG["lam"] is not None and (CONFIG["p"], CONFIG["lam"]) != (p, lam):
        raise RuntimeError("configure() called twice with different parameters")
    cls = make_pred_class(lam)
    hypergraphs.EveryPSetHasEdge = cls
    common.EveryPSetHasEdge = cls
    common.reps6.cache_clear()
    CONFIG.update(p=p, lam=lam)
    return cls


def admissible7_lam(masks, p: int, lam: int) -> np.ndarray:
    import common
    masks = np.asarray(masks, dtype=np.int64)
    ok = np.ones(masks.shape, dtype=bool)
    for pm in common.all_pset_masks7(p):
        ok &= np.bitwise_count(masks & pm) >= lam
    return ok


def data_dir(lam: int) -> str:
    return os.path.join(DATA, f"l{lam}")


def classes_path(p: int, lam: int) -> str:
    """Same file name as n7 (classes_p{p}.npy) inside data/l{lam}/, so n7's drivers can be pointed at it."""
    return os.path.join(data_dir(lam), f"classes_p{p}.npy")


# ------------------------------------------------------------------------------- raw-scan kernels (threshold)
def make_raw_kernels(lam: int):
    from kernels import _eval_Z, _merge, _unpack, popcnt64
    LAMC = np.int64(lam)

    @njit(parallel=True, cache=False)
    def scan_raw_exact_lam(reps, psets_new, nb, nrs, nsig, nU, nfree, npairs, wint, sigpos, freepos, pairs, type_of,
                           flag_of, qoff, qdim, aux, Qi, lchunk):
        nper = (1 << 20) // lchunk
        ntask = reps.size * nper
        mz = np.full((ntask, 36, 71), INT64_MIN, dtype=np.int64)
        ah = np.zeros((ntask, 36, 71), dtype=np.int64)
        cnt = np.zeros((ntask, 36, 71), dtype=np.int64)
        nbad = np.zeros(ntask, dtype=np.int64)
        for task in prange(ntask):
            ri = task // nper
            l0 = (task % nper) * lchunk
            Rm = reps[ri]
            bits = np.empty(35, dtype=np.uint8)
            fl = np.empty(20, dtype=np.int64)
            for L in range(l0, l0 + lchunk):
                h = Rm | (L << 15)
                ok = True
                for j in range(psets_new.size):
                    if popcnt64(h & psets_new[j]) < LAMC:
                        ok = False
                        break
                if not ok:
                    continue
                _unpack(h, bits)
                Z, bd, e, cd = _eval_Z(h, bits, nb, nrs, nsig, nU, nfree, npairs, wint, sigpos, freepos, pairs,
                                       type_of, flag_of, qoff, qdim, aux, Qi, fl)
                if bd:
                    nbad[task] += 1
                    continue
                cnt[task, e, cd] += 1
                if Z > mz[task, e, cd]:
                    mz[task, e, cd] = Z
                    ah[task, e, cd] = h
        return _merge(mz, ah, cnt, nbad)

    return scan_raw_exact_lam
