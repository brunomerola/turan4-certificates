"""Isomorphism classes of admissible 4-graphs on 7 vertices ("every p-set contains an edge"), p = 5, 6, 7, by
one-vertex extension of the 6-vertex representatives, with canonical forms from nauty (pynauty).

Completeness: admissibility is hereditary, so H - {6} is admissible and isomorphic to a representative R; hence
every class contains some R | (L << 15) with L a link on {0..5}.  Links in the same Aut(R)-orbit give isomorphic
extensions, so only the orbit minima are kept (exact orbit test: L <= g(L) for every g in Aut(R)).

Canonical form: the 4-graph (or its edge-complement in K_7^(4) when it has more than 17 edges, an isomorphism-
invariant choice that keeps the incidence graph small) is encoded as the vertex-coloured bipartite incidence graph
(7 vertex nodes | one node per edge); nauty's canonical labelling restricted to the vertex cell gives the relabelling
v -> inv[v] (inv[lab[i]] = i), and the canonical 35-bit mask is the relabelled edge set (complemented back).

Output: data/classes_p{p}.npy (sorted uint64 canonical masks) and results/enum_p{p}.json.
Optional --labelled-check: sum over classes of 7!/|Aut H| (nauty group orders) must equal the number of labelled
admissible graphs, computed independently as sum_R (6!/|Aut R|) * #{admissible links of R}.

Usage: python enum7.py P [--workers 2] [--labelled-check]
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import time
from itertools import combinations, permutations
from math import comb, factorial
from multiprocessing import Pool

import numpy as np

from common import DATA, E7, EIDX, RESULTS, apply_perm, pset_masks7, reps6

FULL = (1 << 35) - 1
T6 = sorted(combinations(range(6), 3), key=lambda e: e[::-1])     # link bits, colex
T6I = {t: i for i, t in enumerate(T6)}
assert all(EIDX[t + (6,)] == 15 + i for i, t in enumerate(T6))
PERMS6 = list(permutations(range(6)))
# edge -> vertex tuple, and per-relabelling colex lookup through a dict on sorted tuples
EDGE_VERTS = [list(e) for e in E7]


def canon7(mask: int) -> int:
    import pynauty
    comp = bin(mask).count("1") > 17
    m = (FULL ^ mask) if comp else mask
    es = [i for i in range(35) if (m >> i) & 1]
    k = len(es)
    if k == 0:
        c = 0
    else:
        g = pynauty.Graph(7 + k, adjacency_dict={7 + j: EDGE_VERTS[i] for j, i in enumerate(es)},
                          vertex_coloring=[set(range(7)), set(range(7, 7 + k))])
        lab = pynauty.canon_label(g)
        inv = [0] * 7
        for i in range(7):
            inv[lab[i]] = i
        c = 0
        for i in es:
            e = E7[i]
            c |= 1 << EIDX[tuple(sorted((inv[e[0]], inv[e[1]], inv[e[2]], inv[e[3]])))]
    return (FULL ^ c) if comp else c


def aut_order7(mask: int) -> int:
    import pynauty
    comp = bin(mask).count("1") > 17
    m = (FULL ^ mask) if comp else mask
    es = [i for i in range(35) if (m >> i) & 1]
    k = len(es)
    if k == 0:
        return factorial(7)
    g = pynauty.Graph(7 + k, adjacency_dict={7 + j: EDGE_VERTS[i] for j, i in enumerate(es)},
                      vertex_coloring=[set(range(7)), set(range(7, 7 + k))])
    _, grpsize1, grpsize2, _, _ = pynauty.autgrp(g)
    # the incidence graph's group acts faithfully on the vertex cell (distinct edges are distinct 4-sets)
    return int(round(grpsize1 * 10 ** grpsize2))


def automorphisms6(R: int) -> list:
    return [g for g in PERMS6 if apply_perm(R, 6, g) == R]


def link_orbit_minima(auts: list) -> np.ndarray:
    """All links L (20-bit) with L <= g(L) for every g in auts."""
    L = np.arange(1 << 20, dtype=np.int64)
    keep = np.ones(L.size, dtype=bool)
    lo, hi = L & 1023, L >> 10
    for g in auts:
        img = [T6I[tuple(sorted(g[v] for v in t))] for t in T6]
        tl = np.zeros(1024, dtype=np.int64)
        th = np.zeros(1024, dtype=np.int64)
        for j in range(10):
            tl |= ((np.arange(1024) >> j) & 1) << img[j]
            th |= ((np.arange(1024) >> j) & 1) << img[10 + j]
        keep &= (tl[lo] | th[hi]) >= L
    return L[keep]


def admissible_links(R: int, p: int) -> np.ndarray:
    L = np.arange(1 << 20, dtype=np.int64)
    H = np.int64(R) | (L << 15)
    ok = np.ones(L.size, dtype=bool)
    for pm in pset_masks7(p):
        ok &= (H & pm) != 0
    return ok


def work(args):
    R, p, labelled = args
    t0 = time.time()
    auts = automorphisms6(int(R))
    L = link_orbit_minima(auts)
    H = np.int64(R) | (L << 15)
    ok = np.ones(H.size, dtype=bool)
    for pm in pset_masks7(p):
        ok &= (H & pm) != 0
    H = H[ok]
    can = np.fromiter((canon7(int(h)) for h in H), dtype=np.int64, count=H.size)
    out = np.unique(can)
    lab = 0
    if labelled:
        lab = (factorial(6) // len(auts)) * int(admissible_links(int(R), p).sum())
    return int(R), len(auts), int(L.size), int(H.size), out, lab, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("p", type=int)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--labelled-check", action="store_true")
    ap.add_argument("--limit-reps", type=int, default=0, help="profiling: only the first k representatives")
    a = ap.parse_args()
    t0 = time.time()
    reps = reps6(a.p)
    order = sorted(range(reps.size), key=lambda i: -factorial(6) / len(automorphisms6(int(reps[i]))))
    if a.limit_reps:
        order = order[:a.limit_reps]
    print(f"# enum7 p={a.p}: {reps.size} reps on 6 vertices; workers {a.workers}", flush=True)
    acc, pending = np.zeros(0, dtype=np.int64), []
    stats = {"candidates_orbit_min": 0, "admissible_candidates": 0, "labelled_admissible": 0}
    per_rep = []
    with Pool(a.workers) as pool:
        for k, (R, na, nl, nh, can, lab, dt) in enumerate(
                pool.imap_unordered(work, [(int(reps[i]), a.p, a.labelled_check) for i in order])):
            stats["candidates_orbit_min"] += nl
            stats["admissible_candidates"] += nh
            stats["labelled_admissible"] += lab
            per_rep.append({"R": R, "aut": na, "orbit_min_links": nl, "admissible": nh, "classes_from_R": int(can.size),
                            "s": round(dt, 2)})
            pending.append(can)
            if sum(x.size for x in pending) > 20_000_000:
                acc = np.unique(np.concatenate([acc] + pending))
                pending = []
            if k % 10 == 0 or k == len(order) - 1:
                print(f"  {k + 1}/{len(order)} reps, {stats['admissible_candidates']:,} candidates, "
                      f"{time.time() - t0:.0f}s, RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.0f} MB",
                      flush=True)
    acc = np.unique(np.concatenate([acc] + pending))
    t_enum = time.time() - t0
    out = {"p": a.p, "n_reps6": int(reps.size), "n_reps_used": len(order), "classes": int(acc.size),
           "time_s": t_enum, **stats}
    if a.labelled_check:
        t1 = time.time()
        with Pool(a.workers) as pool:
            orders = pool.map(aut_order7, [int(x) for x in acc], chunksize=20000)
        s = sum(factorial(7) // o for o in orders)
        out.update({"labelled_from_classes": s, "labelled_check_ok": s == stats["labelled_admissible"],
                    "labelled_check_s": time.time() - t1})
    out["peak_rss_mb_main"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    out["peak_rss_mb_children"] = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024.0
    if not a.limit_reps:
        np.save(os.path.join(DATA, f"classes_p{a.p}.npy"), acc.astype(np.uint64))
        with open(os.path.join(RESULTS, f"enum_p{a.p}.json"), "w") as f:
            json.dump({**out, "per_rep": per_rep}, f, indent=1)
    print(json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
