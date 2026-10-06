"""Isomorphism classes of 7-vertex 4-graphs in which every p-set spans >= LAM edges (gen.py), by one-vertex
extension of the 6-vertex representatives with nauty canonical forms -- the method of ../n7/enum7.py (its canon7,
automorphisms6, link_orbit_minima and aut_order7 are reused; only the admissibility filter differs).
Completeness: the predicate is hereditary, so H - {6} is admissible, isomorphic to a representative R, and links in
one Aut(R)-orbit give isomorphic extensions.
--labelled-check: sum over classes of 7!/|Aut H| must equal sum_R (6!/|Aut R|) * #{admissible links of R}.
Output: data/l{lam}/classes_p{p}.npy, results/enum_p{p}_l{lam}.json.
Usage: python enum_gen.py P LAM [--workers 2] [--labelled-check]
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import sys
import time
from math import factorial
from multiprocessing import Pool

import numpy as np

import gen

P_, LAM_ = int(sys.argv[1]), int(sys.argv[2])
gen.configure(P_, LAM_)

import common  # noqa: E402
from enum7 import aut_order7, automorphisms6, canon7, link_orbit_minima  # noqa: E402


def admissible_ext(R, L, p, lam):
    H = np.int64(R) | (L << 15)
    ok = np.ones(H.size, dtype=bool)
    for pm in common.pset_masks7(p):          # p-sets through vertex 6 (the others lie in R, admissible)
        ok &= np.bitwise_count(H & pm) >= lam
    return H, ok


def work(args):
    R, p, lam, labelled = args
    t0 = time.time()
    auts = automorphisms6(int(R))
    L = link_orbit_minima(auts)
    H, ok = admissible_ext(R, L, p, lam)
    H = H[ok]
    can = np.fromiter((canon7(int(h)) for h in H), dtype=np.int64, count=H.size)
    out = np.unique(can)
    lab = 0
    if labelled:
        _, okall = admissible_ext(R, np.arange(1 << 20, dtype=np.int64), p, lam)
        lab = (factorial(6) // len(auts)) * int(okall.sum())
    return int(R), len(auts), int(L.size), int(H.size), out, lab, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("p", type=int)
    ap.add_argument("lam", type=int)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--labelled-check", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    reps = common.reps6(a.p)
    # sanity: the representatives satisfy the threshold predicate
    assert gen.admissible7_lam(np.array([0]), a.p, a.lam).size == 1
    print(f"# enum_gen p={a.p} lam={a.lam}: {reps.size} reps on 6 vertices; workers {a.workers}", flush=True)
    acc, pending = np.zeros(0, dtype=np.int64), []
    stats = {"candidates_orbit_min": 0, "admissible_candidates": 0, "labelled_admissible": 0}
    per_rep = []
    with Pool(a.workers) as pool:
        for k, (R, na, nl, nh, can, lab, dt) in enumerate(
                pool.imap_unordered(work, [(int(r), a.p, a.lam, a.labelled_check) for r in reps])):
            stats["candidates_orbit_min"] += nl
            stats["admissible_candidates"] += nh
            stats["labelled_admissible"] += lab
            per_rep.append({"R": R, "aut": na, "orbit_min_links": nl, "admissible": nh,
                            "classes_from_R": int(can.size), "s": round(dt, 2)})
            pending.append(can)
            if sum(x.size for x in pending) > 20_000_000:
                acc = np.unique(np.concatenate([acc] + pending))
                pending = []
            if k % 10 == 0 or k == reps.size - 1:
                print(f"  {k + 1}/{reps.size} reps, {stats['admissible_candidates']:,} candidates, "
                      f"{time.time() - t0:.0f}s", flush=True)
    acc = np.unique(np.concatenate([acc] + pending))
    assert gen.admissible7_lam(acc, a.p, a.lam).all()
    out = {"p": a.p, "lam": a.lam, "n_reps6": int(reps.size), "classes": int(acc.size), "time_s": time.time() - t0,
           **stats}
    if a.labelled_check:
        t1 = time.time()
        with Pool(a.workers) as pool:
            orders = pool.map(aut_order7, [int(x) for x in acc], chunksize=20000)
        s = sum(factorial(7) // o for o in orders)
        out.update({"labelled_from_classes": s, "labelled_check_ok": s == stats["labelled_admissible"],
                    "labelled_check_s": time.time() - t1})
    out["peak_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    os.makedirs(gen.data_dir(a.lam), exist_ok=True)
    np.save(gen.classes_path(a.p, a.lam), acc.astype(np.uint64))
    with open(os.path.join(gen.RESULTS, f"enum_p{a.p}_l{a.lam}.json"), "w") as f:
        json.dump({**out, "per_rep": per_rep}, f, indent=1)
    print(json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
