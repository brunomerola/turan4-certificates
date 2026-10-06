"""Throughput measurements for the N = 7 kernels (pricing scans, exact scans, LP row features), for the cloud
pre-launch estimate.  Usage: python profile7.py P [--n 1000000]  (run with taskset -c 0,1)."""
import argparse
import json
import os
import resource
import time

import numba
import numpy as np

from common import DATA, RESULTS, build_model, reps6
from kernels import coef_columns, qsum_flat, scan_exact, scan_float, scan_raw_float_min, tables, terms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("p", type=int)
    ap.add_argument("--n", type=int, default=1_000_000)
    a = ap.parse_args()
    model = build_model(a.p)
    tab = tables(model)
    cls = np.load(os.path.join(DATA, f"classes_p{a.p}.npy")).astype(np.int64)
    rng = np.random.default_rng(0)
    H = np.sort(rng.choice(cls, size=min(a.n, cls.size), replace=False))
    Qs = []
    for bi, t in model.keys:
        n = int(model.qdim[bi, t])
        B = rng.normal(size=(min(n, 64), n))
        Qs.append(B.T @ B / n)
    t = time.time()
    qs = qsum_flat(model, Qs)
    t_qsum = time.time() - t
    qi = qsum_flat(model, [np.rint(Q * 1e6).astype(np.int64) for Q in Qs], dtype=np.int64)
    out = {"p": a.p, "classes": int(cls.size), "sample": int(H.size), "qsum_s": t_qsum}
    # warm-up (compile)
    scan_float(H[:1000], *tab, qs, 0.1, 0.1, 0.1, 8)
    scan_exact(H[:1000], *tab, qi, 8)
    reps = reps6(a.p).astype(np.int64)
    scan_raw_float_min(reps[:1], model.psets_new, *tab, qs, 0.1, 0.1, 0.1, 1 << 18)
    for nt in (1, 2):
        numba.set_num_threads(nt)
        t = time.time()
        scan_float(H, *tab, qs, 0.1, 0.1, 0.1, 64)
        dt = time.time() - t
        out[f"float_scan_us_per_graph_{nt}thr"] = dt / H.size * 1e6
        t = time.time()
        scan_exact(H, *tab, qi, 64)
        dt = time.time() - t
        out[f"exact_scan_us_per_graph_{nt}thr"] = dt / H.size * 1e6
        t = time.time()
        mn, am, na = scan_raw_float_min(reps[:2], model.psets_new, *tab, qs, 0.1, 0.1, 0.1, 1 << 18)
        dt = time.time() - t
        out[f"raw_float_scan_us_per_extension_{nt}thr"] = dt / (2 << 20) * 1e6
        out["raw_admissible_frac_first2reps"] = float(na.sum()) / (2 << 20)
    numba.set_num_threads(2)
    # LP row features: terms + coefficients of 80 columns of a 1024-dim key on 1000 rows
    keyid = np.full((len(model.blocks), 6), -1, dtype=np.int64)
    for k, (bi, t_) in enumerate(model.keys):
        keyid[bi, t_] = k
    Hs = H[:1000]
    terms(Hs[:10], *tab[:12], keyid, model.cdd_pos)
    t = time.time()
    rows, keys, aa, bb, ww, ee, cc = terms(Hs, *tab[:12], keyid, model.cdd_pos)
    out["terms_us_per_row"] = (time.time() - t) / Hs.size * 1e6
    out["terms_per_row"] = rows.size / Hs.size
    big = [k for k, (bi, t_) in enumerate(model.keys) if model.qdim[bi, t_] >= 800]
    for k in (big[0], big[-1]):
        bi, t_ = model.keys[k]
        P = model.blocks[bi].aut[t_]
        sel = keys == k
        V = rng.normal(size=(80, int(model.qdim[bi, t_])))
        coef_columns(rows[sel][:10], aa[sel][:10], bb[sel][:10], ww[sel][:10], P, V[:2], Hs.size)
        t = time.time()
        coef_columns(rows[sel], aa[sel], bb[sel], ww[sel], P, V, Hs.size)
        out[f"coef_us_per_col_per_row_key{k}_aut{P.shape[0]}"] = (time.time() - t) / 80 / Hs.size * 1e6
        out[f"terms_per_row_key{k}"] = int(sel.sum()) / Hs.size
    out["peak_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    print(json.dumps(out, indent=1))
    with open(os.path.join(RESULTS, f"profile_kernels_p{a.p}.json"), "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
