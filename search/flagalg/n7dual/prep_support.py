"""Build the dual-SDP data for a candidate support set of admissible 7-vertex 4-graphs (p = 5, all six blocks).

For each key k (block (s, m), type sigma) and each graph H the moment matrix is
    M_k(H) = sum_{h in Aut(sigma)} Pi_h^T Mt_k(H) Pi_h,     Mt_k(H) = sum_terms c (E_ab + E_ba),  c = w / (2 * 5040),
with the canonical-root-labelling terms (a, b, w) of the n7 kernels (common.py: <Q, M(H)> = sum w Qsum[a,b] / 5040).
This module computes the terms (numba kernel `terms` of ../n7/kernels.py, as in lpcg7.add_graphs), performs the
facial reduction of facial.py and writes, for the alive graphs and the live flags only,
    masks, e (edge counts), and per key: live flag indices, Aut permutations on live flags, terms (i, a, b, w).
Float-free: all stored weights are integers w (the matrix entry is w / 10080).

Usage: python prep_support.py --masks FILE.npy|FILE.npz[:key] --out results/support_X.npz
"""
import argparse
import json
import time

import numpy as np

import paths
from common import build_model
from facial import live_support
from kernels import tables, terms
from lpcg7 import ALL_BLOCKS


def graph_terms(model, H):
    """Per key: (rows, a, b, w) with a <= b (canonical labelling), duplicates merged (w summed, int64)."""
    tab = tables(model)
    nk = len(model.keys)
    keyid = np.full((len(model.blocks), 6), -1, dtype=np.int64)
    for k, (bi, t) in enumerate(model.keys):
        keyid[bi, t] = k
    out = [[] for _ in range(nk)]
    ee = np.zeros(H.size, dtype=np.int64)
    chunk = 20000
    for c0 in range(0, H.size, chunk):
        Hc = H[c0:c0 + chunk]
        rows, keys, aa, bb, ww, e, _ = terms(Hc, *tab[:12], keyid, model.cdd_pos)
        ee[c0:c0 + Hc.size] = e
        code = ((rows * nk + keys) * 1024 + aa) * 1024 + bb
        u, inv = np.unique(code, return_inverse=True)
        wsum = np.bincount(inv, weights=ww).astype(np.int64)
        bb_ = u % 1024
        aa_ = (u // 1024) % 1024
        kk_ = (u // (1024 * 1024)) % nk
        rr_ = u // (1024 * 1024 * nk) + c0
        for k in range(nk):
            s = kk_ == k
            out[k].append((rr_[s], aa_[s], bb_[s], wsum[s]))
    T = []
    for k in range(nk):
        T.append([np.concatenate([x[j] for x in out[k]]) for j in range(4)])
    return T, ee


def build(masks, model=None, verbose=True):
    model = build_model(5, ALL_BLOCKS) if model is None else model
    H = np.unique(np.asarray(masks, dtype=np.int64))
    T, ee = graph_terms(model, H)
    P = [model.blocks[bi].aut[t] for bi, t in model.keys]
    dims = [int(model.qdim[bi, t]) for bi, t in model.keys]
    alive, live = live_support(T, P, dims, H.size, verbose=verbose)
    idx = -np.ones(H.size, dtype=np.int64)
    idx[alive] = np.arange(int(alive.sum()))
    data = {"masks": H[alive], "e": ee[alive], "nkeys": len(dims)}
    for k in range(len(dims)):
        lv = np.nonzero(live[k])[0]
        lmap = -np.ones(dims[k], dtype=np.int64)
        lmap[lv] = np.arange(lv.size)
        r, a, b, w = T[k]
        sel = alive[r]
        assert (lmap[a[sel]] >= 0).all() and (lmap[b[sel]] >= 0).all()
        Pl = lmap[P[k][:, lv]]
        assert (Pl >= 0).all()
        data[f"live{k}"] = lv
        data[f"P{k}"] = Pl
        data[f"ti{k}"] = idx[r[sel]]
        data[f"ta{k}"] = lmap[a[sel]]
        data[f"tb{k}"] = lmap[b[sel]]
        data[f"tw{k}"] = w[sel]
        data[f"dim{k}"] = dims[k]
    return data, alive, live


def load_masks(spec):
    if ":" in spec:
        f, key = spec.split(":")
        return np.load(f)[key]
    x = np.load(spec)
    return x["W"] if hasattr(x, "files") else x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--masks", action="append", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    t0 = time.time()
    masks = np.concatenate([load_masks(s).astype(np.int64) for s in a.masks])
    data, alive, live = build(masks)
    np.savez(a.out, **data)
    info = {"candidates": int(np.unique(masks).size), "alive": int(alive.sum()),
            "live": [int(l.sum()) for l in live], "time_s": time.time() - t0}
    print(json.dumps(info), flush=True)
    with open(a.out.replace(".npz", ".json"), "w") as f:
        json.dump(info, f, indent=1)


if __name__ == "__main__":
    main()
