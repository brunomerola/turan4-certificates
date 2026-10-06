"""Sanity check of an H44 certificate on Sidorenko's K_5^{4-}-free construction (density 3/7), complement form.

Construction (Sidorenko, JCTB 2024, arXiv 2205.02006 Section 1): vertices = L-bit strings; a 4-set is an edge of G
iff, at the first bit position where its four strings do not all agree, they split 2 + 2.  (Recursive: two copies
plus all 4-sets with two vertices in each copy; density x = 6/16 + (2/16) x -> 3/7.)  Every 5-set spans <= 3 edges
(checked here on every sampled 7-set).  Its complement Gc is admissible for t_2(5,4) (every 5-set spans >= 2
edges) with density -> 4/7.
The flag-algebra inequality says E_H[d(H) - c(H)] >= b for random 7-sets H of a large admissible graph, and
E_H[c(H)] >= -o(1) (c = the SOS term).  We sample 7-sets of Gc (L = 40, distinct strings), evaluate d and c with
the n7 float kernel and the certificate's Gram matrices.  Float, sanity only.
Usage: python plugin_check.py RUN_PREFIX [--samples 200000]
"""
import argparse
from itertools import combinations

import numpy as np

import gen

gen.configure(5, 2)
import common  # noqa: E402
from kernels import qsum_flat, scan_float, tables  # noqa: E402

E7 = common.E7


def edge_rule(strings):
    """strings: (4,) uint64 -> is the 4-set an edge of the construction."""
    a = strings
    diff = np.bitwise_or.reduce(a ^ a[0])
    top = int(diff).bit_length() - 1          # highest differing bit (first position from the top)
    ones = int(((a >> np.uint64(top)) & np.uint64(1)).sum())
    return ones == 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix")
    ap.add_argument("--samples", type=int, default=200000)
    ap.add_argument("--bits", type=int, default=40)
    a = ap.parse_args()
    src = np.load(a.prefix + ".best.npz")
    blocks = [tuple(int(x) for x in b) for b in src["blocks"]]
    model = common.build_model(5, blocks)
    Fs = [src[f"F{k}"] for k in range(len(model.keys))]
    Qs = [F.T @ F for F in Fs]
    qs = qsum_flat(model, Qs)
    rng = np.random.default_rng(0)
    H = np.empty(a.samples, dtype=np.int64)
    for i in range(a.samples):
        while True:
            v = rng.integers(0, 1 << a.bits, size=7, dtype=np.uint64)
            if np.unique(v).size == 7:
                break
        m = 0
        for j, e in enumerate(E7):
            if not edge_rule(v[list(e)]):      # complement
                m |= 1 << j
        H[i] = m
    ok = gen.admissible7_lam(H, 5, 2)
    sl, bad = scan_float(H, *tables(model), qs, 0.0, 0.0, 0.0, 8)
    d = np.array([bin(int(h)).count("1") for h in H]) / 35.0
    print(f"samples {a.samples}: all admissible {bool(ok.all())}, kernel bad {int(bad.sum())}")
    print(f"E[d] = {d.mean():.6f} (4/7 = {4/7:.6f}), E[c] = {(d - sl).mean():.6f}, E[d - c] = {sl.mean():.6f} "
          f"(+- {sl.std() / np.sqrt(sl.size):.6f}); certificate float b* = {float(src['b_star']):.9f}")


if __name__ == "__main__":
    main()
