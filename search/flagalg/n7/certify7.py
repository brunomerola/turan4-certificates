"""Exact certificate for t(p, 4) >= b* from a float LP-CUT-CG solution (lpcg7.py --out X writes X.best.npz).

1. Rounding (Jeong et al. 2026, Appendix B.3 style; independent implementation): for every Gram Q_sigma = F^T F,
   eigendecompose, drop eigenvalues <= drop * max, and round the scaled factors to integers, a_i = round(M sqrt(w_i) u_i);
   Qint = sum_i a_i a_i^T (integer, PSD by construction), Q = Qint / M^2.  tau_j = floor(tau_j M^2) / M^2 >= 0.
2. Exact evaluation, independent of floats: for every admissible H the integer
        Z(H) = sum_terms wint * Qsum[a, b]   (int64, with an a-priori overflow bound)
   so that <Q, M(H)> = Z(H) / (5040 M^2).  The slack is affine in (e, cdd) and decreasing in Z:
        a_H  = d - Z/(5040 M^2) + tau1 c_dd + tau2 (d - c_dd),    beta_H = tau1 d + tau2 (1 - d),
        slack(H; b0) = a_H - b0 beta_H,
   so the kernel returns max Z per group (e, cdd), and Python computes in exact rationals
        b* = min_H a_H / (1 + beta_H).
   Proof of t >= b*: with b0 = b*, every H has slack(H; b*) >= b*(1 + beta_H) - b* beta_H = b*, so the flag-algebra
   inequality d - b* >= sum <Q, M> + tau1 (b* d - [[d.d]]) + tau2 (b*(1-d) - [[d(1-d)]]) holds on every large
   admissible graph up to o(1); on degree-regular extremal graphs with t < b* the right side is >= 0 (Q PSD, tau >= 0,
   (t - b*) [[F]] <= 0), a contradiction.  Hence t >= b*.
3. Two exact scans: over the nauty class list (data/classes_p{p}.npy) and over ALL raw one-vertex extensions
   R | (L << 15) of the 6-vertex representatives (no canonical forms at N = 7; completeness by heredity).  Both must
   give the same group maxima.

Usage: python certify7.py RUN_PREFIX [--M-bits auto] [--no-raw]
Writes RUN_PREFIX.cert.npz (integer factor rows) and RUN_PREFIX.cert.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from fractions import Fraction

import numpy as np

from common import DATA, build_model, reps6
from kernels import qsum_flat, scan_exact, scan_raw_exact, tables

I63 = (1 << 63) - 1


def round_factors(F, M, drop):
    """F (k x n) float factor rows -> integer rows A (r x n), Q ~ A^T A / M^2."""
    n = F.shape[1]
    if F.shape[0] == 0:
        return np.zeros((0, n), dtype=np.int64)
    Q = F.T @ F
    w, U = np.linalg.eigh(Q)
    wmax = max(float(w.max()), 0.0)
    keep = w > drop * max(wmax, 1e-300)
    A = np.rint(M * (np.sqrt(w[keep])[:, None] * U[:, keep].T)).astype(np.int64)
    A = A[np.abs(A).sum(axis=1) > 0]
    return A


def exact_gram(A):
    """A^T A exactly (int64 after an overflow check on the diagonal, which bounds every entry of a PSD Gram)."""
    diag = (A.astype(object) ** 2).sum(axis=0) if A.size else np.zeros(A.shape[1], dtype=object)
    if A.size and max(int(x) for x in diag) > I63 // 2:
        raise OverflowError("Gram diagonal exceeds int64")
    return A.T @ A if A.size else np.zeros((A.shape[1], A.shape[1]), dtype=np.int64)


def bound_from_groups(MZ, CN, M2, t1, t2):
    """Exact b* = min over groups (e, cdd) of a/(1 + beta), see the module docstring."""
    best, arg = None, None
    D = 5040 * M2
    for e in range(36):
        for c in range(71):
            if CN[e, c] == 0:
                continue
            Zm = int(MZ[e, c])
            A = 144 * e * M2 - Zm + t1 * 72 * c + t2 * (144 * e - 72 * c)
            Bt = t1 * 144 * e + t2 * (5040 - 144 * e)
            val = Fraction(A, D + Bt)
            if best is None or val < best:
                best, arg = val, (e, c)
    return best, arg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix")
    ap.add_argument("--M-bits", type=int, default=0, help="0 = largest power of 2 passing the overflow bound")
    ap.add_argument("--drop", type=float, default=1e-11)
    ap.add_argument("--no-raw", action="store_true")
    ap.add_argument("--nchunks", type=int, default=64)
    ap.add_argument("--lchunk", type=int, default=1 << 18)
    a = ap.parse_args()
    t0 = time.time()
    src = np.load(a.prefix + ".best.npz")
    p = int(src["p"])
    blocks = [tuple(int(x) for x in b) for b in src["blocks"]]
    model = build_model(p, blocks)
    assert [tuple(k) for k in src["keys"]] == [tuple(k) for k in model.keys]
    Fs = [src[f"F{k}"] for k in range(len(model.keys))]
    # scale M: |Z| <= 5040 * nblocks * max_key max|Qint| and max|Qint| <= max diag ~ M^2 * max diag(Q)
    qmax = max([float((F ** 2).sum(axis=0).max()) for F in Fs if F.size] + [1e-300])
    if a.M_bits:
        mb = a.M_bits
    else:
        mb = 1
        while (2 ** (2 * (mb + 1))) * qmax * 1.05 * 5040 * len(blocks) < I63 / 2 and mb < 40:
            mb += 1
    M = 1 << mb
    M2 = M * M
    As = [round_factors(F, M, a.drop) for F in Fs]
    Qi = [exact_gram(A) for A in As]
    zbound = 5040 * sum(max(int(np.abs(q).max()) if q.size else 0 for q, (bi, t) in zip(Qi, model.keys) if bi == b)
                        for b in range(len(blocks)))
    if zbound > I63:
        raise OverflowError(f"a-priori |Z| bound {zbound} exceeds int64; lower --M-bits")
    qi = qsum_flat(model, Qi, dtype=np.int64)
    t1 = int(np.floor(float(src["tau1"]) * M2))
    t2 = int(np.floor(float(src["tau2"]) * M2))
    t1, t2 = max(t1, 0), max(t2, 0)
    tab = tables(model)
    cls = np.load(os.path.join(DATA, f"classes_p{p}.npy")).astype(np.int64)
    t = time.time()
    MZ, AH, CN, nbad = scan_exact(cls, *tab, qi, a.nchunks)
    t_cls = time.time() - t
    assert nbad == 0
    bstar, arg = bound_from_groups(MZ, CN, M2, t1, t2)
    out = {"p": p, "blocks": blocks, "M": M, "M_bits": mb, "tau1_num": t1, "tau2_num": t2, "tau_den": M2,
           "float_b_star": float(src["b_star"]), "float_b_lp": float(src["b_lp"]),
           "bound": str(bstar), "bound_float": float(bstar), "argmin_group": arg,
           "argmin_graph": int(AH[arg]), "z_bound_apriori": str(zbound),
           "factor_rows": [int(A.shape[0]) for A in As], "classes_scanned": int(cls.size),
           "class_scan_s": t_cls}
    print(f"# exact (class scan, {cls.size:,} graphs, {t_cls:.1f}s): b* = {float(bstar):.12f} "
          f"[float LP b* {float(src['b_star']):.12f}]", flush=True)
    if not a.no_raw:
        reps = reps6(p).astype(np.int64)
        t = time.time()
        MZr, AHr, CNr, nbr = scan_raw_exact(reps, model.psets_new, *tab, qi, a.lchunk)
        t_raw = time.time() - t
        same = bool(np.array_equal(MZr[CNr > 0], MZ[CN > 0]) and np.array_equal(CNr > 0, CN > 0))
        braw, _ = bound_from_groups(MZr, CNr, M2, t1, t2)
        out.update({"raw_extensions_scanned": int(CNr.sum()), "raw_scan_s": t_raw, "raw_bad": int(nbr),
                    "raw_bound": str(braw), "raw_groups_identical": same})
        print(f"# exact (raw scan, {int(CNr.sum()):,} extensions, {t_raw:.1f}s): b* = {float(braw):.12f}; "
              f"group maxima identical to class scan: {same}", flush=True)
    np.savez_compressed(a.prefix + ".cert.npz", **{f"A{k}": A for k, A in enumerate(As)})
    with open(a.prefix + ".cert.npz", "rb") as f:
        out["cert_npz_sha256"] = hashlib.sha256(f.read()).hexdigest()
    out["groups"] = [[e, c, str(int(MZ[e, c])), int(AH[e, c]), int(CN[e, c])] for e in range(36) for c in range(71)
                     if CN[e, c]]
    out["keys"] = [{"s": model.blocks[bi].s, "m": model.blocks[bi].m, "type_index": t,
                    "sigma": model.blocks[bi].types[t], "n_flags": int(model.qdim[bi, t]),
                    "flags": [int(f) for f in model.blocks[bi].flags[t]]} for bi, t in model.keys]
    out["time_s"] = time.time() - t0
    with open(a.prefix + ".cert.json", "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k not in ("groups",)}), flush=True)


if __name__ == "__main__":
    main()
