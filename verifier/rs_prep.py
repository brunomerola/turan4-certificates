"""Independent checker (own code): rebuild the certificate's Gram matrices EXACTLY, check them, check the flag lists, and write
the multi-limb blob read by rs_eval.c.  No producer module is imported; only producer DATA files are read.

Inputs
  --sharp NPZ        the sharp certificate (results/sharp_klp5.cert.npz): kb, q, D (decimal strings), W_k (decimal
                     strings, integer factor rows), B_k (int64), corr (JSON: key, X_times_q).  Claimed
                         Q_k = W_k^T W_k / 2^(2 kb) + B_k^T (X_k) B_k,   X_k = X_times_q / q,   D = 2^(2 kb) q,
                     so the integer Gram  G_k = D Q_k = q W_k^T W_k + 2^(2 kb) B_k^T X_times_q B_k.
  --cloud NPZ        alternatively the non-sharp cloud certificate (A_k int64, M): G_k = A_k^T A_k, D = M^2.
  --keys JSON        key list with flag masks: a producer cert.json whose 'keys' carry (s, m, sigma, flags) in COLEX
                     (certificates/sigma_K5_4/sharp_klp5.cert.json, taken from certificates/K5_4/lpcg_p5_conv_full.cert.json);
  --keys-check JSON  a second key list (certificates/K5_4/lpcg_p5_conv_full.cert.json) -- (s, m, sigma, n_flags) must agree.
  --target FRAC      the claimed bound b; the evaluator reports S(H) = 5040 D (slack(H) - b) exactly.
Checks done here (all exact):
  * every flag decodes to a flag whose root part == sigma (labelled); flags pairwise non-isomorphic (canonical form =
    min over non-root permutations); all admissible (every 5-set has an edge); list complete (= all admissible
    sigma-flag classes); types: admissible, distinct per block;
  * G_k symmetric; X_k symmetric and PSD (own fraction LDL^T with the zero-pivot rule);
  * limb decomposition G = sum_j L_j 2^(45 j), |L_j| <= 2^44 (balanced), recombination checked; int64 accumulation
    bound checked.
Usage: python rs_prep.py (--sharp NPZ | --cloud NPZ) --keys JSON --keys-check JSON --target FRAC --out BLOB
"""
import argparse
import json
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial

import numpy as np

from rs_lib import admissible, canon_flag, colex_decode, dec, enc, flag_edge_order

LB = 45                      # limb bits
HALF = 1 << (LB - 1)


def ldl_psd(M):
    """exact PSD test of a symmetric rational matrix (list of lists of Fractions)"""
    A = [row[:] for row in M]
    n = len(A)
    for i in range(n):
        for j in range(n):
            if A[i][j] != A[j][i]:
                return False, "not symmetric"
    for k in range(n):
        d = A[k][k]
        if d < 0:
            return False, f"negative pivot at {k}"
        if d == 0:
            if any(A[k][j] != 0 for j in range(k + 1, n)):
                return False, f"zero pivot with nonzero row at {k}"
            continue
        for i in range(k + 1, n):
            if A[i][k] == 0:
                continue
            f = A[i][k] / d
            for j in range(k + 1, n):
                A[i][j] -= f * A[k][j]
            A[i][k] = Fraction(0)
    return True, "ok"


def int_gram(W):
    """exact W^T W for an integer matrix (Python ints, |W| < 2^40) via an int64 two-part split"""
    W = np.array(W, dtype=object)
    r, n = W.shape
    if r == 0:
        return np.zeros((n, n), dtype=object)
    mx = max(abs(int(x)) for x in W.ravel())
    assert mx < (1 << 40)
    lo = np.array([[((int(x) + (1 << 19)) % (1 << 20)) - (1 << 19) for x in row] for row in W], dtype=np.int64)
    hi_o = (W - lo.astype(object))
    assert all(int(x) % (1 << 20) == 0 for x in hi_o.ravel())
    hi = np.array([[int(x) >> 20 for x in row] for row in hi_o], dtype=np.int64)
    assert np.abs(hi).max() < (1 << 21) and np.abs(lo).max() <= (1 << 19)
    # each product sum: r * 2^42 < 2^51 -> exact in int64
    assert r * (1 << 42) < (1 << 62)
    hh = (hi.T @ hi).astype(object)
    hl = (hi.T @ lo).astype(object)
    ll = (lo.T @ lo).astype(object)
    G = hh * (1 << 40) + (hl + hl.T) * (1 << 20) + ll
    # spot check against a direct Python-int computation on random entries
    rng = np.random.default_rng(1)
    for _ in range(50):
        a, c = rng.integers(0, n, 2)
        assert G[a, c] == sum(int(W[i, a]) * int(W[i, c]) for i in range(r))
    return G


def limbs_of(x, NL):
    out = []
    for j in range(NL - 1):
        d = ((x + HALF) % (1 << LB)) - HALF
        out.append(d)
        x = (x - d) >> LB
    out.append(x)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sharp")
    ap.add_argument("--cloud")
    ap.add_argument("--keys", required=True)
    ap.add_argument("--keys-check", required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sym-out", help="also write the Aut(sigma)-symmetrised blob for rs_eval.c -DSYM")
    a = ap.parse_args()
    target = Fraction(a.target)
    kj = json.load(open(a.keys))
    keys = kj["keys"]
    kc = json.load(open(a.keys_check))["keys"]
    assert len(keys) == len(kc)
    for x, y in zip(keys, kc):
        assert (x["s"], x["m"], x["sigma"], len(x["flags"])) == (y["s"], y["m"], y["sigma"], y["n_flags"]), (x, y)
    rep = {"keys_file": a.keys, "target": str(target), "keys": []}

    # ---------------------------------------------------------------- Gram matrices
    Gs = []
    if a.sharp:
        z = np.load(a.sharp)
        kb, q, D = int(z["kb"]), int(str(z["q"])), int(str(z["D"]))
        D0 = 1 << (2 * kb)
        assert D == D0 * q
        corr = {c["key"]: c for c in json.loads(str(z["corr"]))}
        rep.update({"kb": kb, "q": str(q), "D": str(D), "D_factor": "57 * 2^152" if D == 57 * (1 << 152) else "?"})
        for k, key in enumerate(keys):
            n = len(key["flags"])
            Wraw = z[f"W{k}"]
            W = [[int(x) for x in row] for row in Wraw] if Wraw.size else []
            if W:
                assert len(W[0]) == n
            G = int_gram(np.array(W, dtype=object).reshape(-1, n)) * q
            info = {"k": k, "W_rows": len(W)}
            if k in corr:
                B = z[f"B{k}"].astype(object)
                Xq = [[int(x) for x in row] for row in corr[k]["X_times_q"]]
                assert B.shape == (len(Xq), n)
                X = [[Fraction(v, q) for v in row] for row in Xq]
                psd, why = ldl_psd(X)
                Xo = np.array(Xq, dtype=object)
                G = G + (B.T @ (Xo @ B)) * D0
                info.update({"B_rows": int(B.shape[0]), "B_max": int(max(abs(int(v)) for v in B.ravel())),
                             "X_dim": len(Xq), "X_psd_exact": psd, "X_psd_note": why,
                             "X_min_diag": str(min(X[i][i] for i in range(len(X)))),
                             "X_max_offdiag": str(max([abs(X[i][j]) for i in range(len(X)) for j in range(len(X))
                                                       if i != j] + [Fraction(0)]))})
                assert psd, (k, why)
            Gs.append(G)
            rep["keys"].append(info)
        extra = sorted(set(corr) - set(range(len(keys))))
        assert not extra
        # every stored array is used
        used = {"kb", "q", "D", "eps", "corr"} | {f"W{k}" for k in range(len(keys))} | {f"B{k}" for k in corr}
        rep["npz_unused_fields"] = sorted(set(z.files) - used)
    else:
        z = np.load(a.cloud)
        cj = json.load(open(a.cloud.replace(".npz", ".json")))
        M = int(cj["M"])
        D = M * M
        rep.update({"M": M, "D": str(D)})
        for k, key in enumerate(keys):
            A = z[f"A{k}"].astype(object)
            n = len(key["flags"])
            G = (A.T @ A) if A.shape[0] else np.zeros((n, n), dtype=object)
            Gs.append(G)
            rep["keys"].append({"k": k, "A_rows": int(A.shape[0])})
    # symmetry
    for k, G in enumerate(Gs):
        assert G.shape[0] == G.shape[1] == len(keys[k]["flags"])
        assert (G == G.T).all(), f"G{k} not symmetric"
    gmax = max([abs(int(x)) for G in Gs for x in G.ravel()] + [1])
    NL = 1
    while (gmax >> (LB * (NL - 1))) >= HALF:
        NL += 1
    rep["G_max_bits"] = gmax.bit_length()
    rep["limbs"] = NL

    # ---------------------------------------------------------------- flags and types
    blocks = []
    for k, key in enumerate(keys):
        sm = (key["s"], key["m"])
        if not blocks or blocks[-1][0] != sm:
            assert all(b[0] != sm for b in blocks)
            blocks.append((sm, []))
        blocks[-1][1].append(k)
    blob_blocks = []
    canon_all, flags_all = {}, {}
    for (s, m), ks in blocks:
        assert 2 * m - s <= 7
        order = flag_edge_order(s, m)
        nE, nsig = len(order), comb(s, 4)
        sig_of = {}
        canon_lists = {}
        for k in ks:
            key = keys[k]
            sig_es = colex_decode(key["sigma"], s)
            sx = enc(sig_es, order)
            assert sx < (1 << nsig) or (nsig == 0 and sx == 0)
            assert sx not in sig_of
            assert s < 5 or admissible(sig_es, s)
            sig_of[sx] = k
            canon = []
            n_adm = 0
            for f in key["flags"]:
                es = colex_decode(f, m)
                assert {e for e in es if e <= frozenset(range(s))} == sig_es, "root part != sigma"
                n_adm += admissible(es, m)
                canon.append(canon_flag(es, s, m, order))
            allc = set()
            for fb in range(1 << (nE - nsig)):
                es = sig_es | dec(fb << nsig, order)
                if admissible(es, m):
                    allc.add(canon_flag(es, s, m, order))
            info = rep["keys"][k]
            info.update({"s": s, "m": m, "sigma_colex": key["sigma"], "sigma_edges": len(sig_es),
                         "n_flags": len(canon), "distinct_classes": len(set(canon)), "all_admissible": n_adm == len(canon),
                         "admissible_classes": len(allc), "missing": len(allc - set(canon)),
                         "extra": len(set(canon) - allc), "G_zero": not any(Gs[k].ravel())})
            assert len(set(canon)) == len(canon), "duplicate flag classes"
            assert n_adm == len(canon) and not (set(canon) - allc)
            canon_lists[k] = {c: i for i, c in enumerate(canon)}
            canon_all[k] = canon_lists[k]
            flags_all[k] = [colex_decode(f, m) for f in key["flags"]]
        # types of one block must be pairwise NON-isomorphic (needed by the symmetrised evaluator only)
        tys = [colex_decode(keys[k]["sigma"], s) for k in ks]
        for i in range(len(tys)):
            for j in range(i + 1, len(tys)):
                assert not any({frozenset(g[v] for v in e) for e in tys[i]} == tys[j] for g in permutations(range(s)))
        T = np.full(1 << nE, -1, dtype=np.int64)
        V = np.full(1 << nsig, -1, dtype=np.int64)
        for sx, k in sig_of.items():
            V[sx] = k
        for x in range(1 << nE):
            r = x & ((1 << nsig) - 1)
            if r not in sig_of:
                continue
            k = sig_of[r]
            c = canon_flag(dec(x, order), s, m, order)
            T[x] = canon_lists[k].get(c, -1)
        conf = factorial(7) // factorial(7 - s) * comb(7 - s, m - s) * comb(7 - m, m - s)
        assert 5040 % conf == 0
        blob_blocks.append((s, m, nsig, nE, conf, ks, order, V, T))

    # ---------------------------------------------------------------- objective constants
    # S(H) = 5040 D (slack - b) = 12 D F - T - Zsum,  F = 24 e - P,  slack = F/420 - Zsum/(5040 D)
    c12 = 12 * D
    Tt = 5040 * D * target
    assert Tt.denominator == 1
    Tt = int(Tt)
    # accumulation bound: sum_b (5040/conf_b) * conf_b * max|limb| < 2^62
    acc_bound = sum(5040 for _ in blob_blocks) * HALF
    assert acc_bound < (1 << 62)
    rep["int64_acc_bound_bits"] = acc_bound.bit_length()
    c12L, TtL = limbs_of(c12, NL), limbs_of(Tt, NL)
    assert all(abs(v) * 840 < (1 << 62) for v in c12L + TtL)

    write_blob(a.out, blob_blocks, keys, Gs, NL, c12L, TtL)
    rep["c12_limbs"] = [str(v) for v in c12L]
    rep["T_limbs"] = [str(v) for v in TtL]
    if a.sym_out:
        # Gsym_k[a, c] = sum_{h in Aut(sigma_k)} G_k[pi_h a, pi_h c],  pi_h a = index of the flag a with roots relabelled
        # by h (non-roots fixed).  For a root set S with first valid ordering theta0, the valid orderings are theta0 o h
        # and flag(theta0 o h, U) = pi_{h^-1} flag(theta0, U), so the sum over valid orderings is one Gsym entry.
        Gsym, auts = [], {}
        for (s, m), ks in blocks:
            order = flag_edge_order(s, m)
            for k in ks:
                sig = colex_decode(keys[k]["sigma"], s)
                aut = [g for g in permutations(range(s)) if {frozenset(g[v] for v in e) for e in sig} == sig]
                auts[k] = len(aut)
                pis = []
                for g in aut:
                    gg = list(g) + list(range(s, m))
                    pi = [canon_all[k][canon_flag({frozenset(gg[v] for v in e) for e in es}, s, m, order)]
                          for es in flags_all[k]]
                    assert sorted(pi) == list(range(len(pi)))
                    pis.append(np.array(pi))
                G = Gs[k]
                acc = np.zeros(G.shape, dtype=object)
                for pi in pis:
                    acc = acc + G[np.ix_(pi, pi)]
                assert (acc == acc.T).all()
                Gsym.append(acc)
        gmax2 = max([abs(int(x)) for G in Gsym for x in G.ravel()] + [1])
        NL2 = 1
        while (gmax2 >> (LB * (NL2 - 1))) >= HALF:
            NL2 += 1
        assert NL2 == NL, (NL2, NL)
        write_blob(a.sym_out, blob_blocks, keys, Gsym, NL, c12L, TtL, magic=0x53594D31)
        rep["sym"] = {"out": a.sym_out, "aut_sizes": auts, "Gsym_max_bits": gmax2.bit_length()}
    json.dump(rep, open(a.out + ".json", "w"), indent=1)
    for info in rep["keys"]:
        print(json.dumps(info))
    print("D", rep["D"], "limbs", NL, "G_max_bits", rep["G_max_bits"], "unused", rep.get("npz_unused_fields"),
          "sym", rep.get("sym"))


def write_blob(path, blob_blocks, keys, Gs, NL, c12L, TtL, magic=0x52375347):
    with open(path, "wb") as f:
        f.write(np.array([magic, len(blob_blocks), len(keys), NL, LB], dtype=np.int64).tobytes())
        f.write(np.array(c12L, dtype=np.int64).tobytes())
        f.write(np.array(TtL, dtype=np.int64).tobytes())
        for (s, m, nsig, nE, conf, ks, order, V, T) in blob_blocks:
            f.write(np.array([s, m, nsig, nE, conf, len(ks)], dtype=np.int64).tobytes())
            f.write(np.array(order, dtype=np.int64).ravel().tobytes())
            f.write(V.tobytes())
            f.write(T.tobytes())
        for k, G in enumerate(Gs):
            n = G.shape[0]
            nz = int(any(G.ravel()))
            f.write(np.array([k, n, nz], dtype=np.int64).tobytes())
            L = np.zeros((n, n, NL), dtype=np.int64)
            flat = G.ravel()
            for idx in range(flat.size):
                x = int(flat[idx])
                if x:
                    ls = limbs_of(x, NL)
                    assert sum(v << (LB * j) for j, v in enumerate(ls)) == x
                    assert all(abs(v) <= HALF for v in ls)
                    L[idx // n, idx % n, :] = ls
            f.write(L.tobytes())


if __name__ == "__main__":
    main()
