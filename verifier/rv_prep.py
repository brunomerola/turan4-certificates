"""Independent checker: build evaluation tables for an N = 7 flag-algebra certificate (no producer imports).

Reads only the certificate FILES (RUN.cert.json for the key list (s, m, sigma, flags) and RUN.cert.npz for the integer
factor rows A_k).  Everything else is re-derived here from definitions:
  * own edge order on [n]: lexicographic order of 4-subsets, BUT on a flag vertex set [m] the root-only 4-sets
    (subsets of [s]) come first (lex), then the remaining ones (lex).  Producer masks use colex; they are decoded with
    an own colex rank function, so any encoding mismatch shows up as a disagreement in the cross-checks.
  * flag classes: canonical form = min over permutations of the non-root vertices of the own mask.
  * Qint_k = A_k^T A_k in exact Python integers (PSD by construction), checked to fit int64.
Writes a binary blob read by rv_eval.c:
  header, then per block: s, m, conf, nkeys, key list; table T[x] for every labelled m-vertex mask x (own order) ->
  global flag id or -1; offsets of every key's Q; then the Q arrays (int64, row-major).
Usage: python rv_prep.py RUN_PREFIX OUT.bin
"""
import json
import sys
from itertools import combinations, permutations
from math import comb, factorial

import numpy as np

R = 4


def colex_rank(e):
    return sum(comb(v, i + 1) for i, v in enumerate(sorted(e)))


def colex_decode(mask, n):
    """Producer mask (colex bit order on [n]) -> set of frozensets."""
    out = set()
    for e in combinations(range(n), R):
        if (mask >> colex_rank(e)) & 1:
            out.add(frozenset(e))
    return out


def flag_edge_order(s, m):
    allE = list(combinations(range(m), R))
    root = [e for e in allE if set(e) <= set(range(s))]
    rest = [e for e in allE if not set(e) <= set(range(s))]
    return root + rest


def enc(es, order):
    idx = {frozenset(e): i for i, e in enumerate(order)}
    x = 0
    for e in es:
        x |= 1 << idx[frozenset(e)]
    return x


def dec(x, order):
    return {frozenset(order[i]) for i in range(len(order)) if (x >> i) & 1}


def canon_flag(es, s, m, order):
    best = None
    for pu in permutations(range(s, m)):
        g = list(range(s)) + list(pu)
        x = enc({frozenset(g[v] for v in e) for e in es}, order)
        if best is None or x < best:
            best = x
    return best


LAM = int(__import__("os").environ.get("RV_LAM", "1"))   # admissible: every p-set spans >= LAM edges


def admissible(es, n, p):
    return all(sum(frozenset(f) in es for f in combinations(P, R)) >= LAM for P in combinations(range(n), p))


def main():
    prefix, out = sys.argv[1], sys.argv[2]
    cert = json.load(open(prefix + ".cert.json"))
    Z = np.load(prefix + ".cert.npz")
    p = cert["p"]
    keys = cert["keys"]
    assert cert["tau1_num"] == 0 and cert["tau2_num"] == 0, "stationarity terms present: not handled here"
    M = int(cert["M"])
    blocks = []
    for k, key in enumerate(keys):
        sm = (key["s"], key["m"])
        if not blocks or blocks[-1][0] != sm:
            assert all(b[0] != sm for b in blocks)
            blocks.append((sm, []))
        blocks[-1][1].append(k)
    report = {"p": p, "M": M, "blocks": []}
    gid_base = 0
    qs, qoffs = [], []
    blob_blocks = []
    for (s, m), ks in blocks:
        assert 2 * m - s <= 7, (s, m)
        order = flag_edge_order(s, m)
        nE = len(order)
        nsig = comb(s, R)
        sig_of = {}
        flagidx = {}
        binfo = {"s": s, "m": m, "keys": []}
        for k in ks:
            key = keys[k]
            sig_es = colex_decode(key["sigma"], s)
            sx = enc(sig_es, order)          # root-only edges are the first nsig bits of the own order
            assert sx < (1 << nsig) or nsig == 0 and sx == 0
            assert sx not in sig_of, "two keys with the same labelled type"
            sig_of[sx] = k
            # Qint
            A = Z[f"A{k}"]
            n = len(key["flags"])
            assert A.ndim == 2 and A.shape[1] == n and A.dtype == np.int64
            Ao = A.astype(object)
            Q = Ao.T @ Ao if A.shape[0] else np.zeros((n, n), dtype=object)
            qmax = max([abs(int(v)) for v in Q.ravel()] + [0])
            assert qmax < 2 ** 62
            # flags: decode, check root part == sigma, canonical class
            canon = []
            for f in key["flags"]:
                es = colex_decode(f, m)
                assert {e for e in es if e <= frozenset(range(s))} == sig_es, "flag root part != sigma"
                canon.append(canon_flag(es, s, m, order))
            dup = len(canon) - len(set(canon))
            # completeness of the flag list (informational; not needed for soundness)
            allc = set()
            for fb in range(1 << (nE - nsig)):
                es = sig_es | dec(fb << nsig, order)
                if admissible(es, m, p):
                    allc.add(canon_flag(es, s, m, order))
            missing = len(allc - set(canon))
            extra_inadm = len(set(canon) - allc)
            for i, c in enumerate(canon):
                flagidx.setdefault((k, c), []).append(i)
            binfo["keys"].append({"k": k, "sigma_colex": key["sigma"], "n_flags": n, "rows": int(A.shape[0]),
                                  "dup_classes": dup, "admissible_classes": len(allc), "missing": missing,
                                  "inadmissible_listed": extra_inadm, "Q_zero": qmax == 0, "qmax": qmax})
            # If a class is listed twice (dup), the vector entry is duplicated; we sum over all listed copies by
            # folding Q: v^T Q v with v having equal entries i, j  ==  merge rows/cols i and j.  Only done if needed.
            if dup:
                raise SystemExit("duplicate flag classes: not expected; handle explicitly")
            qoffs.append((k, n))
            qs.append(np.array(Q, dtype=np.int64))
        # table: labelled m-vertex mask -> global flag id
        T = np.full(1 << nE, -1, dtype=np.int32)
        V = np.full(1 << nsig, -1, dtype=np.int32)
        for sx, k in sig_of.items():
            V[sx] = k
        canon_cache = {}
        for x in range(1 << nE):
            r = x & ((1 << nsig) - 1)
            if r not in sig_of:
                continue
            k = sig_of[r]
            es = dec(x, order)
            c = canon_flag(es, s, m, order)
            ids = flagidx.get((k, c))
            if ids:
                T[x] = ids[0]
        blob_blocks.append((s, m, nsig, nE, order, ks, T, V))
        report["blocks"].append(binfo)
    # write blob
    with open(out, "wb") as f:
        hdr = np.array([p, len(blob_blocks), len(keys), M], dtype=np.int64)
        f.write(hdr.tobytes())
        for (s, m, nsig, nE, order, ks, T, V) in blob_blocks:
            conf = factorial(7) // factorial(7 - s) * comb(7 - s, m - s) * comb(7 - m, m - s)
            f.write(np.array([s, m, nsig, nE, conf, len(ks)], dtype=np.int64).tobytes())
            # edge order on [m]: nE x 4 vertex tuples
            f.write(np.array(order, dtype=np.int64).ravel().tobytes())
            f.write(V.astype(np.int64).tobytes())
            f.write(T.astype(np.int64).tobytes())
        for (k, n), Q in zip(qoffs, qs):
            f.write(np.array([k, n, int(np.any(Q != 0))], dtype=np.int64).tobytes())
            f.write(Q.astype(np.int64).ravel().tobytes())
    json.dump(report, open(out + ".json", "w"), indent=1)
    for b in report["blocks"]:
        print(b["s"], b["m"], [(x["k"], x["sigma_colex"], x["n_flags"], x["rows"], x["admissible_classes"],
                                x["missing"], x["inadmissible_listed"], x["Q_zero"]) for x in b["keys"]])


if __name__ == "__main__":
    main()
