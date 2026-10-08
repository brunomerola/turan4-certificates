"""Independent checker: compare the producer's key lists (sigma_catalogue/<dir>/<prefix>.keys.json) with
the reviewer's OWN types and flags (s8_prep.own_types_flags, s8_lib definitions).  No producer module is imported;
keys.json, cert.json and cert.npz are read as data only.

usage: s8_keys_cmp.py OUT.json
Per certificate:
  (1) sha256 of the three certificate files = the reviewer's recorded values of the originals
      (s8_review_hashes.txt) = keys.json "copied_files_sha256";
  (2) keys.json scalars = cert.json (problem, objective, r, blocks, M, M_bits, bound) and sigma_upper = 1 - bound;
  (3) predicate: s3 -- the copy set in K_7 of the listed forbidden graph(s) (reviewer's copies()) = the reviewer's
      problem_copies(); c4 -- "every p-set spans >= lam edges" evaluated directly = the reviewer's admissible_np on all
      2^15 masks on 6 vertices and 200,000 random + 40,000 sparse masks on 7 vertices;
  (4) every key in matrix order: s, m, type_index, sigma, |Aut(sigma)| (reviewer: perms of [s] fixing the mask),
      n_flags, the FULL flag list in order, factor_rows = rows of cert.npz A<k> (reviewer's read), and
      A<k>.shape[1] = n_flags; n_keys and total_flags;
  (5) format string: names M = 2^M_bits and the forbidden graph / (p, lam); its objective formula (cosig3:
      F = 210 - sum_Q C(e_H(Q), 2), f = F/210; sigma4: F = 24 e(H) - sum_T c_T (c_T - 1), f = F/420) evaluated
      directly = the reviewer's definitional objective_value on 300 random 7-vertex graphs;
  (6) negative controls (in memory, first certificate): two flags swapped, |Aut| changed, a flag replaced, factor_rows
      changed -- each must be detected.
"""
from __future__ import annotations

import hashlib
import json
import os
import random
import sys
from fractions import Fraction
from itertools import combinations, permutations
from math import comb

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s8_lib as L  # noqa: E402
from s8_prep import own_types_flags  # noqa: E402

PKG = os.path.normpath(os.path.join(HERE, "..", "certificates", "sigma_catalogue"))   # the package's copies
DIRS = {"C5_f": "C5", "K5eq_A": "K5eq", "K5m_f": "K5_3minus", "K63_A": "K6_3", "K6m4_f": "K6_4minus",
        "K64_A": "K6_4", "K74_A": "K7_4"}
CERTS = [  # (name, problem, results dir, prefix, objective)
    ("C5_f", "s3_C5", "cloud_a5", "cosig3_s3_C5_f", "cosig3"),
    ("K5eq_A", "s3_K5eq", "cloud_a5", "cosig3_s3_K5eq", "cosig3"),
    ("K5m_f", "s3_K5m", "cloud_a6", "cosig3_s3_K5m_f", "cosig3"),
    ("K63_A", "s3_K6", "cloud_a6", "cosig3_s3_K6", "cosig3"),
    ("K6m4_f", "c4_K6m", "cloud_a5", "sigma4_c4_K6m_f", "sigma4"),
    ("K64_A", "c4_K6", "cloud_a5", "sigma4_c4_K6", "sigma4"),
    ("K74_A", "c4_K7", "cloud_a6", "sigma4_c4_K7", "sigma4"),
]


def sha(p):
    with open(p, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def my_inputs():
    d = {}
    for line in open(os.path.join(HERE, "s8_review_hashes.txt")):
        if line.startswith("#") or not line.strip():
            continue
        h, p = line.rstrip("\n").split("  ", 1)
        d[p] = h
    return d


def aut_size(sig, s, r):
    return sum(1 for g in permutations(range(s)) if L.apply_perm(sig, s, r, g) == sig)


def own_keys(prob, blocks):
    r = L.PROBLEMS[prob][0]
    out = []
    for (s, m) in blocks:
        types, flags = own_types_flags(prob, s, m)
        for i, (t, fl) in enumerate(zip(types, flags)):
            out.append({"s": s, "m": m, "type_index": i, "sigma": int(t), "aut_size": aut_size(int(t), s, r),
                        "n_flags": int(fl.size), "flags": [int(x) for x in fl]})
    return out


def compare_keys(theirs, mine, rows):
    """Returns a list of mismatch descriptions (empty = all equal)."""
    bad = []
    if len(theirs) != len(mine):
        return [f"n_keys {len(theirs)} != {len(mine)}"]
    for k, (a, b) in enumerate(zip(theirs, mine)):
        if a.get("k") != k:
            bad.append(f"key {k}: k field {a.get('k')}")
        for f in ("s", "m", "type_index", "sigma", "aut_size", "n_flags"):
            if a.get(f) != b[f]:
                bad.append(f"key {k}: {f} {a.get(f)} != {b[f]}")
        if list(a.get("flags", [])) != b["flags"]:
            bad.append(f"key {k}: flag list differs")
        if a.get("factor_rows") != rows[k][0]:
            bad.append(f"key {k}: factor_rows {a.get('factor_rows')} != {rows[k][0]}")
        if rows[k][1] != b["n_flags"]:
            bad.append(f"key {k}: A cols {rows[k][1]} != n_flags {b['n_flags']}")
    return bad


def popcount(x):
    x = x.astype(np.uint64)
    c = np.zeros(x.shape, dtype=np.int64)
    while np.any(x):
        c += (x & np.uint64(1)).astype(np.int64)
        x >>= np.uint64(1)
    return c


def c4_direct(masks, n, p, lam):
    """every p-subset of [n] spans >= lam edges (4-uniform), evaluated directly from the colex encoding."""
    ok = np.ones(masks.shape, dtype=bool)
    for P in combinations(range(n), p):
        sub = 0
        for e in combinations(P, 4):
            sub |= 1 << L.cidx(e)
        ok &= popcount(masks & np.int64(sub)) >= lam
    return ok


def fmt_objective(obj, mask):
    if obj == "cosig3":
        F = 210 - sum(comb(sum(1 for t in combinations(Q, 3) if (mask >> L.cidx(t)) & 1), 2)
                      for Q in combinations(range(7), 4))
        return Fraction(F, 210)
    e = bin(mask).count("1")
    P = 0
    for T in combinations(range(7), 3):
        c = sum(1 for x in range(7) if x not in T and (mask >> L.cidx(tuple(sorted(T + (x,))))) & 1)
        P += c * (c - 1)
    return Fraction(24 * e - P, 420)


def main():
    outp = sys.argv[1]
    mine_sha = my_inputs()
    rng = random.Random(81)
    res = {}
    controls = None
    for name, prob, d, pre, obj in CERTS:
        r = L.PROBLEMS[prob][0]
        K = json.load(open(os.path.join(PKG, DIRS[name], pre + ".keys.json")))
        cj = os.path.join(PKG, DIRS[name], f"{pre}.{obj}.cert.json")
        J = json.load(open(cj))
        Z = np.load(cj.replace(".cert.json", ".cert.npz"), allow_pickle=False)
        c = {}
        # (1) hashes
        hs = {}
        for suf in ("cert.json", "cert.npz", "verify.json"):
            fn = f"{pre}.{obj}.{suf}"
            rel = f"tools/flagalg/sigma_cat/{d}/results_vm/{fn}"
            hs[fn] = {"copy": sha(os.path.join(PKG, DIRS[name], fn)), "original": sha(os.path.join(PKG, DIRS[name],
                      fn)), "mine": mine_sha.get(rel), "keys_json": K["copied_files_sha256"].get(fn)}
        c["hashes_equal"] = all(len(set(v.values())) == 1 and v["mine"] is not None for v in hs.values())
        c["hashes"] = {k: v["copy"] for k, v in hs.items()}
        # (2) scalars
        b = Fraction(J["bound"])
        c["scalars_equal"] = (K["problem"] == J["problem"] == prob and K["objective"] == J["objective"] == obj
                              and K["r"] == r and [list(x) for x in K["blocks"]] == [list(x) for x in J["blocks"]]
                              and K["M"] == J["M"] == 1 << K["M_bits"] and K["M_bits"] == J["M_bits"]
                              and Fraction(K["bound"]) == b and Fraction(K["sigma_upper"]) == 1 - b)
        # (3) predicate
        pr = K["predicate"]
        if r == 3:
            cp = set()
            for spec in pr["forbidden"]:
                cp |= set(L.copies(L.parse(spec), 7, 3))
            c["predicate_equal"] = cp == set(L.problem_copies(prob, 7))
            c["predicate_detail"] = f"forbidden {pr['forbidden']}; copies in K_7 {len(cp)}"
            fmt_names = all(spec in K["format"] for spec in pr["forbidden"])
        else:
            p, lam = pr["p"], pr["lam"]
            a6 = np.arange(1 << 15, dtype=np.int64)
            ok6 = np.array_equal(c4_direct(a6, 6, p, lam), L.admissible_np(prob, a6, 6)) if p <= 6 else \
                bool(np.all(L.admissible_np(prob, a6, 6)))
            nr = np.random.default_rng(81)
            a7 = nr.integers(0, 1 << 35, 200000, dtype=np.int64)
            sp = np.zeros(40000, dtype=np.int64)
            for i in range(sp.size):
                for _ in range(int(nr.integers(0, 12))):
                    sp[i] |= np.int64(1) << np.int64(nr.integers(0, 35))
            a7 = np.concatenate([a7, sp, np.array([0, (1 << 35) - 1], dtype=np.int64)])
            ok7 = np.array_equal(c4_direct(a7, 7, p, lam), L.admissible_np(prob, a7, 7))
            c["predicate_equal"] = bool(ok6 and ok7)
            c["predicate_detail"] = f"p = {p}, lam = {lam}; n = 6 all masks {ok6}, n = 7 sample {ok7}"
            fmt_names = f">= {lam} edges" in K["format"]
        # (4) keys
        mine = own_keys(prob, [tuple(x) for x in J["blocks"]])
        rows = [tuple(int(x) for x in Z[f"A{k}"].shape) for k in range(len(J["keys"]))]
        mism = compare_keys(K["keys"], mine, rows)
        cert_keys = [(k["s"], k["m"], k["type_index"], k["sigma"], k["n_flags"]) for k in J["keys"]]
        c["keys_equal"] = not mism
        c["keys_vs_cert_json"] = cert_keys == [(x["s"], x["m"], x["type_index"], x["sigma"], x["n_flags"]) for x in mine]
        c["n_keys"] = len(mine)
        c["total_flags"] = sum(x["n_flags"] for x in mine)
        c["counts_equal"] = K["n_keys"] == len(mine) and K["total_flags"] == c["total_flags"]
        c["aut_sizes"] = sorted(set(x["aut_size"] for x in mine))
        c["mismatches"] = mism[:20]
        # (5) format string
        fm = K["format"]
        objs_ok = all(fmt_objective(obj, m) == L.objective_value(obj, m, 7)
                      for m in [L.rand_mask(7, r, rng.choice([0.2, 0.5, 0.8, 0.95]), rng) for _ in range(300)])
        c["format_ok"] = bool(f"M = 2^{K['M_bits']}" in fm and fmt_names and "colex" in fm and objs_ok
                              and ("F(H)/210" in fm if obj == "cosig3" else "F(H)/420" in fm))
        c["format_objective_equals_definition"] = objs_ok
        c["PASS"] = all(c[x] for x in ("hashes_equal", "scalars_equal", "predicate_equal", "keys_equal",
                                       "keys_vs_cert_json", "counts_equal", "format_ok"))
        res[name] = c
        print(name, {k: v for k, v in c.items() if k not in ("hashes", "mismatches")}, flush=True)
        # (6) negative controls once
        if controls is None:
            import copy
            controls = {}
            ki = next(i for i, x in enumerate(mine) if x["n_flags"] >= 2 and x["aut_size"] >= 1)
            t1 = copy.deepcopy(K["keys"]); t1[ki]["flags"][0], t1[ki]["flags"][1] = t1[ki]["flags"][1], t1[ki]["flags"][0]
            t2 = copy.deepcopy(K["keys"]); t2[ki]["aut_size"] += 1
            t3 = copy.deepcopy(K["keys"]); t3[ki]["flags"][-1] += 1 << 40
            t4 = copy.deepcopy(K["keys"]); t4[ki]["factor_rows"] += 1
            for lab, t in (("flags_swapped", t1), ("aut_changed", t2), ("flag_replaced", t3), ("rows_changed", t4)):
                controls[lab] = bool(compare_keys(t, mine, rows))
            controls["key_used"] = ki
            print("controls (True = detected)", controls, flush=True)
    out = {"certificates": res, "negative_controls_detected": controls,
           "ALL_PASS": all(v["PASS"] for v in res.values()) and all(v for k, v in controls.items() if k != "key_used")}
    json.dump(out, open(outp, "w"), indent=1)
    print("ALL_PASS", out["ALL_PASS"])


if __name__ == "__main__":
    main()
