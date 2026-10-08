"""Independent checker: the key-list comparison of s8_keys_cmp.py for the three certificates of the first review:
cosig3_s3_J4_r30 (sigma(J_4)), cosig3_s3_K5lt_f (sigma(K_5^<)), sigma4_c4_K5m_f (sigma(K_5^{4-})).
Reviewer's own definitions only (s8_lib, s8_prep.own_types_flags, helpers of s8_keys_cmp.py); keys.json, cert.json and
cert.npz read as data; reference hashes = those recorded by the first review (s7_review_hashes.txt).

usage: s8_keys_cmp_stage1.py OUT.json
Checks (1)-(6) as in s8_keys_cmp.py, plus predicate negative controls: the K_5^< spec replaced by K_5^='s, and
(p, lam) = (5, 2) replaced by (5, 1) -- both must be detected.
"""
from __future__ import annotations

import copy
import json
import os
import random
import sys
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s8_lib as L  # noqa: E402
from s8_keys_cmp import PKG, sha, own_keys, compare_keys, c4_direct, fmt_objective  # noqa: E402

CERTS = [  # (name, problem, package_keys subdir, original dir (relative to sigma_cat), prefix, objective)
    ("J4_r30", "s3_J4", "S3", "j4_dual_scope/S3/results", "cosig3_s3_J4_r30", "cosig3"),
    ("K5lt_f", "s3_K5lt", "cloud_a4", "cloud_a4/results_vm", "cosig3_s3_K5lt_f", "cosig3"),
    ("K5m4_f", "c4_K5m", "cloud_a4", "cloud_a4/results_vm", "sigma4_c4_K5m_f", "sigma4"),
]
DIRS = {"J4_r30": "J4", "K5lt_f": "K5lt", "K5m4_f": "K5_4minus"}   # directories under certificates/sigma_catalogue


def sig7_inputs():
    d = {}
    for line in open(os.path.join(HERE, "s7_review_hashes.txt")):
        if line.startswith("#") or "  " not in line:
            continue
        h, p = line.rstrip("\n").split("  ", 1)
        d[p] = h
    return d


def c4_masks():
    nr = np.random.default_rng(81)
    a7 = nr.integers(0, 1 << 35, 200000, dtype=np.int64)
    sp = np.zeros(40000, dtype=np.int64)
    for i in range(sp.size):
        for _ in range(int(nr.integers(0, 12))):
            sp[i] |= np.int64(1) << np.int64(nr.integers(0, 35))
    return np.arange(1 << 15, dtype=np.int64), np.concatenate([a7, sp, np.array([0, (1 << 35) - 1], dtype=np.int64)])


def pred_check(prob, pr, r):
    if r == 3:
        cp = set()
        for spec in pr["forbidden"]:
            cp |= set(L.copies(L.parse(spec), 7, 3))
        return cp == set(L.problem_copies(prob, 7)), f"forbidden {pr['forbidden']}; copies in K_7 {len(cp)}"
    p, lam = pr["p"], pr["lam"]
    a6, a7 = c4_masks()
    ok6 = np.array_equal(c4_direct(a6, 6, p, lam), L.admissible_np(prob, a6, 6)) if p <= 6 else \
        bool(np.all(L.admissible_np(prob, a6, 6)))
    ok7 = np.array_equal(c4_direct(a7, 7, p, lam), L.admissible_np(prob, a7, 7))
    return bool(ok6 and ok7), f"p = {p}, lam = {lam}; n = 6 all masks {ok6}, n = 7 sample {ok7}"


def main():
    outp = sys.argv[1]
    ref = sig7_inputs()
    rng = random.Random(82)
    res, controls = {}, {}
    for name, prob, pk, od, pre, obj in CERTS:
        r = L.PROBLEMS[prob][0]
        K = json.load(open(os.path.join(PKG, DIRS[name], pre + ".keys.json")))
        cj = os.path.join(PKG, DIRS[name], f"{pre}.{obj}.cert.json")
        J = json.load(open(cj))
        Z = np.load(cj.replace(".cert.json", ".cert.npz"), allow_pickle=False)
        c = {}
        hs = {}
        for suf in ("cert.json", "cert.npz", "verify.json"):
            fn = f"{pre}.{obj}.{suf}"
            hs[fn] = {"copy": sha(os.path.join(PKG, DIRS[name], fn)), "original": sha(os.path.join(PKG, DIRS[name], fn)),
                      "R7_SIG7": ref.get(f"tools/flagalg/sigma_cat/{od}/{fn}"),
                      "keys_json": K["copied_files_sha256"].get(fn)}
        c["hashes_equal"] = all(len(set(v.values())) == 1 and v["R7_SIG7"] is not None for v in hs.values())
        c["hashes"] = {k: v["copy"] for k, v in hs.items()}
        b = Fraction(J["bound"])
        c["scalars_equal"] = (K["problem"] == J["problem"] == prob and K["objective"] == J["objective"] == obj
                              and K["r"] == r and [list(x) for x in K["blocks"]] == [list(x) for x in J["blocks"]]
                              and K["M"] == J["M"] == 1 << K["M_bits"] and K["M_bits"] == J["M_bits"]
                              and Fraction(K["bound"]) == b and Fraction(K["sigma_upper"]) == 1 - b)
        c["M_bits"], c["blocks"] = K["M_bits"], K["blocks"]
        pr = K["predicate"]
        c["predicate_equal"], c["predicate_detail"] = pred_check(prob, pr, r)
        fm = K["format"]
        fmt_names = all(spec in fm for spec in pr["forbidden"]) if r == 3 else f">= {pr['lam']} edges" in fm
        mine = own_keys(prob, [tuple(x) for x in J["blocks"]])
        rows = [tuple(int(x) for x in Z[f"A{k}"].shape) for k in range(len(J["keys"]))]
        mism = compare_keys(K["keys"], mine, rows)
        c["keys_equal"] = not mism
        c["keys_vs_cert_json"] = [(k["s"], k["m"], k["type_index"], k["sigma"], k["n_flags"]) for k in J["keys"]] == \
            [(x["s"], x["m"], x["type_index"], x["sigma"], x["n_flags"]) for x in mine]
        c["n_keys"], c["total_flags"] = len(mine), sum(x["n_flags"] for x in mine)
        c["counts_equal"] = K["n_keys"] == len(mine) and K["total_flags"] == c["total_flags"]
        c["aut_sizes"] = sorted(set(x["aut_size"] for x in mine))
        c["mismatches"] = mism[:20]
        objs_ok = all(fmt_objective(obj, m) == L.objective_value(obj, m, 7)
                      for m in [L.rand_mask(7, r, rng.choice([0.2, 0.5, 0.8, 0.95]), rng) for _ in range(300)])
        c["format_objective_equals_definition"] = objs_ok
        c["format_ok"] = bool(f"M = 2^{K['M_bits']}" in fm and fmt_names and "colex" in fm and objs_ok
                              and ("F(H)/210" in fm if obj == "cosig3" else "F(H)/420" in fm))
        c["PASS"] = all(c[x] for x in ("hashes_equal", "scalars_equal", "predicate_equal", "keys_equal",
                                       "keys_vs_cert_json", "counts_equal", "format_ok"))
        res[name] = c
        print(name, {k: v for k, v in c.items() if k not in ("hashes", "mismatches")}, flush=True)
        # negative controls
        if name == "J4_r30":
            ki = next(i for i, x in enumerate(mine) if x["n_flags"] >= 2)
            t1 = copy.deepcopy(K["keys"]); t1[ki]["flags"][0], t1[ki]["flags"][1] = t1[ki]["flags"][1], t1[ki]["flags"][0]
            t2 = copy.deepcopy(K["keys"]); t2[ki]["aut_size"] += 1
            t3 = copy.deepcopy(K["keys"]); t3[ki]["flags"][-1] += 1 << 40
            t4 = copy.deepcopy(K["keys"]); t4[ki]["factor_rows"] += 1
            for lab, t in (("flags_swapped", t1), ("aut_changed", t2), ("flag_replaced", t3), ("rows_changed", t4)):
                controls[lab] = bool(compare_keys(t, mine, rows))
            controls["key_used"] = f"J4_r30 key {ki}"
        if name == "K5lt_f":
            bad = dict(pr, forbidden=[L.FORB["K5eq"].replace(",", " ")])
            controls["K5lt_spec_replaced_by_K5eq"] = not pred_check(prob, bad, r)[0]
        if name == "K5m4_f":
            controls["lam_2_to_1"] = not pred_check(prob, dict(pr, lam=1), r)[0]
    print("controls (True = detected)", controls, flush=True)
    out = {"certificates": res, "negative_controls_detected": controls,
           "ALL_PASS": all(v["PASS"] for v in res.values()) and all(v for k, v in controls.items() if k != "key_used")}
    json.dump(out, open(outp, "w"), indent=1)
    print("ALL_PASS", out["ALL_PASS"])


if __name__ == "__main__":
    main()
