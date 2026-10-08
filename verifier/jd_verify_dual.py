#!/usr/bin/env python3
"""Independent checker: exact verifier of an N = 7 dual point for sigma(J_4) (reviewer code; imports
only the reviewer's jw_lib / jd_lib; reads the producer's JSON as data).

Usage: python -B jd_verify_dual.py DUAL.json [--json OUT] [--bareiss-all] [--crosscheck] [--relabel-check K]
Checks:
  D1 y = num_H / den: integers, every num_H >= 1, sum num_H = den; support classes distinct (own canonical form);
     every support graph J_4-free (own link-K_4 test).
  D2 V = sum_H y_H f(H), f = 1 - sigma_7 (own codegree formula, cross-checked by the 4-set formula) exactly; equals
     the file's value when present; 41/57 - V and 1 - V reported exactly.
  D3 for every block of the 9 and every CANONICAL labelled type (min lex mask over S_s): Y = sum num_H M(H) (integer
     configuration counts); every nonzero entry in O x O (O = flags with positive diagonal; a zero diagonal with a
     nonzero row is a failure); Y_OO PSD exactly: n <= 40 or --bareiss-all: fraction-free symmetric elimination
     (exact PSD criterion, singular allowed); n > 40: Sylvester (all leading principal minors > 0, computed exactly
     by multimodular elimination + CRT with a Hadamard bound), falling back to the exact elimination if a leading
     minor is 0.
  --crosscheck: for the smallest n > 40 key, exact Bareiss minors == multimodular minors.
  --relabel-check K: on K support graphs, every NON-canonical labelled type's count matrix equals the canonical one
     with flags relabelled by the root permutation (the relabel lemma: Y_{sigma^pi} = P Y_sigma P^T).
Exit 0 = all pass; 1 = a check fails.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import Counter
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import jd_lib as R  # noqa: E402

L = R.L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dual")
    ap.add_argument("--json", default=None)
    ap.add_argument("--bareiss-all", action="store_true")
    ap.add_argument("--crosscheck", action="store_true")
    ap.add_argument("--relabel-check", type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    D = json.load(open(a.dual))
    rep = {"dual": os.path.abspath(a.dual), "note": D.get("note")}
    fails = []

    def check(k, ok, info=""):
        rep[k] = bool(ok)
        print(f"  [{'ok' if ok else 'FAIL'}] {k} {info}", flush=True)
        if not ok:
            fails.append(k)

    sup = D["support"]
    den = int(D["den"])
    masks = [L.colex_to_lex(int(e["mask"])) for e in sup]
    nums = [int(e["num"]) for e in sup]
    rep["support"] = len(sup)
    rep["den_is_2^110"] = den == 2 ** 110
    print(f"# {a.dual}: {len(sup)} support graphs, den = {den} (2^110: {den == 2 ** 110})", flush=True)
    check("D1_nums_positive", all(x >= 1 for x in nums))
    check("D1_sum_equals_den", sum(nums) == den)
    can = L.canon_many(masks)
    check("D1_classes_distinct", len(set(can)) == len(can), f"({len(set(can))} distinct classes)")
    j4 = [H for H in masks if L.contains_J4(H)]
    check("D1_all_J4_free", not j4, f"({len(j4)} contain J_4)")
    # D2
    fv = [L.fval(H) for H in masks]
    V = sum(Fraction(x) * f for x, f in zip(nums, fv)) / den
    rep["V"] = str(V)
    rep["V_float"] = float(V)
    if "value" in D:
        check("D2_V_equals_file", V == Fraction(D["value"]), f"V = {float(V):.14f}")
    gap = Fraction(41, 57) - V
    rep["41/57-V"] = str(gap)
    rep["41/57-V_float"] = float(gap)
    rep["1-V"] = str(1 - V)
    rep["1-V_float"] = float(1 - V)
    rep["1-V > 0.280703275315"] = (1 - V) > Fraction(280703275315, 10 ** 12)
    print(f"# V = {V} = {float(V):.14f}; 41/57 - V = {float(gap):.10e}; 1 - V = {float(1 - V):.13f}", flush=True)
    if j4 or sum(nums) != den:
        return finish(rep, fails, a, t0)
    # D3
    Y = {}
    for H, x in zip(masks, nums):
        for key, cnt in R.moment_counts_canon(H).items():
            Yk = Y.setdefault(key, {})
            for fp, c in cnt.items():
                Yk[fp] = Yk.get(fp, 0) + x * c
    print(f"# moments: {len(Y)} (block, canonical type) keys ({time.time() - t0:.1f}s)", flush=True)
    rep["n_keys"] = len(Y)
    per = []
    nfail = 0
    cross_done = not a.crosscheck
    for key in sorted(Y, key=lambda k: (k[0], k[1], k[2])):
        Yk = Y[key]
        O = sorted({F for (F, G), v in Yk.items() if F == G and v > 0})
        Os = set(O)
        bad_rows = any((F not in Os or G not in Os) for (F, G), v in Yk.items() if v)
        n = len(O)
        rec = {"block": list(key[:2]), "type_lex": key[2], "type_edges": bin(key[2]).count("1"), "n": n}
        if bad_rows:
            rec.update(PSD=False, method="zero diagonal with nonzero row")
            per.append(rec)
            nfail += 1
            print(f"   FAIL {key}: zero diagonal with nonzero row", flush=True)
            continue
        ix = {F: i for i, F in enumerate(O)}
        A = [[0] * n for _ in range(n)]
        for (F, G), v in Yk.items():
            A[ix[F]][ix[G]] = v
        t1 = time.time()
        if n <= 40 or a.bareiss_all:
            psd, rk, piv, zero, why = L.psd_bareiss(A)
            rec.update(PSD=psd, method="exact fraction-free elimination", rank=rk, why=why)
        else:
            mins, info = R.leading_minors_exact(A)
            first_np = next((i for i, d in enumerate(mins) if d <= 0), None)
            if first_np is None:
                rec.update(PSD=True, method="Sylvester: all leading minors > 0 (multimodular exact)", rank=n,
                           **info, min_log2_pivot=round(min(
                               (mins[i].bit_length() - (mins[i - 1].bit_length() if i else 0)) for i in range(n))))
            elif mins[first_np] < 0:
                rec.update(PSD=False, method="Sylvester: negative leading minor after positive ones",
                           at=first_np, **info)
            else:
                psd, rk, piv, zero, why = L.psd_bareiss(A)
                rec.update(PSD=psd, method="zero leading minor -> exact fraction-free elimination", rank=rk,
                           why=why)
            if not cross_done:
                bm = R.bareiss_leading_minors(A)
                rec["crosscheck_bareiss_minors_equal"] = bm == mins[:len(bm)] and len(bm) == len(mins)
                check("D3_crosscheck_multimodular_vs_bareiss", rec["crosscheck_bareiss_minors_equal"],
                      f"(key {key}, n = {n})")
                cross_done = True
        rec["time_s"] = round(time.time() - t1, 2)
        per.append(rec)
        if not rec["PSD"]:
            nfail += 1
        print(f"   {key[:2]} type {key[2]:4d} ({rec['type_edges']} edges) n = {n:3d}: "
              f"{'PSD' if rec['PSD'] else 'NOT PSD'} [{rec['method']}] {rec['time_s']}s", flush=True)
    rep["blocks"] = per
    rep["max_n"] = max((r["n"] for r in per), default=0)
    check("D3_all_moment_matrices_PSD", nfail == 0, f"({nfail} of {len(per)} keys fail)")
    if a.relabel_check:
        relabel_check(rep, check, masks[:a.relabel_check])
    return finish(rep, fails, a, t0)


def relabel_check(rep, check, graphs):
    bad = 0
    nkeys = 0
    inv = lambda p: tuple(sorted(range(len(p)), key=lambda i: p[i]))  # noqa: E731
    for H in graphs:
        full = L.moment_counts(H)
        cano = R.moment_counts_canon(H)
        for (s, m, tau), cnt in full.items():
            sig = R.CANON_TYPE[s][tau]
            pi = next(p for p in R._PERMS_S[s] if R.relabel_mask(sig, s, inv(p)) == tau)
            ref = cano.get((s, m, sig), Counter())
            mapped = Counter()
            for (F1, F2), c in ref.items():
                mapped[(R.relabel_flag(sig, F1, s, m, pi), R.relabel_flag(sig, F2, s, m, pi))] += c
            nkeys += 1
            if mapped != cnt:
                bad += 1
    rep["relabel_keys_checked"] = nkeys
    check("D4_relabel_lemma_spotcheck", bad == 0, f"({nkeys} labelled (graph, block, type) keys, {bad} mismatches)")


def finish(rep, fails, a, t0):
    rep["fails"] = fails
    rep["PASS"] = not fails
    rep["time_s"] = round(time.time() - t0, 1)
    print(f"VERDICT: {'PASS' if not fails else 'FAIL'}; V = {rep.get('V_float')}; fails = {fails} ({rep['time_s']}s)",
          flush=True)
    if a.json:
        json.dump(rep, open(a.json, "w"), indent=1, default=str)
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
