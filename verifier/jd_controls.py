#!/usr/bin/env python3
"""Independent checker: controls for the verifier jd_verify_dual.py (reviewer code).
Positive: producer's ctrl_pos_Sn / ctrl_pos_Fanoc; reviewer-built exact laws of S_n (V = 25/27) and of the iterated
Fano complement phi (V = 41/57, singular -> --bareiss-all); the producer's second point dual_W1.
Negative: producer's ctrl_neg_lp / ctrl_neg_dead / ctrl_neg_j4; reviewer-built: B_n law (contains J_4), headline with
one numerator + 1 (not a probability vector), headline with a wrong claimed value, headline with 2^-40 moved onto a
J_4-free graph outside the support found by the reviewer's own search to touch a dead flag (zero diagonal of Y with a
nonzero off-diagonal entry), and (1/2) headline + (1/2) point mass on its heaviest graph (all flags stay live; must be
caught by the exact PSD test itself).  Exit 0 iff every control behaves as expected."""
from __future__ import annotations

import json
import os
import random
import subprocess
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import jd_lib as R  # noqa: E402

L = R.L
PY = sys.executable
S2R = os.path.normpath(os.path.join(HERE, "..", "certificates", "sigma_catalogue", "J4_limit"))
OUT = os.path.join(os.getcwd(), "out", "ctrl")
os.makedirs(OUT, exist_ok=True)


def law_json(name, note):
    law7, prob, desc = L.construction_law7(name)
    cl = L.classes_of(law7)
    den = L.lcm_den(cl.values())
    sup = [{"mask": L.lex_to_colex(H), "num": int(p * den)} for H, p in sorted(cl.items())]
    return {"problem": "J4", "den": den, "note": note, "support": sup}


def write(nm, D):
    p = os.path.join(OUT, nm + ".json")
    json.dump(D, open(p, "w"))
    return p


def run(path, *extra):
    js = os.path.join(OUT, os.path.basename(path).replace(".json", ".rv2.json"))
    r = subprocess.run([PY, "-B", os.path.join(HERE, "jd_verify_dual.py"), path, "--json", js, *extra],
                       capture_output=True, text=True)
    rep = json.load(open(js)) if os.path.exists(js) else {}
    if r.returncode not in (0, 1) or not rep:
        print(r.stdout[-2000:], r.stderr[-2000:])
    return r.returncode, rep


def main():
    W2 = json.load(open(os.path.join(S2R, "dual_W2.json")))
    den = int(W2["den"])
    sup = W2["support"]
    cases = []
    cases.append(("producer ctrl_pos_Sn", os.path.join(S2R, "ctrl_pos_Sn.json"), (), 0))
    cases.append(("producer ctrl_pos_Fanoc", os.path.join(S2R, "ctrl_pos_Fanoc.json"), ("--bareiss-all",), 0))
    cases.append(("producer dual_W1", os.path.join(S2R, "dual_W1.json"), (), 0))
    cases.append(("producer ctrl_neg_lp", os.path.join(S2R, "ctrl_neg_lp.json"), (), 1))
    cases.append(("producer ctrl_neg_dead", os.path.join(S2R, "ctrl_neg_dead.json"), (), 1))
    cases.append(("producer ctrl_neg_j4", os.path.join(S2R, "ctrl_neg_j4.json"), (), 1))
    cases.append(("rv S_n law", write("rv_pos_Sn", law_json("Sn", "reviewer S_n law")), (), 0))
    cases.append(("rv phi law", write("rv_pos_Fanoc", law_json("Fanoc", "reviewer phi law")), ("--bareiss-all",), 0))
    cases.append(("rv B_n law (J_4 inside)", write("rv_neg_Bn", law_json("Bn", "reviewer B_n law")), (), 1))
    s1 = [dict(e) for e in sup]
    s1[0]["num"] = int(s1[0]["num"]) + 1
    cases.append(("rv num+1", write("rv_neg_sum", {**W2, "support": s1, "note": "num+1"}), (), 1))
    cases.append(("rv wrong value", write("rv_neg_value", {**W2, "value": str(Fraction(W2["value"]) + Fraction(1, 2 ** 100)),
                                                         "note": "value + 2^-100"}), (), 1))
    # dead-flag graph: own search
    masks = [L.colex_to_lex(int(e["mask"])) for e in sup]
    nums = [int(e["num"]) for e in sup]
    Y = {}
    for H, x in zip(masks, nums):
        for key, cnt in R.moment_counts_canon(H).items():
            Yk = Y.setdefault(key, {})
            for fp, c in cnt.items():
                Yk[fp] = Yk.get(fp, 0) + x * c
    live = {k: {F for (F, G), v in Yk.items() if F == G and v > 0} for k, Yk in Y.items()}
    supc = set(L.canon_many(masks))
    rnd = random.Random(7)
    found = None
    tries = 0
    while found is None and tries < 20000:
        tries += 1
        h = masks[rnd.randrange(len(masks))] ^ (1 << rnd.randrange(35))     # one-triple neighbour of a support graph
        if L.contains_J4(h) or L.canon_many([h])[0] in supc:
            continue
        for key, cnt in R.moment_counts_canon(h).items():
            if key not in Y:
                continue
            for (F, G), c in cnt.items():
                if c and F != G and F not in live[key] and (F, F) not in cnt:
                    found = (h, key, F)
                    break
            if found:
                break
    print(f"dead-flag graph search: {tries} tries, found {found is not None}")
    eps = den >> 40
    s3 = [dict(e) for e in sup]
    i_max = max(range(len(s3)), key=lambda i: int(s3[i]["num"]))
    s3[i_max]["num"] = int(s3[i_max]["num"]) - eps
    s3.append({"mask": L.lex_to_colex(found[0]), "num": eps})
    D3 = {k: v for k, v in W2.items() if k != "value"}
    cases.append(("rv 2^-40 onto dead-flag graph", write("rv_neg_dead", {**D3, "support": s3, "note": f"dead flag {found[1:]}"}), (), 1))
    s4 = [dict(e) for e in sup]
    for e in s4:
        e["num"] = int(e["num"]) // 2
    s4[i_max]["num"] += den - sum(int(e["num"]) for e in s4)
    cases.append(("rv half mass on heaviest graph", write("rv_neg_half", {**D3, "support": s4, "note": "1/2 W2 + 1/2 point mass"}), (), 1))
    allok = True
    for nm, path, extra, ex in cases:
        rc, rep = run(path, *extra)
        ok = rc == ex
        allok &= ok
        bad = [(b["block"], b["type_lex"], b["n"], b["method"]) for b in rep.get("blocks", []) if not b.get("PSD", True)]
        print(f"{'CONTROL OK ' if ok else 'CONTROL BAD'} {nm:34s} exit {rc} (expected {ex}); V = {rep.get('V')}; "
              f"fails {rep.get('fails')}; non-PSD keys {len(bad)} {bad[:3]}", flush=True)
    print("ALL CONTROLS AS EXPECTED" if allok else "SOME CONTROL NOT AS EXPECTED")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
