"""Independent checker: compare an s7_eval scan with the certificate's claims (exact rationals).

usage: s7_analyze.py SCAN.json CERT.cert.json PREP.prep.json OUT.json [--list]
  min slack (exact), its multiplicity, the next distinct values, per-objective-value group maxima of Z against the
  producer's 'groups', number of admissible extensions against 'raw_extensions_scanned' (raw) or 'classes_scanned'
  (list), labelled count of admissible 7-vertex graphs (raw: sum over reps of (720/|Aut R|) * ext(R)),
  decimals: sigma <= 1 - b, ceiled at 12 digits.
"""
from __future__ import annotations

import json
import os
import sys
from decimal import ROUND_CEILING, ROUND_FLOOR, Decimal, getcontext
from fractions import Fraction

import numpy as np

getcontext().prec = 60
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.getcwd()  # s7_prep.py wrote local/ here


def dec(fr: Fraction, digits: int, mode) -> str:
    q = Decimal(1).scaleb(-digits)
    return str((Decimal(fr.numerator) / Decimal(fr.denominator)).quantize(q, rounding=mode))


def main():
    scan, certp, prepp, outp = sys.argv[1:5]
    S = json.load(open(scan))
    J = json.load(open(certp))
    P = json.load(open(prepp))
    den = int(S["S_den"])
    low = [(Fraction(int(v), den), int(m), int(c)) for v, m, c in S["lowest"]]
    bmin = low[0][0]
    claimed = Fraction(J["bound"])
    mine = {int(g[0]): (int(g[1]), int(g[2]), int(g[3])) for g in S["groups"]}
    theirs = {int(g[0]): (int(g[1]), int(g[2]), int(g[3])) for g in J["groups"]}
    same_keys = sorted(mine) == sorted(theirs)
    same_max = same_keys and all(mine[F][1] == theirs[F][0] for F in mine)
    diffs = [F for F in mine if F in theirs and mine[F][1] != theirs[F][0]]
    out = {"scan": os.path.basename(scan), "mode": S["mode"], "inputs": S["inputs"], "admissible": S["admissible"],
           "rejected": S["rejected"], "bad": S["bad"], "seconds": S["seconds"],
           "min_slack": str(bmin), "min_slack_float": float(bmin), "min_multiplicity": low[0][2],
           "argmin_mask_smallest": low[0][1],
           "claimed_bound": str(claimed), "min_equals_claimed": bmin == claimed,
           "next_values_minus_min": [float(v - bmin) for v, _, _ in low[1:8]],
           "next_masks": [m for _, m, _ in low[1:8]],
           "n_groups_mine": len(mine), "n_groups_cert": len(theirs), "group_keys_identical": same_keys,
           "group_maxZ_identical": same_max, "group_maxZ_diffs": diffs[:10]}
    if S["mode"] == "list":
        out["group_counts_identical"] = same_keys and all(mine[F][0] == theirs[F][2] for F in mine)
        out["classes_scanned_cert"] = J["classes_scanned"]
    else:
        out["raw_extensions_cert"] = J.get("raw_extensions_scanned")
        out["raw_extensions_equal"] = J.get("raw_extensions_scanned") == S["admissible"]
        prob = P["problem"]
        aut = np.load(os.path.join(OUT, "local", f"reps6_{prob}_aut.npy"))
        ext = S["ext"]
        assert len(ext) == aut.size
        out["labelled_admissible_7vertex"] = int(sum((720 // int(a)) * int(e) for a, e in zip(aut, ext)))
        out["reps6"] = int(aut.size)
        out["raw_reps6_cert"] = J.get("raw_reps6")
    one_minus = 1 - bmin
    out["sigma_upper_exact"] = str(one_minus)
    out["sigma_upper_ceil12"] = dec(one_minus, 12, ROUND_CEILING)
    out["sigma_upper_ceil10"] = dec(one_minus, 10, ROUND_CEILING)
    out["b_floor12"] = dec(bmin, 12, ROUND_FLOOR)
    out["PASS"] = bool(bmin == claimed and S["bad"] == 0 and same_max)
    json.dump(out, open(outp, "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
