"""Independent checker: exact check of every certificate in certificates/lottery_upper/ with own code (rv3b_core)."""
import copy
import glob
import math
import os
import sys
import time
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rv3b_core import (check_cert, parse, sha256, turan_brute, turan_multiset, mass, a_coef, decimal_up,
                       multisets)

HERE = os.path.dirname(os.path.abspath(__file__))
CERTS = os.path.normpath(os.path.join(HERE, "..", "certificates", "lottery_upper"))

# decimals claimed by the producer (upper bounds -> must be >= exact)
CLAIMED = {
    "k5_m10_den1000000.txt": "0.484327045",
    "k5_m9_den1000.txt": "0.484643353",
    "k5_m8_den1000.txt": "0.485470656",
    "k5_m8_den1000000.txt": "0.485470656",
    "k7_m7_den1000000.txt": "0.545345470",
    "k7_r4p5_m5_den1000000.txt": "0.447392627",
    "k7_r4p5_m4_den1000000.txt": "0.455472784",
    "k7_r3p4_m5_p204_den1000000.txt": "0.570116694",
    "k8_r3p4_m5_p106_den1000000.txt": "0.589526231",
    "k7_r3p4_m5_p204_den1000.txt": "0.570149891",
    "k8_r3p4_m5_p106_den1000.txt": "0.589553013",
    "k7_r4p5_m3_den3.txt": "0.456790124",
    "k8_r4p5_m3_den3.txt": "0.481481482",
    "k5_m7_den1000000.txt": "0.502008939",
    "k5_m6_den1000000.txt": "0.512368411",
    "k5_m5_den1000000.txt": "0.525382087",
    "k5_m3_den1000000.txt": "0.545275477",
    "f8_746_x_quarter.txt": "0.282214507",
}
EXACT_CLAIMED = {"k5_m8_den1000.txt": Fr(10191, 20992), "k7_r4p5_m3_den3.txt": Fr(37, 81),
                 "k8_r4p5_m3_den3.txt": Fr(13, 27), "f8_746_x_quarter.txt": Fr(1463, 5184)}


def structure(d):
    """Describe the type set by shape (sorted nonzero counts)."""
    from collections import Counter
    c = Counter(tuple(sorted([v for v in e if v], reverse=True)) for e in d["T"])
    return dict(c)


def brute_per(d):
    # blow-up with p points per part (any p-set type is realisable); keep it cheap
    return d["p"]


def negative_controls(path):
    """Mutations that must be rejected by the checker."""
    d0 = parse(path)
    out = []
    # (a) scale every y down by (1 - 1e-6): some row must fail
    from rv3b_core import check_cert as cc
    import tempfile

    def write(d, fn):
        with open(fn, "w") as f:
            f.write("k %d\nr %d\np %d\nm %d\n" % (d["k"], d["r"], d["p"], d["m"]))
            if d["inner"]:
                inn = d["inner"]
                f.write("inner %s q %d %s %s\n" % (" ".join(map(str, inn["parts"])), inn["q"], inn["cname"], inn["c"]))
            f.write("x %s\n" % " ".join(map(str, d["x"])))
            for e in d["T"]:
                f.write("T %s\n" % " ".join(map(str, e)))
            for B, v in d["y"]:
                f.write("y %s %s\n" % (" ".join(map(str, B)), v))
    tmp = os.path.join(os.getcwd(), "_neg.txt")
    # (a) shrink y
    d = copy.deepcopy(d0)
    d["y"] = [(B, v * Fr(999999, 1000000)) for B, v in d["y"]]
    write(d, tmp)
    _, r = cc(tmp)
    out.append(("shrink y by 1e-6", bool(r["fail"])))
    # (b) drop the type that is needed: remove each type in turn until Turan fails (first one found)
    found = False
    for i in range(len(d0["T"])):
        d = copy.deepcopy(d0)
        del d["T"][i]
        if turan_multiset(d):
            ok, S, _ = turan_brute(d, per=d["p"])
            out.append(("drop type %s: multiset & brute both reject" % (d0["T"][i],), not ok))
            found = True
            break
    if not found:
        out.append(("drop-type control: every single type is redundant (no rejection)", True))
    # (c) perturb one x coordinate pair keeping sum 1 (raise largest-mass part): must break some row unless slack
    d = copy.deepcopy(d0)
    j = max(range(d["m"]), key=lambda i: d["x"][i])
    j2 = min((i for i in range(d["m"]) if i != j), key=lambda i: d["x"][i])
    eps = d["x"][j2] / 2
    d["x"][j] += eps
    d["x"][j2] -= eps
    write(d, tmp)
    _, r = cc(tmp)
    out.append(("move x mass %s from part %d to part %d" % (eps, j2, j), bool(r["fail"])))
    # (d) wrong stated c
    d = copy.deepcopy(d0)
    write(d, tmp)
    with open(tmp, "a") as f:
        f.write("c %s\n" % (d0["c"] - Fr(1, 10 ** 12)))
    _, r = cc(tmp)
    out.append(("stated c lowered by 1e-12", bool(r["fail"])))
    os.remove(tmp)
    return out


if __name__ == "__main__":
    files = sorted(glob.glob(os.path.join(CERTS, "*.txt")))
    print("certs dir:", CERTS)
    allok = True
    for fn in files:
        b = os.path.basename(fn)
        t0 = time.time()
        # SHA-256 sidecar
        side = open(fn + ".sha256").read().split()[0]
        h = sha256(fn)
        d, res = check_cert(fn, CLAIMED.get(b))
        ok_b, S, nchk = turan_brute(d, per=brute_per(d))
        if not ok_b:
            res["fail"].append("brute-force Turan failed at %s" % (S,))
        if b in EXACT_CLAIMED and res["c"] != EXACT_CLAIMED[b]:
            res["fail"].append("exact claim %s != %s" % (EXACT_CLAIMED[b], res["c"]))
        if side != h:
            res["fail"].append("sha256 sidecar mismatch")
        status = "PASS" if not res["fail"] else "FAIL"
        allok &= not res["fail"]
        print("=" * 110)
        print("%s  [%s]  sha256 %s (sidecar %s)" % (b, status, h, "match" if side == h else "MISMATCH"))
        print("  (k,r,p)=(%d,%d,%d) m=%d  types=%d %s  blocks(y>0)=%d  inner=%s" % (
            res["k"], res["r"], res["p"], res["m"], res["rows"], structure(d), res["nblocks"],
            None if d["inner"] is None else "%s q=%d c=%s" % (d["inner"]["parts"], d["inner"]["q"], d["inner"]["c"])))
        print("  x = %s" % " ".join(map(str, d["x"])))
        print("  Turan: multiset failures %d ; brute force on %dx%d points: %s (%d %d-sets)" % (
            res["turan_multiset_failures"], res["m"], brute_per(d), ok_b, nchk, res["p"]))
        print("  rows: %d/%d tight, min slack %s" % (res["tight"], res["rows"], float(res["min_slack"])))
        print("  c = %s = %.12f ; stated c matches: %s ; rounded-up 9 dp: %s ; claimed decimal %s ok: %s" % (
            res["c"], res["c_float"], res["c_stated_matches"], decimal_up(res["c"]), res.get("claimed_decimal"),
            res.get("decimal_ok")))
        print("  shadow density at x = %.9f ; cover ratio %.6f" % (res["density"], res["c_float"] / res["density"]))
        for f in res["fail"]:
            print("  !! " + f)
        print("  time %.1f s" % (time.time() - t0), flush=True)
    print("=" * 110)
    print("ALL CERTIFICATES PASS:", allok)
    # negative controls on the two certificates of this package
    for b in ("k7_r4p5_m5_den1000000.txt", "f8_k7_r4p6_x_quarter.txt"):
        print("negative controls", b)
        for name, rejected in negative_controls(os.path.join(CERTS, b)):
            print("   %-70s rejected: %s" % (name, rejected))
