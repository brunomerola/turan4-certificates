"""Independent checker: compare s8_pycheck values with s8_eval list-mode per-mask values (exact integers S).
usage: s8_cmp.py PY.json EACH.txt CERT.cert.json OUT.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction


def main():
    pyp, eachp, certp, outp = sys.argv[1:5]
    P = json.load(open(pyp))
    b = Fraction(json.load(open(certp))["bound"])
    E = {}
    for line in open(eachp):
        m, st, F, Z, S = line.split()
        E[m] = (int(st), int(S))
    rows = []
    ok = True
    for m, v in P.items():
        st, S = E[m]
        same = st == 0 and int(v["S"]) == S
        ok &= same
        sl = Fraction(v["slack"])
        ok &= sl >= b
        rows.append({"mask": int(m), "equal": same, "slack_minus_b": float(sl - b), "at_min": sl == b})
    out = {"n": len(rows), "all_equal": ok, "rows": rows}
    json.dump(out, open(outp, "w"), indent=1)
    print(json.dumps({"n": len(rows), "all_equal_and_ge_b": ok, "n_at_min": sum(r["at_min"] for r in rows)}))


if __name__ == "__main__":
    main()
