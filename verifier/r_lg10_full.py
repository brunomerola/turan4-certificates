"""Independent checker: C_10 over ALL 604,426 labelled R in G_9 (no isomorphism reduction).

For every labelled Giraud system R on range(9) the exact join of r_lg.py lists all links L of vertex 9 such that
every 7-set through 9 induces a Giraud system (so H = R + L is locally Giraud).  Every Giraud extension of R is
locally Giraud, hence among the survivors.  The number of Giraud extensions of R with (unique, CF3) partition sizes
(a, b) is 2^(a-1) + 2^(b-1) for b >= 1 and 2 for b = 0 (new vertex joins A with a new row modulo complementation,
or joins B with a new column modulo complementation; the two kinds have different partitions, unique at 10
vertices by Lemma 3).  So  #survivors(R) == 2^(a-1) + 2^(b-1)  for every R  proves that every survivor is a Giraud
system, i.e. C_10.  The count formula itself is checked by explicit construction (r_lg.direct_giraud_extensions,
set equality with the survivors) on the 33 representatives (results_lg_n10.json) and on a random sample of 2,000
labelled R here, where every survivor is also tested with st_lib.is_giraud_struct.
Output: results_lg_n10_full.json
"""
from __future__ import annotations

import json
import os
import random
import sys
import time
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.getcwd()   # outputs (results_*.json, local/) go to the working directory
import r_lg  # noqa: E402
import st_lib as L  # noqa: E402


def main():
    t0 = time.time()
    G7 = r_lg.load_set("G7.txt")
    tabs, _ = r_lg.build_tables(9, G7)
    fam9 = L.enum_giraud(9, normalised=True)            # mask -> set of partitions (unique: CF3)
    assert len(fam9) == 604426 and all(len(p) == 1 for p in fam9.values())
    bad = []
    tot = 0
    by_sizes = {}
    for i, (R, parts) in enumerate(sorted(fam9.items())):
        (A,) = tuple(parts)
        a = len(A)
        b = 9 - a
        want = 2 if b == 0 else (1 << (a - 1)) + (1 << (b - 1))
        S = r_lg.extensions(R, 9, tabs)
        k = len(S)
        if k != want or len(set(S)) != k:
            bad.append({"R": str(R), "survivors": k, "expected": want})
        tot += k
        key = f"{min(a, b)},{max(a, b)}"
        by_sizes[key] = by_sizes.get(key, 0) + 1
        if i % 50000 == 0:
            print(i, k, want, tot, round(time.time() - t0, 1), flush=True)
    # |G_10| from the definition (unique partitions at 10 vertices): 1 + sum over A containing 0, B nonempty
    g10 = 1 + sum(comb(9, a - 1) * (1 << ((a - 1) * (10 - a - 1))) for a in range(1, 10))
    # sample: explicit construction and the structural Giraud test
    rng = random.Random(10_20261006)
    smp = rng.sample(sorted(fam9), 2000)
    smp_ok = True
    for R in smp:
        S = set(r_lg.extensions(R, 9, tabs))
        smp_ok &= S == r_lg.direct_giraud_extensions(R, 9)
        smp_ok &= all(L.is_giraud_struct(h, 10) for h in S)
    res = {"labelled_R": len(fam9), "R_by_part_sizes": by_sizes, "survivors_total": tot,
           "labelled_giraud_on_10_vertices_by_formula": g10, "total_equals_G10": tot == g10,
           "R_with_count_mismatch": bad[:20], "n_mismatch": len(bad),
           "sample_2000_explicit_and_struct_ok": smp_ok,
           "C10_FULL_PASS": bool(not bad and tot == g10 and smp_ok), "time_s": round(time.time() - t0, 1)}
    json.dump(res, open(os.path.join(OUT, "results_lg_n10_full.json"), "w"), indent=1)
    print(json.dumps(res))


if __name__ == "__main__":
    main()
