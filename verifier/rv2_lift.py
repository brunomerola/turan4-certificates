"""Independent checker: exact check of the lifting identity (rv1_lift.check_graph, used unchanged) on the
heaviest K support graphs of each of the six headline dual points plus NR random 7-vertex graphs.  The identity
Mmix = T Mtop T^T (every mixed or smaller flag-size block of a type of size s is a congruence image of the top block
(s, mt(s)), mt(s) = max{m : 2m - s <= 7}) is a statement about one graph H, independent of the predicate.
Usage: python rv2_lift.py K NR"""
import json
import os
import random
import sys

import rv1_lib as L
from rv1_lift import check_graph

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certificates", "n7_dual_points")
TAGS = ["p5no56_W", "p5no56_Wodd", "p6_f3", "p7_f1", "h44_nt3", "k6m_W"]
K, NR = int(sys.argv[1]), int(sys.argv[2])
graphs = []
for t in TAGS:
    D = json.load(open(f"{RES}/dual_{t}.json"))
    sup = sorted(D["support"], key=lambda x: -int(x["num"]))[:K]
    graphs += [(t, L.colex_to_own(int(x["mask"]))) for x in sup]
rng = random.Random(20261006 + 2)
graphs += [("random", rng.getrandbits(35)) for _ in range(NR)]
tot, ok_all = 0, True
per = {}
for src, h in graphs:
    ok, info = check_graph(h)
    if not ok:
        print("MISMATCH", src, h, info, flush=True)
        ok_all = False
    else:
        tot += info
        per[src] = per.get(src, 0) + 1
print(f"graphs checked {len(graphs)} {per}; (graph, s, ma, mb, labelled type) identities checked {tot}; "
      f"all equal: {ok_all}")
print("tops (s, mt(s)):", [(s, (7 + s) // 2) for s in range(6)], "; s = 6, 7 only allow m = s (1x1)")
