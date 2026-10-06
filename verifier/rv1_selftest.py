"""Independent checker: self-tests of rv1_lib (own code): encodings, flag universes vs the paper's Table (p = 5), timing."""
import time
from itertools import combinations

import rv1_lib as L

# colex round trip and convention: bit i of a producer mask = i-th 4-subset of {0..6} in colex order
colex = sorted(combinations(range(7), 4), key=lambda e: e[::-1])
for i, e in enumerate(colex):
    assert L.colex_rank(e) == i
    h = L.colex_to_own(1 << i)
    assert L.edges_of(h) == [e]
    assert L.own_to_colex(h) == 1 << i
print("colex convention OK (first 15 colex 4-sets are the subsets of {0..5}:",
      all(max(e) <= 5 for e in colex[:15]), ")")

# flag universes for one representative of each type class, p = 5; paper Table 2 row t(5,4):
# (1,4) 2, (2,4) 2, (3,4) 2, (3,5) 23, (4,5) 15/16, (5,6) 809/854/904/960/1024 for 1..5-edge types
exp = {(1, 4): {0: 2}, (2, 4): {0: 2}, (3, 4): {0: 2}, (3, 5): {0: 23}, (4, 5): {0: 15, 1: 16},
       (5, 6): {0: 0, 1: 809, 2: 854, 3: 904, 4: 960, 5: 1024}}
for (s, m), byE in exp.items():
    got = {}
    for tau in L.all_labelled_types(s):
        got.setdefault(L.popcount(tau), set()).add(len(L.flag_universe(s, m, tau)))
    print((s, m), {k: sorted(v) for k, v in got.items()}, "expected", byE,
          all(got[k] == {v} for k, v in byE.items()))

# conf numbers
print("conf:", {sm: L.conf(*sm) for sm in exp}, "(paper: 140, 1260, 2520, 1260, 5040, 5040)")

# timing of the moment counts on the complete graph and a check that every count matrix is symmetric
h = (1 << 35) - 1
t = time.time()
for sm in exp:
    mc = L.moment_counts(h, *sm)
    for tau, d in mc.items():
        assert all(d.get((b, a)) == c for (a, b), c in d.items())
print(f"moment counts of K_7^(4) for all six blocks: {time.time() - t:.2f}s; symmetric")
# K_7: every theta has the complete type; one flag per block; count = conf
for sm in exp:
    mc = L.moment_counts(h, *sm)
    assert len(mc) == 1 and sum(mc[list(mc)[0]].values()) == L.conf(*sm) and len(mc[list(mc)[0]]) == 1
print("K_7 sanity OK")
