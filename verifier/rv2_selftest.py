"""Independent checker: self-tests of rv2_lib (own code): universes, canonical forms, the three exact PSD criteria."""
import random
import time
from fractions import Fraction

import numpy as np

import rv1_lib as L
import rv2_lib as R

rng = random.Random(7)
nrng = np.random.default_rng(7)

# 1. universes: p = 5, lam = 1 must reproduce rv1_lib / Table 2; other predicates printed
for (s, m) in [(1, 4), (2, 4), (3, 4), (3, 5), (4, 5), (5, 6)]:
    for tau in L.all_labelled_types(s):
        assert R.flag_universe_pl(s, m, tau, 5, 1) == L.flag_universe(s, m, tau), (s, m, tau)
print("flag_universe_pl(p=5, lam=1) == rv1_lib.flag_universe for every labelled type of all six blocks")
for (p, lam) in [(5, 1), (6, 1), (7, 1), (5, 2), (6, 2)]:
    row = {}
    for (s, m) in [(1, 4), (2, 4), (3, 4), (3, 5), (4, 5), (5, 6)]:
        sz = {}
        for tau in L.all_labelled_types(s):
            sz.setdefault(L.popcount(tau), set()).add(len(R.flag_universe_pl(s, m, tau, p, lam)))
        row[(s, m)] = {k: sorted(v) for k, v in sz.items()}
    print(f"universe sizes p={p} lam={lam}:", row)

# 2. admissibility / canonical forms against brute force
for _ in range(300):
    h = rng.getrandbits(35)
    for (p, lam) in [(5, 1), (6, 1), (7, 1), (5, 2), (6, 2)]:
        bf = all(sum(1 for e in L.edges_of(h) if set(e) <= set(P)) >= lam
                 for P in __import__("itertools").combinations(range(7), p))
        assert R.admissible(h, p, lam) == bf
    assert R.admissible(h, 5, 1) == L.admissible5(h)
pm = L.perm_edge_maps()
hs = [rng.getrandbits(35) for _ in range(60)]
hs += [L.canonical_graph(h, pm) for h in hs[:10]]
perm = list(range(7))
for h in hs[:10]:
    rng.shuffle(perm)
    hs.append(sum(1 << L.IX[((perm[a] * 7 + perm[b]) * 7 + perm[c]) * 7 + perm[d]] for (a, b, c, d) in L.edges_of(h)))
cf = R.canonical_forms(hs)
assert cf == [L.canonical_graph(h, pm) for h in hs]
assert cf[:10] == cf[60:70] == cf[70:80]
print("canonical_forms == rv1_lib.canonical_graph on 80 graphs (incl. relabelled copies); admissible() == brute force")


# 3. PSD criteria on exact integer matrices with entries ~2^120 and tiny scaled eigenvalues
def make(n, lam_target, neg=False):
    Q, _ = np.linalg.qr(nrng.standard_normal((n, n)))
    ev = np.exp(nrng.uniform(np.log(lam_target), np.log(5.0), n))
    ev[0] = -lam_target if neg else lam_target
    M = (Q * ev) @ Q.T
    d = 2.0 ** nrng.integers(0, 31, n)                  # wildly different diagonal scales (entries up to ~2^122)
    A = np.zeros((n, n), dtype=object)
    for i in range(n):
        for j in range(i, n):
            v = int(Fraction(M[i, j]) * (1 << 60)) * int(d[i]) * int(d[j])
            A[i, j] = A[j, i] = v
    return A


for n, lt in [(8, 1e-3), (30, 1e-8), (60, 1e-11), (120, 1e-12), (300, 1e-12)]:
    for neg in (False, True):
        A = make(n, lt, neg)
        t = time.time()
        rr = R.pd_residual(A)
        tr = time.time() - t
        bb = R.pd_bareiss(A)[0] if n <= 160 else None
        ll = L.psd_exact([[A[i, j] for j in range(n)] for i in range(n)])[0] if n <= 64 else None
        wit = R.neg_witness(A)[0] if neg else None
        print(f"n={n:>3} lam~{lt:.0e} {'INDEFINITE' if neg else 'PD        '}: residual {rr['PASS']} "
              f"(worst {rr.get('worst_ratio')}, {tr:.1f}s)  Bareiss {bb}  LDL-Fraction {ll}  neg_witness {wit}")
        if neg:
            assert not rr["PASS"] and bb in (False, None) and ll in (False, None) and wit
        else:
            assert rr["PASS"] and bb in (True, None) and ll in (True, None)
print("self-tests OK")
