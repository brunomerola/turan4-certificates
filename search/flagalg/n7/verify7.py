"""Independent sample verifier of an N = 7 certificate written by certify7.py (RUN.cert.json + RUN.cert.npz).

Uses neither the engine (../flags.py, ../hypergraphs.py) nor the n7 tables/kernels; pure Python + exact integers:
 1. for every key (s, m, sigma) the sigma-flags are re-enumerated here (all labelled m-vertex 4-graphs extending
    sigma, every p-subset contains an edge, canonical form = min colex mask over permutations fixing the roots)
    and must coincide with the flag lists rebuilt from the certificate's key list by this script's own code;
 2. Qint = A^T A from the integer factor rows (PSD by construction), exact Python integers;
 3. for sampled H (the argmin graph of every (e, cdd) group of the certificate, plus random admissible graphs) the
    slack is recomputed from the definition: an explicit loop over all ordered root tuples theta with H[theta] = sigma
    and all ordered pairs (U1, U2) of disjoint (m-s)-sets, M_sigma entries as exact Fractions;
 4. checks a_H/(1 + beta_H) >= claimed bound for every sampled H, with equality at the certificate's argmin graph,
    and that every argmin graph's recomputed Z equals the certificate's group maximum.
The completeness of the scan over all admissible H is NOT re-checked here (it is checked by certify7.py's two
independent scans: nauty classes and raw one-vertex extensions).
Usage: python verify7.py RUN_PREFIX [--random 100]
"""
import argparse
import json
import random
import time
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial

import numpy as np

R = 4


def colex(e):
    return sum(comb(v, i + 1) for i, v in enumerate(sorted(e)))


def elist(n):
    return sorted(combinations(range(n), R), key=lambda e: e[::-1])


def admissible(es, n, p):
    return all(any(frozenset(f) in es for f in combinations(P, R)) for P in combinations(range(n), p))


def canon_fix(es, m, s):
    best = None
    for pu in permutations(range(s, m)):
        g = tuple(range(s)) + pu
        mask = 0
        for e in es:
            mask |= 1 << colex(tuple(g[v] for v in e))
        if best is None or mask < best:
            best = mask
    return best


def flags_for(s, m, sigma, p):
    sig = {frozenset(e) for i, e in enumerate(elist(s)) if (sigma >> i) & 1}
    free = [e for e in elist(m) if not set(e) <= set(range(s))]
    out = set()
    for bits in range(1 << len(free)):
        es = sig | {frozenset(free[i]) for i in range(len(free)) if (bits >> i) & 1}
        if admissible(es, m, p):
            out.add(canon_fix(es, m, s))
    return sorted(out)


def graph_edges(h):
    return {frozenset(e) for i, e in enumerate(elist(7)) if (h >> i) & 1}


def inner_exact(h, keys, flags, Qi, p):
    """sum over keys of <Qint, M_sigma(H)> * (denominator-free): returns Fraction = sum <Q, M>  * M^2."""
    g = graph_edges(h)
    tot = Fraction(0)
    by_sm = {}
    for k, key in enumerate(keys):
        by_sm.setdefault((key["s"], key["m"]), []).append(k)
    for (s, m), ks in by_sm.items():
        sig_of = {key_sigma: k for k in ks for key_sigma in [keys[k]["sigma"]]}
        index = {k: {f: i for i, f in enumerate(flags[k])} for k in ks}
        acc = {k: 0 for k in ks}
        nconf = factorial(7) // factorial(7 - s) * comb(7 - s, m - s) * comb(7 - m, m - s)
        for th in permutations(range(7), s):
            sig = 0
            for i, e in enumerate(elist(s)):
                if frozenset(th[x] for x in e) in g:
                    sig |= 1 << i
            if sig not in sig_of:
                continue
            k = sig_of[sig]
            rest = [v for v in range(7) if v not in th]
            for U1 in combinations(rest, m - s):
                rest2 = [v for v in rest if v not in U1]
                for U2 in combinations(rest2, m - s):
                    fl = []
                    for U in (U1, U2):
                        order = th + U
                        es = [frozenset(e) for e in combinations(range(m), R) if frozenset(order[x] for x in e) in g]
                        fl.append(index[k][canon_fix(es, m, s)])
                    acc[k] += int(Qi[k][fl[0], fl[1]])
        for k in ks:
            tot += Fraction(acc[k], nconf)
    return tot


def cdd_exact(h):
    g = graph_edges(h)
    c = 0
    for v in range(7):
        O = [u for u in range(7) if u != v]
        for U in combinations(O, 3):
            Uc = tuple(u for u in O if u not in U)
            c += (frozenset((v,) + U) in g) and (frozenset((v,) + Uc) in g)
    return Fraction(c, 140)         # [[d.d]]_1 at N = 7: c counts ordered (v, U1, U2), 7 * 20 = 140 configurations


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix")
    ap.add_argument("--random", type=int, default=50)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    t0 = time.time()
    cert = json.load(open(a.prefix + ".cert.json"))
    Z = np.load(a.prefix + ".cert.npz")
    p, M2 = cert["p"], cert["M"] ** 2
    keys = cert["keys"]
    flags = []
    for k, key in enumerate(keys):
        fl = flags_for(key["s"], key["m"], key["sigma"], p)
        assert fl == key["flags"], f"key {k}: flag list differs from the certificate's"
        flags.append(fl)
    print(f"flags re-enumerated ({time.time() - t0:.0f}s): {[len(f) for f in flags]}", flush=True)
    Qi = []
    for k in range(len(keys)):
        A = Z[f"A{k}"].astype(object)
        Qi.append(A.T @ A if A.shape[0] else np.zeros((len(flags[k]),) * 2, dtype=object))
    claim = Fraction(cert["bound"])
    t1, t2 = Fraction(cert["tau1_num"], M2), Fraction(cert["tau2_num"], M2)
    sample = [(int(g[3]), (g[0], g[1]), int(g[2])) for g in cert["groups"]]
    rnd = random.Random(a.seed)
    # random admissible graphs (rejection sampling at random densities)
    while len(sample) < len(cert["groups"]) + a.random:
        dens = rnd.uniform(0.05, 0.6)
        h = sum(1 << i for i in range(35) if rnd.random() < dens)
        if admissible(graph_edges(h), 7, p):
            sample.append((h, None, None))
    worst = None
    for h, grp, zmax in sample:
        S = inner_exact(h, keys, flags, Qi, p)          # = sum <Qint, M>  (Q = Qint / M^2)
        e = bin(h).count("1")
        d = Fraction(e, 35)
        c = cdd_exact(h)
        if grp is not None:
            assert (e, int(c * 70)) == tuple(grp), (h, grp)
            assert S * 5040 == zmax, (h, S * 5040, zmax)           # kernel's Z equals the definition
        aH = d - S / M2 + t1 * c + t2 * (d - c)
        beta = t1 * d + t2 * (1 - d)
        val = aH / (1 + beta)
        assert val >= claim, (h, float(val), float(claim))
        worst = val if worst is None or val < worst else worst
    print(f"verified {len(sample)} graphs ({len(cert['groups'])} group argmins + {a.random} random): "
          f"min a/(1+beta) = {float(worst):.12f}; claim = {float(claim):.12f}; attained: {worst == claim} "
          f"({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
