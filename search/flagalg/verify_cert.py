"""Independent verifier of a Turan flag-algebra certificate written by turan.py --cert (predicate: every p-set
spans a number of edges in cert["allowed"]; default = at least one edge).

Does not use the engine's enumeration, flag or density code:
 1. admissible r-graphs on N vertices: ALL 2^C(N,r) labelled graphs are tested ("every p-set contains an edge"),
    and the isomorphism classes are obtained by a brute-force canonical form written here (pure Python);
 2. for every block: the sigma-flags are re-enumerated here and must coincide with the certificate's list;
    p(F_a, F_b; H) is recomputed by a direct loop over (theta, U1, U2);
 3. every Q_sigma is checked PSD exactly (symmetric elimination over the rationals);
 4. bound = min_H ( e(H)/C(N,r) - sum_sigma <Q_sigma, M_sigma(H)> ) in exact arithmetic, compared with the claim.
Feasible for C(N, r) <= 20 (N = 6, r = 3).  Usage: python verify_cert.py CERT.json
"""
import json
import sys
import time
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial


def colex(e):
    return sum(comb(v, i + 1) for i, v in enumerate(sorted(e)))


def edge_list(n, r):
    return sorted(combinations(range(n), r), key=lambda e: e[::-1])


def canon(edgeset, n, r, fixed=0):
    """min over permutations fixing 0..fixed-1 of the colex mask (edgeset = set of frozensets on range(n))."""
    best = None
    E = edge_list(n, r)
    for pu in permutations(range(fixed, n)):
        g = tuple(range(fixed)) + pu
        mask = 0
        for e in edgeset:
            mask |= 1 << colex(tuple(g[v] for v in e))
        if best is None or mask < best:
            best = mask
    return best


def is_psd_exact(Q):
    n = len(Q)
    A = [row[:] for row in Q]
    alive = list(range(n))
    while alive:
        piv = next((i for i in alive if A[i][i] > 0), None)
        if any(A[i][i] < 0 for i in alive):
            return False
        if piv is None:
            return all(A[i][j] == 0 for i in alive for j in alive)
        alive.remove(piv)
        p = A[piv][piv]
        for i in alive:
            f = A[i][piv] / p
            if f:
                for j in alive:
                    A[i][j] -= f * A[piv][j]
    return True


def main():
    cert = json.load(open(sys.argv[1]))
    r, p, N = cert["r"], cert["p"], cert["N"]
    claim = Fraction(cert["bound"])
    allowed = set(cert.get("allowed") or range(1, comb(p, r) + 1))   # edge counts allowed inside every p-set
    t0 = time.time()
    E = edge_list(N, r)
    pmasks = [sum(1 << i for i, e in enumerate(E) if set(e) <= set(P)) for P in combinations(range(N), p)]
    perm_maps = [[colex(tuple(g[v] for v in e)) for e in E] for g in permutations(range(N))]
    seen = bytearray(1 << len(E))
    classes = []
    for mask in range(1 << len(E)):          # orbit marking: one representative per isomorphism class
        if seen[mask] or not all(bin(mask & pm).count("1") in allowed for pm in pmasks):
            continue
        bits = [i for i in range(len(E)) if (mask >> i) & 1]
        for pm in perm_maps:
            seen[sum(1 << pm[i] for i in bits)] = 1
        classes.append([frozenset(E[i]) for i in bits])
    print(f"admissible classes on {N} vertices: {len(classes)} ({time.time() - t0:.1f}s)", flush=True)
    graphs = [set(v) for v in classes]
    slack = [Fraction(len(g), comb(N, r)) for g in graphs]
    for bi, blk in enumerate(cert["blocks"]):
        s, m, sigma = blk["s"], blk["m"], blk["sigma"]
        sig_edges = {frozenset(e) for i, e in enumerate(edge_list(s, r)) if (sigma >> i) & 1}
        # re-enumerate sigma-flags
        free = [e for e in edge_list(m, r) if not set(e) <= set(range(s))]
        fl = set()
        for bits in range(1 << len(free)):
            es = sig_edges | {frozenset(free[i]) for i in range(len(free)) if (bits >> i) & 1}
            ok = all(sum(frozenset(t) in es for t in combinations(P, r)) in allowed for P in combinations(range(m), p))
            if ok:
                fl.add(canon(es, m, r, fixed=s))
        assert sorted(fl) == sorted(blk["flags"]), f"block {bi}: flag list mismatch"
        index = {f: i for i, f in enumerate(blk["flags"])}
        den = Fraction(blk["den"])
        Q = [[Fraction(int(x)) / den for x in row] for row in blk["Qnum"]]
        assert all(Q[i][j] == Q[j][i] for i in range(len(Q)) for j in range(len(Q)))
        assert is_psd_exact(Q), f"block {bi} not PSD"
        total = factorial(N) // factorial(N - s) * comb(N - s, m - s) * comb(N - m, m - s)
        for gi, g in enumerate(graphs):
            acc = Fraction(0)
            for th in permutations(range(N), s):
                if {frozenset(th[i] for i in e) for e in combinations(range(s), r)} & g != \
                        {frozenset(th[i] for i in e) for e in sig_edges}:
                    continue
                rest = [v for v in range(N) if v not in th]
                for U1 in combinations(rest, m - s):
                    rest2 = [v for v in rest if v not in U1]
                    for U2 in combinations(rest2, m - s):
                        fs = []
                        for U in (U1, U2):
                            order = th + U
                            es = [frozenset(e) for e in combinations(range(m), r)
                                  if frozenset(order[i] for i in e) in g]
                            fs.append(index[canon(es, m, r, fixed=s)])
                        acc += Q[fs[0]][fs[1]]
            slack[gi] -= acc / total
        print(f"block {bi} (s={s}, m={m}, {len(blk['flags'])} flags): flags match, PSD ok ({time.time() - t0:.1f}s)",
              flush=True)
    b = min(slack)
    print(f"verified bound = {float(b):.12f}; claim = {float(claim):.12f}; equal: {b == claim}; "
          f"verified >= claim: {b >= claim}")


if __name__ == "__main__":
    main()
