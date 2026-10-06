"""Independent check of the N = 6 UPPER bound (dual point) written by n6_control.py.  Pure Python + Fractions;
imports nothing from the engine (flag enumeration, types, densities and the PSD test are re-implemented here).

Claim checked: y (rational, >= 0, sum 1) on admissible 6-vertex 4-graphs, with
  (i)  every graph in the support admissible (every p-set contains an edge),
  (ii) for EVERY type sigma on s vertices and every m with m > s, 2m - s <= 6 (all admissible sigma, all
       sigma-flags on m vertices, flags re-enumerated here), Y_sigma = sum_H y_H M_sigma(H) is PSD, where
       M_sigma(H)[a, b] = P[(H[theta+U1], theta) ~ F_a and (H[theta+U2], theta) ~ F_b] over uniform injective
       theta: [s] -> V and ordered disjoint (U1, U2),
  (iii) sum_H y_H e(H)/15 = claimed value.
Then by weak duality no plain flag-algebra certificate at N = 6 proves t(p, 4) > value.
Usage: python dual_check.py results/n6_control_p5.json
"""
import json
import sys
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial

N, R = 6, 4


def canon(es, m, s):
    best = None
    for pu in permutations(range(s, m)):
        g = tuple(range(s)) + pu
        key = tuple(sorted(tuple(sorted(g[v] for v in e)) for e in es))
        if best is None or key < best:
            best = key
    return best


LAM = 1


def admissible(es, verts, p):
    """every p-subset spans >= LAM edges (LAM from the json, default 1)"""
    return all(sum(frozenset(e) in es for e in combinations(P, R)) >= LAM for P in combinations(verts, p))


def is_psd(Q):
    n = len(Q)
    A = [row[:] for row in Q]
    alive = list(range(n))
    while alive:
        if any(A[i][i] < 0 for i in alive):
            return False
        piv = next((i for i in alive if A[i][i] > 0), None)
        if piv is None:
            return all(A[i][j] == 0 for i in alive for j in alive)
        alive.remove(piv)
        for i in alive:
            f = A[i][piv] / A[piv][piv]
            if f:
                for j in alive:
                    A[i][j] -= f * A[piv][j]
    return True


def main():
    global LAM
    d = json.load(open(sys.argv[1]))
    p = d["p"]
    LAM = int(d.get("lam", 1))
    sup = [(Fraction(x["y"]), {frozenset(e) for e in x["graph_edges"]}) for x in d["dual"]["support"]]
    assert all(y > 0 for y, _ in sup) and sum(y for y, _ in sup) == 1
    for _, g in sup:
        assert admissible(g, range(N), p)
    val = sum(y * Fraction(len(g), comb(N, R)) for y, g in sup)
    print(f"p={p} lam={LAM}: support {len(sup)}, sum y = 1, all admissible, value = {val}")
    ntypes = 0
    for s in range(0, N):
        for m in range(s + 1, N + 1):
            if 2 * m - s > N:
                continue
            E_s = list(combinations(range(s), R))
            types = set()
            for bits in range(1 << len(E_s)):
                es = {frozenset(E_s[i]) for i in range(len(E_s)) if (bits >> i) & 1}
                if admissible(es, range(s), p):
                    types.add(canon(es, s, 0))
            free = [e for e in combinations(range(m), R) if not set(e) <= set(range(s))]
            for sig in sorted(types):
                sig_es = {frozenset(e) for e in sig}
                flags = set()
                for bits in range(1 << len(free)):
                    es = sig_es | {frozenset(free[i]) for i in range(len(free)) if (bits >> i) & 1}
                    if admissible(es, range(m), p):
                        flags.add(canon(es, m, s))
                flags = sorted(flags)
                idx = {f: i for i, f in enumerate(flags)}
                n = len(flags)
                Y = [[Fraction(0)] * n for _ in range(n)]
                total = factorial(N) // factorial(N - s) * comb(N - s, m - s) * comb(N - m, m - s)
                for y, g in sup:
                    for th in permutations(range(N), s):
                        if {frozenset(th[i] for i in e) for e in E_s if frozenset(th[i] for i in e) in g} != \
                                {frozenset(th[i] for i in e) for e in sig}:
                            continue
                        rest = [v for v in range(N) if v not in th]
                        for U1 in combinations(rest, m - s):
                            rest2 = [v for v in rest if v not in U1]
                            for U2 in combinations(rest2, m - s):
                                fs = []
                                for U in (U1, U2):
                                    order = th + U
                                    es = [e for e in combinations(range(m), R)
                                          if frozenset(order[i] for i in e) in g]
                                    fs.append(idx[canon(es, m, s)])
                                Y[fs[0]][fs[1]] += y / total
                ok = is_psd(Y)
                ntypes += 1
                assert ok, f"Y not PSD for s={s} m={m} sigma={sig}"
            print(f"  (s,m)=({s},{m}): {len(types)} types, all moment matrices PSD", flush=True)
    print(f"PASS p={p}: plain flag algebra at N = 6 (all blocks 2m-s<=6, {ntypes} (type, m) pairs) cannot prove "
          f"t_{LAM}({p},4) > {val} = {float(val):.12f}")


if __name__ == "__main__":
    main()
