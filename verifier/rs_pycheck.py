"""Independent checker: pure-Python exact slack f(H) - sum_k <Q_k, M_k(H)> for a few graphs, independent of rs_prep.py /
rs_eval.c: Q entries as Fractions computed ON DEMAND from the certificate file,
    Q_k[a, c] = sum_i W[i,a] W[i,c] / 2^(2 kb) + sum_{i,j} B[i,a] (Xq[i][j] / q) B[j,c];
flag identity by explicit isomorphism search (bijections fixing the roots) against the decoded flag list; every
ordered root tuple theta with H[theta] == sigma (labelled); every ordered pair of disjoint (m-s)-sets; f from the
definition (explicit (T, x, y) loop).
Usage: python rs_pycheck.py CERT_NPZ KEYS_JSON TARGET MASK_COLEX [MASK_COLEX ...]
"""
import json
import sys
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial

import numpy as np

from rs_lib import colex_decode, f_value


def same_flag(A, B, s, m):
    for pu in permutations(range(s, m)):
        g = list(range(s)) + list(pu)
        if {frozenset(g[v] for v in e) for e in A} == B:
            return True
    return False


def main():
    npz, kjs, target = sys.argv[1], sys.argv[2], Fraction(sys.argv[3])
    z = np.load(npz)
    kb, q = int(z["kb"]), int(str(z["q"]))
    corr = {c["key"]: c for c in json.loads(str(z["corr"]))}
    keys = json.load(open(kjs))["keys"]
    W, Bm, Xq = [], {}, {}
    for k in range(len(keys)):
        Wr = z[f"W{k}"]
        W.append([[int(x) for x in row] for row in Wr] if Wr.size else [])
        if k in corr:
            Bm[k] = z[f"B{k}"].astype(int).tolist()
            Xq[k] = [[int(x) for x in row] for row in corr[k]["X_times_q"]]
    cache = {}

    def Q(k, a, c):
        key = (k, a, c)          # no symmetry assumed
        if key not in cache:
            v = Fraction(sum(r[a] * r[c] for r in W[k]), 1 << (2 * kb))
            if k in Bm:
                B, X = Bm[k], Xq[k]
                v += Fraction(sum(B[i][a] * X[i][j] * B[j][c] for i in range(len(B)) if B[i][a]
                                  for j in range(len(B)) if B[j][c]), q)
            cache[key] = v
        return cache[key]

    dec_keys = []
    for key in keys:
        s, m = key["s"], key["m"]
        dec_keys.append((s, m, colex_decode(key["sigma"], s), [colex_decode(f, m) for f in key["flags"]]))
    for hm in map(int, sys.argv[4:]):
        H = colex_decode(hm, 7)
        tot = f_value(H, 7)
        for k, (s, m, sig, fl) in enumerate(dec_keys):
            if not W[k] and k not in Bm:
                continue
            conf = factorial(7) // factorial(7 - s) * comb(7 - s, m - s) * comb(7 - m, m - s)
            acc = Fraction(0)
            for th in permutations(range(7), s):
                if {e for e in map(frozenset, combinations(range(s), 4)) if frozenset(th[v] for v in e) in H} != sig:
                    continue
                rest = [v for v in range(7) if v not in th]
                idx = {}
                for U in combinations(rest, m - s):
                    lab = list(th) + list(U)
                    F = {e for e in map(frozenset, combinations(range(m), 4)) if frozenset(lab[v] for v in e) in H}
                    idx[U] = next((i for i, f in enumerate(fl) if same_flag(F, f, s, m)), None)
                for U1 in idx:
                    for U2 in idx:
                        if set(U1) & set(U2) or idx[U1] is None or idx[U2] is None:
                            continue
                        acc += Q(k, idx[U1], idx[U2])
            tot -= acc / conf
        print("PY", hm, "slack", tot, "slack-b", tot - target, "S", (tot - target) * 5040 * int(str(z["D"])),
              flush=True)


if __name__ == "__main__":
    main()
