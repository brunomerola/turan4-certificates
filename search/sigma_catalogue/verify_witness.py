#!/usr/bin/env python3
"""Self-contained EXACT verifier of a sharpness-obstruction witness for the plain N = 7 flag-algebra method
(3-graphs, the 9 blocks of the r = 3 pipeline).  Own code: Python fractions / integers and numpy integer arrays only;
no import from the flag-algebra engine or from any other file of this repository.  Producer code (not reviewed).

CLAIM checked for a witness file W = (construction c, problem F, rational vector mu on 7-vertex classes):
  Let phi be the exact 7-vertex law of construction c (recomputed here from its definition), t = E_phi[f] where
  f(H) = 1 - (1/210) sum_{4-sets Q of [7]} C(e_H(Q), 2) (= cosig3, the x != y codegree-squared objective).
  If every check below passes, then every plain N = 7 certificate (Q, b) -- Q_{beta,sigma} real symmetric PSD for
  every block beta in BLOCKS and every LABELLED type sigma, slack(H) = f(H) - sum <Q_{beta,sigma}, M_{beta,sigma}(H)>
  >= b for every F-free 7-vertex 3-graph H -- with <Q_{beta,sigma}, M_{beta,sigma}(phi)> = 0 for all (beta, sigma)
  ("Q kills range M(phi)") has  b <= E_mu[f] = t - delta.  If delta > 0, in particular NO certificate attains b = t.
Checks (each printed):
  C1  phi: exact recursion, total mass 1, t = E_phi f equals the construction's stated value; every class F-free.
  C2  mu >= 0, sum mu = 1, supp(mu) contained in supp(phi) (classes compared by min-over-S_7 canonical forms).
  C3  for every block beta and every labelled type sigma occurring in phi (types not occurring in phi do not occur in
      any class of supp(mu)): (a) M_{beta,sigma}(phi) restricted to its occurring flags O is PSD (exact LDL^T);
      (b) every nonzero entry of M_{beta,sigma}(H) for H in supp(mu) lies in O x O; (c) an exact basis X of
      ker M_{beta,sigma}(phi)|_O (from the LDL^T; M X = 0 re-checked); (d) X^T M_{beta,sigma}(mu)|_O X = 0 exactly.
  C4  delta = t - E_mu f, exact; must equal the file's value.
PROOF (given C1-C4): for (Q, b) as in the claim and each (beta, sigma): Q, M(phi) PSD and <Q, M(phi)> = 0 give
  Q M(phi) = 0; M(phi) lives on O x O, so the principal submatrix Q_OO satisfies Q_OO M_OO(phi) = 0, hence (symmetric)
  Q_OO = X Y X^T for a symmetric Y.  By (b), <Q, M(H)> = <Q_OO, M_OO(H)> for H in supp(mu), so
  sum_H mu_H <Q, M(H)> = <Y, X^T M_OO(mu) X> = 0 by (d).  As every H in supp(mu) is F-free (C1, C2),
  b <= min_{H in supp mu} slack(H) <= sum_H mu_H slack(H) = E_mu f - 0 = t - delta.            QED
  (M_{beta,sigma}(H) here = integer COUNTS of configurations; certificates use counts / (number of configurations of
  beta), a positive factor per block that Q absorbs, so the claim is unchanged.)
NOT claimed: anything about certificates that do not kill range M(phi) other than "b = t is impossible when
  delta > 0" (b = t forces <Q, M(phi)> = 0, since E_phi slack = t - sum <Q, M(phi)> >= b); in particular NOT that the
  N = 7 optimum is < t (that needs dual attainment, e.g. a strictly feasible pseudo-density, or an explicit dual
  point); nothing about stationarity/tau terms, other extra constraints, or N >= 8.
Definitions used for M: H a 3-graph on [7] (35-bit mask, colex 3-subsets: {a<b<c} -> a + C(b,2) + C(c,3)); block
  beta = (s, m); a configuration = (theta, U1, U2): theta an injective map [s] -> [7] (ordered tuple), U1, U2
  disjoint (m - s)-subsets of [7] - theta (ordered pair); its type = H[theta] relabelled theta_i -> i (labelled 3-graph on
  [s]); the flag of U = the 3-graph on [m] obtained from H[theta + U] (theta_i -> i, U -> s..m-1), up to permutations of
  s..m-1 (canonical code = min over them of the mask of the triples that meet {s..m-1}); M_{beta,sigma}(H)[F1, F2] =
  #{configurations with type sigma, flag(U1) = F1, flag(U2) = F2}.
Usage: python verify_witness.py WITNESS.json [--json OUT.json]
Exit code: 0 = all checks pass and delta > 0 (obstruction proved); 2 = all checks pass, delta = 0 (nothing proved);
1 = some check fails.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from fractions import Fraction
from math import comb

import numpy as np

BLOCKS = [(0, 3), (1, 3), (2, 3), (1, 4), (2, 4), (3, 4), (3, 5), (4, 5), (5, 6)]

# ------------------------------------------------------------------------------------------------ constructions
FANO = [(0, 1, 3), (1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 0), (5, 6, 1), (6, 0, 2)]     # lines {i, i+1, i+3} mod 7
S6 = [(0, 1, 2), (1, 2, 3), (2, 3, 4), (3, 4, 0), (4, 0, 1), (0, 2, 5), (1, 3, 5), (2, 4, 5), (1, 4, 5), (0, 3, 5)]
CONSTRUCTIONS = {   # name: (m parts of weight 1/m, pattern of part-triples, iterated?, problem, sigma = 1 - t)
    "Fanoc": (7, [t for t in itertools.combinations(range(7), 3) if set(t) not in [set(f) for f in FANO]], True,
              "J4", Fraction(16, 57)),
    "S6it": (6, S6, True, "K4m", Fraction(4, 43)),
    "edgeit": (3, [(0, 1, 2)], True, "C5m", Fraction(1, 13)),
    "Sn": (3, [(0, 1, 2)], False, "F5", Fraction(2, 27)),
    "Cn": (3, [(0, 1, 2), (0, 0, 1), (1, 1, 2), (0, 2, 2)], False, "K4", Fraction(1, 3)),
    "Bn": (2, [(0, 0, 1), (0, 1, 1)], False, "K5", Fraction(5, 8)),
}
FORBIDDEN = {       # BCL Table 2 edge lists, 1-based
    "J4": "123 124 125 134 135 145", "K4m": "123 124 134", "C5m": "123 234 345 145", "F5": "123 124 345",
    "K4": "123 124 134 234", "K5": "123 124 125 134 135 145 234 235 245 345",
}


def cidx(t) -> int:
    a, b, c = sorted(t)
    return a + b * (b - 1) // 2 + c * (c - 1) * (c - 2) // 6


def triples(k):
    return sorted(itertools.combinations(range(k), 3), key=cidx)


def _rgs(A):
    n, k = A.shape
    ids = np.zeros((n, k), dtype=np.int64)
    nxt = np.zeros(n, dtype=np.int64)
    for i in range(k):
        found = np.zeros(n, dtype=bool)
        for j in range(i):
            msk = (~found) & (A[:, j] == A[:, i])
            ids[msk, i] = ids[msk, j]
            found |= msk
        ids[~found, i] = nxt[~found]
        nxt[~found] += 1
    return (ids << (3 * np.arange(k))).sum(axis=1)


def _blocks_of(code, k):
    bl = {}
    for v in range(k):
        bl.setdefault((code >> (3 * v)) & 7, []).append(v)
    return [tuple(b) for _, b in sorted(bl.items())]


def law(name):
    """Exact law of 7 i.i.d. points of the construction: dict canonical mask -> Fraction."""
    m, pat, iterated, _, _ = CONSTRUCTIONS[name]
    P = {tuple(sorted(t)) for t in pat}
    tab = np.zeros((m, m, m), dtype=bool)
    for t in itertools.product(range(m), repeat=3):
        tab[t] = (tuple(sorted(t)) in P) and not (t[0] == t[1] == t[2])
    D, DEN = {1: {0: 1}, 2: {0: 1}}, {1: 1, 2: 1}
    for k in range(3, 8):
        A = np.array(list(itertools.product(range(m), repeat=k)), dtype=np.int64)
        if iterated:                                  # constant assignments recurse onto the same law
            A = A[~(A == A[:, :1]).all(axis=1)]
        FE = np.zeros(A.shape[0], dtype=np.int64)
        for t in triples(k):
            FE |= tab[A[:, t[0]], A[:, t[1]], A[:, t[2]]].astype(np.int64) << cidx(t)
        if not iterated:                              # a triple inside one part is a non-edge
            u, c = np.unique(FE, return_counts=True)
            D[k], DEN[k] = {int(a): int(b) for a, b in zip(u, c)}, m ** k
            continue
        base = 1
        for j in range(3, k):
            base *= DEN[j] ** (k // j)
        key = (_rgs(A) << 36) | FE
        u, c = np.unique(key, return_counts=True)
        tot = {}
        for kk, cnt in zip(u.tolist(), c.tolist()):
            bl = _blocks_of(kk >> 36, k)
            mult = base
            for B in bl:
                if len(B) >= 3:
                    mult //= DEN[len(B)]
            cur = {kk & ((1 << 36) - 1): cnt * mult}
            for B in bl:
                if len(B) < 3:
                    continue
                loc = [(cidx(t), cidx((B[t[0]], B[t[1]], B[t[2]]))) for t in triples(len(B))]
                emb = []
                for msk, num in D[len(B)].items():
                    g = 0
                    for li, gi in loc:
                        if (msk >> li) & 1:
                            g |= 1 << gi
                    emb.append((g, num))
                new = {}
                for g, w in cur.items():
                    for gm, num in emb:
                        new[g | gm] = new.get(g | gm, 0) + w * num
                cur = new
            for h, w in cur.items():
                tot[h] = tot.get(h, 0) + w
        D[k], DEN[k] = tot, (m ** k - m) * base
    assert sum(D[7].values()) == DEN[7]
    masks = np.array(sorted(D[7]), dtype=np.int64)
    can = canon(masks)
    out = {}
    for h, cc in zip(masks.tolist(), can.tolist()):
        out[cc] = out.get(cc, 0) + D[7][h]
    return {c: Fraction(v, DEN[7]) for c, v in out.items()}


_W = None


def canon(masks):
    """min over S_7 of the relabelled masks (exact int64 arithmetic; masks < 2^35)."""
    global _W
    if _W is None:
        T = triples(7)
        _W = np.zeros((5040, 35), dtype=np.int64)
        for gi, g in enumerate(itertools.permutations(range(7))):
            for t in T:
                _W[gi, cidx(t)] = 1 << cidx(tuple(g[v] for v in t))
    masks = np.asarray(masks, dtype=np.int64)
    out = np.empty(masks.size, dtype=np.int64)
    for i in range(0, masks.size, 1024):
        Bt = (masks[i:i + 1024, None] >> np.arange(35)[None, :]) & 1
        out[i:i + 1024] = (Bt @ _W.T).min(axis=1)
    return out


def F_value(h) -> Fraction:
    s = 0
    for Q in itertools.combinations(range(7), 4):
        e = sum((h >> cidx(t)) & 1 for t in itertools.combinations(Q, 3))
        s += e * (e - 1) // 2
    return Fraction(210 - s, 210)


def free_of(h, prob) -> bool:
    es = [tuple(int(ch) - 1 for ch in w) for w in FORBIDDEN[prob].split()]
    nv = 1 + max(max(e) for e in es)
    for g in itertools.permutations(range(7), nv):
        if all((h >> cidx((g[a], g[b], g[c]))) & 1 for a, b, c in es):
            return False
    return True


# ------------------------------------------------------------------------------------------------ moments
def _flag_tables(s, m):
    LT = [t for t in triples(m) if max(t) >= s]       # local triples meeting the free vertices s..m-1
    pos = {t: i for i, t in enumerate(LT)}
    perms = []
    for p in itertools.permutations(range(s, m)):
        mp = list(range(s)) + list(p)
        perms.append([pos[tuple(sorted(mp[v] for v in t))] for t in LT])
    return LT, perms


FT = {b: _flag_tables(*b) for b in BLOCKS}


def moments(h):
    """dict (s, m, sigma) -> dict (F1, F2) -> count, for the 3-graph h on [7] (definition in the docstring)."""
    E = [[[False] * 7 for _ in range(7)] for _ in range(7)]
    for t in itertools.combinations(range(7), 3):
        if (h >> cidx(t)) & 1:
            for a, b, c in itertools.permutations(t):
                E[a][b][c] = True
    out = {}
    for (s, m) in BLOCKS:
        LT, perms = FT[(s, m)]
        tl = triples(s)
        for th in itertools.permutations(range(7), s):
            sig = 0
            for t in tl:
                if E[th[t[0]]][th[t[1]]][th[t[2]]]:
                    sig |= 1 << cidx(t)
            rest = [v for v in range(7) if v not in th]
            cache = {}

            def flag(U):
                if U not in cache:
                    L = th + U
                    bits = [E[L[t[0]]][L[t[1]]][L[t[2]]] for t in LT]
                    best = None
                    for pm in perms:
                        c = 0
                        for i, bt in enumerate(bits):
                            if bt:
                                c |= 1 << pm[i]
                        best = c if best is None or c < best else best
                    cache[U] = best
                return cache[U]
            d = out.setdefault((s, m, sig), {})
            for U1 in itertools.combinations(rest, m - s):
                f1 = flag(U1)
                for U2 in itertools.combinations([v for v in rest if v not in U1], m - s):
                    key = (f1, flag(U2))
                    d[key] = d.get(key, 0) + 1
    return out


# ------------------------------------------------------------------------------------------------ exact linear algebra
def ldl(M):
    """Exact symmetric elimination.  Returns (psd, L, D, msg); L unit lower triangular (list of lists)."""
    n = len(M)
    A = [row[:] for row in M]
    L = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    D = [Fraction(0)] * n
    for i in range(n):
        p = A[i][i]
        if p < 0:
            return False, L, D, f"negative pivot at {i}"
        if p == 0:
            if any(A[i][j] != 0 for j in range(i + 1, n)):
                return False, L, D, f"zero pivot with nonzero row at {i}"
            continue
        D[i] = p
        rowi = A[i]
        nz = [k for k in range(i + 1, n) if rowi[k] != 0]
        for j in nz:                          # row j -= (A[j][i] / p) row i, on the whole trailing block
            f = A[j][i] / p
            L[j][i] = f
            Aj = A[j]
            for k in nz:
                Aj[k] -= f * rowi[k]
            Aj[i] = Fraction(0)
    return True, L, D, "ok"


def kernel_from_ldl(L, D):
    """Basis of ker(L D L^T): x with L^T x = e_i for each zero pivot i."""
    n = len(D)
    out = []
    for i in range(n):
        if D[i] != 0:
            continue
        x = [Fraction(0)] * n
        x[i] = Fraction(1)
        for r in range(i - 1, -1, -1):
            x[r] = -sum((L[k][r] * x[k] for k in range(r + 1, i + 1) if L[k][r] != 0 and x[k] != 0),
                        Fraction(0))
        out.append(x)
    return out


def matvec(M, x):
    return [sum((a * b for a, b in zip(row, x) if a != 0 and b != 0), Fraction(0)) for row in M]


# ------------------------------------------------------------------------------------------------ the checks
def analyse(name, mu_classes=None, mu=None, verbose=True):
    """Recompute phi and its moment data; if mu is given, run C2-C4.  Returns a result dict."""
    t0 = time.time()
    m, pat, it, prob, sig = CONSTRUCTIONS[name]
    phi = law(name)
    C = sorted(phi)
    P = [phi[c] for c in C]
    t = sum((p * F_value(c) for p, c in zip(P, C)), Fraction(0))
    res = {"construction": name, "problem": prob, "n_classes": len(C), "t": str(t),
           "C1_mass_1": sum(P) == 1, "C1_t_equals_1_minus_sigma": t == 1 - sig,
           "C1_all_free": all(free_of(c, prob) for c in C)}
    if verbose:
        print(f"# {name}: {len(C)} classes, t = E_phi f = {t}; checks {res}  ({time.time() - t0:.1f}s)", flush=True)
    mom = [moments(c) for c in C]
    keys = sorted(set().union(*[set(x) for x in mom]))
    pos = {c: i for i, c in enumerate(C)}
    if mu is not None:
        idx = [pos.get(int(c), -1) for c in mu_classes]
        res["C2_support_in_phi"] = all(i >= 0 for i in idx)
        res["C2_mu_nonneg"] = all(x >= 0 for x in mu)
        res["C2_mu_sum_1"] = sum(mu) == 1
        if not res["C2_support_in_phi"]:
            res["PASS"] = False
            return res
        muv = [Fraction(0)] * len(C)
        for i, x in zip(idx, mu):
            muv[i] += x
    data = []
    n_psd_fail = n_outside = n_kernel_bad = n_cond_bad = 0
    for key in keys:
        cnts = [mm.get(key, {}) for mm in mom]
        diag = {}
        for p, d in zip(P, cnts):
            for (a, b), v in d.items():
                if a == b:
                    diag[a] = diag.get(a, 0) + p * v
        O = sorted(f for f, v in diag.items() if v > 0)
        op = {f: i for i, f in enumerate(O)}
        n = len(O)
        Mphi = [[Fraction(0)] * n for _ in range(n)]
        outside = 0
        for p, d in zip(P, cnts):
            for (a, b), v in d.items():
                if a in op and b in op:
                    Mphi[op[a]][op[b]] += p * v
                else:
                    outside += 1
        psd, L, Dg, msg = ldl(Mphi)
        n_psd_fail += (not psd)
        n_outside += outside
        X = kernel_from_ldl(L, Dg) if psd else []
        for x in X:
            if any(v != 0 for v in matvec(Mphi, x)):
                n_kernel_bad += 1
        entry = {"key": key, "occurring": n, "rank": n - len(X), "kernel": len(X), "psd": psd, "msg": msg,
                 "outside_entries": outside, "O": O, "X": X, "cnts": cnts}
        if mu is not None and X:
            Mmu = [[Fraction(0)] * n for _ in range(n)]
            for wv, d in zip(muv, cnts):
                if wv == 0:
                    continue
                for (a, b), v in d.items():
                    Mmu[op[a]][op[b]] += wv * v
            Y = [matvec(Mmu, x) for x in X]
            bad = sum(1 for x in X for y in Y if sum((a * b for a, b in zip(x, y) if a and b), Fraction(0)) != 0)
            n_cond_bad += bad
            entry["cond_bad"] = bad
        data.append(entry)
    res.update({"n_keys": len(keys), "C3a_psd_fail": n_psd_fail, "C3b_entries_outside_O": n_outside,
                "C3c_kernel_check_fail": n_kernel_bad, "sum_kernel_dims": sum(e["kernel"] for e in data),
                "n_keys_with_kernel": sum(1 for e in data if e["kernel"])})
    if mu is not None:
        Emu = sum((w * F_value(c) for w, c in zip(muv, C) if w), Fraction(0))
        res["C3d_condition_fail"] = n_cond_bad
        res["E_mu_f"] = str(Emu)
        res["delta"] = str(t - Emu)
        res["delta_float"] = float(t - Emu)
        res["PASS"] = bool(res["C1_mass_1"] and res["C1_t_equals_1_minus_sigma"] and res["C1_all_free"]
                           and res["C2_support_in_phi"] and res["C2_mu_nonneg"] and res["C2_mu_sum_1"]
                           and n_psd_fail == 0 and n_outside == 0 and n_kernel_bad == 0 and n_cond_bad == 0)
    res["time_s"] = round(time.time() - t0, 1)
    res["_data"] = data
    res["_C"], res["_P"] = C, P
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("witness")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    W = json.load(open(a.witness))
    mu = [Fraction(x) for x in W["mu"]]
    res = analyse(W["construction"], W["classes"], mu)
    stated = Fraction(W["delta"])
    res["C4_delta_equals_file"] = Fraction(res.get("delta", "nan")) == stated if "delta" in res else False
    res["PASS"] = bool(res.get("PASS")) and res["C4_delta_equals_file"]
    public = {k: v for k, v in res.items() if not k.startswith("_")}
    print(json.dumps(public, indent=1), flush=True)
    if res["PASS"]:
        d = Fraction(res["delta"])
        t = Fraction(res["t"])
        print(f"VERDICT: PASS.  delta = {d} = {float(d):.6e} (support {sum(1 for x in mu if x)} classes).\n"
              f"Every plain N = 7 certificate (9 blocks, all labelled types, Q PSD) for the {res['problem']}-free "
              f"problem whose Q kills range M(phi) has b <= {t - d} = t - delta, t = {t};"
              + (" hence no certificate attains b = t." if d > 0 else " (delta = 0: no obstruction)."), flush=True)
    else:
        print("VERDICT: FAIL", flush=True)
    if a.json:
        with open(a.json, "w") as f:
            json.dump(public, f, indent=1)
    # exit 0: all checks pass and delta > 0 (obstruction); 2: all checks pass, delta = 0 (no obstruction); 1: FAIL
    sys.exit(1 if not res["PASS"] else (0 if Fraction(res["delta"]) > 0 else 2))


if __name__ == "__main__":
    main()
