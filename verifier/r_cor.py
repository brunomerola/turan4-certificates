"""Independent checker: exact checks of the hand-proved parts (own code, st_lib only).
 A. exact finite formula for f(Gamma) (Corollary 3.2) against brute force on own random instances;
 B. the O(1/n) numerics of Corollary 3.3 (phi = 1/2 and Hadamard-type phi = 1/2 - 1/(n-2), balanced parts);
 C. single flips in Giraud systems: the n = 14 doubly-degenerate example (Delta Lambda = -(n-4)/5), a large-n
    near-optimal Sylvester-Hadamard instance, the (n-4)/10 bound for non-degenerate systems and the exact
    characterisation of admissible cross-edge removals;
 D. Cauchy-Schwarz numbers;
 E. frame counts of 10-vertex Giraud systems (rho(Gamma) = [(a)_5 b + (b)_5 a]/(10)_6, minimum 1/126);
 F. sanity test of Claim R1 on perturbed Giraud systems (n = 13), not part of any proof.
Output: results_cor.json
"""
from __future__ import annotations

import json
import os
import random
import sys
import time
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, sqrt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.getcwd()   # outputs (results_*.json, local/) go to the working directory
import st_lib as L  # noqa: E402


def even_minors(A, B, M):
    return sum(1 for a, a2 in combinations(A, 2) for b, b2 in combinations(B, 2)
               if (M[(a, b)] + M[(a, b2)] + M[(a2, b)] + M[(a2, b2)]) % 2 == 0)


def f_formula(n, a, b, ev):
    """Producer's Corollary 3.2 formula (as stated in REPORT.txt), re-implemented."""
    t = Fraction(comb(a, 5) + comb(b, 5)) + Fraction(2, 5) * (comb(a, 4) * b + a * comb(b, 4))
    X = comb(a, 2) * comb(b, 2)
    if X:
        phi = Fraction(ev, X)
        t += Fraction(2, 15) * (n - 4) * X + Fraction(n - 4, 12) * X * (3 * phi - 1)
    return t / comb(n, 5)


def f_formula_mine(n, a, b, ev):
    """Reviewer's own derivation: (3,2)/(2,3) 5-sets contribute (n-4) X/20 + (n-4) Ev/4."""
    t = Fraction(comb(a, 5) + comb(b, 5)) + Fraction(2, 5) * (comb(a, 4) * b + a * comb(b, 4))
    X = comb(a, 2) * comb(b, 2)
    t += Fraction((n - 4) * X, 20) + Fraction((n - 4) * ev, 4)
    return t / comb(n, 5)


def rand_M(A, B, rng):
    return {(x, y): rng.randrange(2) for x in A for y in B}


def partA(rng):
    ok = True
    cnt = 0
    for n in range(8, 17):
        for a in sorted({0, 1, 2, 3, n // 2 - 1, n // 2, n // 2 + 1, n - 3, n - 1, n}):
            for trial in range(2):
                A = list(range(a))
                B = list(range(a, n))
                M = rand_M(A, B, rng)
                if trial == 1 and a >= 2 and n - a >= 2:      # structured: two equal rows and two equal columns
                    for y in B:
                        M[(1, y)] = M[(0, y)]
                    for x in A:
                        M[(x, B[1])] = M[(x, B[0])]
                g = L.giraud(n, A, M)
                ev = even_minors(A, B, M)
                fb = L.f5(g, n)
                ok &= fb == f_formula(n, a, n - a, ev) == f_formula_mine(n, a, n - a, ev)
                if n <= 10:
                    ok &= fb == L.f_def(g, n)
                cnt += 1
    return {"instances": cnt, "formula_equals_bruteforce": ok}


def f_bal(n, phi):
    a = b = n // 2
    X = comb(a, 2) * comb(b, 2)
    return f_formula(n, a, b, phi * X)


def partB():
    out = {}
    for n in (10 ** 3, 10 ** 4, 10 ** 6):
        out[str(n)] = {"n(f-33/64)_phi=1/2": float(n * (f_bal(n, Fraction(1, 2)) - Fraction(33, 64))),
                       "n(f-33/64)_phi=1/2-1/(n-2)": float(n * (f_bal(n, Fraction(1, 2) - Fraction(1, n - 2))
                                                                - Fraction(33, 64)))}
    # limit of the asymptotic expansion f(s, phi) = 1 - 3 s + s^2/2 + 15/2 phi s^2 at (1/4, 1/2)
    s, ph = Fraction(1, 4), Fraction(1, 2)
    out["limit_formula_value"] = str(1 - 3 * s + s * s / 2 + Fraction(15, 2) * ph * s * s)
    # (1/4 - s)(3 - 17/4 (s + 1/4)) == 31/64 - 3 s + 17/4 s^2 identity, and the 7/8 bound on [0, 1/4]
    ok = True
    for k in range(0, 101):
        s = Fraction(k, 400)
        ok &= (Fraction(1, 4) - s) * (3 - Fraction(17, 4) * (s + Fraction(1, 4))) == Fraction(31, 64) - 3 * s + \
            Fraction(17, 4) * s * s
        ok &= 3 - Fraction(17, 4) * (s + Fraction(1, 4)) >= Fraction(7, 8)
        for ph in (Fraction(0), Fraction(1, 2), Fraction(1)):
            lhs = 1 - 3 * s + s * s / 2 + Fraction(15, 2) * ph * s * s - Fraction(33, 64)
            rhs = (Fraction(1, 4) - s) * (3 - Fraction(17, 4) * (s + Fraction(1, 4))) + Fraction(15, 2) * s * s * (
                ph - Fraction(1, 2))
            ok &= lhs == rhs
    out["expansion_identities_ok"] = ok
    return out


def five_map(g, n):
    return {S: c for S, c in zip(combinations(range(n), 5), L.five_counts(g, n))}


def flip_effect(g, n, T, fm):
    """(admissible after the flip, Delta Lambda) for flipping the 4-set T."""
    e = L.bit(g, T)
    d = Fraction(0)
    adm = True
    for z in range(n):
        if z in T:
            continue
        S = tuple(sorted(T + (z,)))
        c0 = fm[S]
        c1 = c0 - 1 if e else c0 + 1
        if c1 == 0:
            adm = False
        d += L.G5[c1] - L.G5[c0]
    return adm, d


def row_deg(M, x, y, B):
    return len({(M[(x, c)] + M[(y, c)]) % 2 for c in B}) == 1


def col_deg(M, x, y, A):
    return len({(M[(r, x)] + M[(r, y)]) % 2 for r in A}) == 1


def partC(rng):
    out = {}
    # C1: the n = 14 example (several random completions)
    ex = []
    for seed in range(5):
        r = random.Random(1000 + seed)
        n, a = 14, 7
        A, B = list(range(a)), list(range(a, n))
        M = rand_M(A, B, r)
        for y in B:
            M[(1, y)] = M[(0, y)]
        for x in A:
            M[(x, a + 1)] = M[(x, a)]
        g = L.giraud(n, A, M)
        T = (0, 1, a, a + 1)
        assert L.bit(g, T) == 1
        h = g ^ (1 << L.cidx(T))
        dl = L.Lambda(h, n) - L.Lambda(g, n)
        ex.append({"seed": seed, "admissible_after": L.admissible(h, n), "delta_Lambda": str(dl),
                   "equals_-(n-4)/5": dl == -Fraction(n - 4, 5), "f_before": str(L.f5(g, n)),
                   "f_after": str(L.f5(h, n)), "odd_after": L.is_odd(h, n),
                   "giraud_after_forced": bool(L.is_giraud_forced(h, n)) if seed == 0 else None,
                   "five_counts_through_T_before": sorted({five_map(g, n)[tuple(sorted(T + (z,)))]
                                                           for z in range(n) if z not in T})})
    out["n14_examples"] = ex
    # C2: large n: Sylvester-Hadamard 16 x 16 (0/1 form) with row 1 := row 0 and column 1 := column 0, n = 32
    k = 16
    H01 = [[bin(i & j).count("1") % 2 for j in range(k)] for i in range(k)]
    n = 2 * k
    A, B = list(range(k)), list(range(k, n))
    M = {(i, k + j): H01[i][j] for i in range(k) for j in range(k)}
    ev_h = even_minors(A, B, M)
    f_h = f_formula(n, k, k, ev_h)
    for j in range(k):
        M[(1, k + j)] = M[(0, k + j)]
    for i in range(k):
        M[(i, k + 1)] = M[(i, k)]
    ev = even_minors(A, B, M)
    fM = f_formula(n, k, k, ev)
    # effect of removing T = {0,1,k,k+1}: the n-4 five-sets through T, counted directly from the definition
    T = (0, 1, k, k + 1)

    def e_of(S):
        Aset = set(A)
        cnt = 0
        for q in combinations(S, 4):
            ins = [v for v in q if v in Aset]
            if len(ins) in (0, 4):
                cnt += 1
            elif len(ins) == 2:
                outs = [v for v in q if v not in Aset]
                if sum(M[(x, y)] for x in ins for y in outs) % 2 == 0:
                    cnt += 1
        return cnt
    cs = [e_of(tuple(sorted(T + (z,)))) for z in range(n) if z not in T]
    out["n32_hadamard_degenerate"] = {"phi_hadamard": str(Fraction(ev_h, comb(k, 2) ** 2)),
                                      "f_hadamard-33/64": float(f_h - Fraction(33, 64)),
                                      "phi_after_making_degenerate": str(Fraction(ev, comb(k, 2) ** 2)),
                                      "f-33/64": float(fM - Fraction(33, 64)),
                                      "counts_through_T": sorted(set(cs)),
                                      "delta_Lambda_removal": str(sum(L.G5[c - 1] - L.G5[c] for c in cs)),
                                      "equals_-(n-4)/5": sum(L.G5[c - 1] - L.G5[c] for c in cs) == -Fraction(n - 4, 5),
                                      "admissible_after": all(c - 1 >= 1 for c in cs)}
    # C3: all single flips of random balanced systems (n = 12, 14), plus forced-degenerate ones
    res = []
    allok = True
    for n in (12, 14):
        for trial in range(6):
            a = n // 2
            A, B = list(range(a)), list(range(a, n))
            M = rand_M(A, B, rng)
            if trial >= 3:                      # force some degeneracy
                if trial in (3, 5):
                    for y in B:
                        M[(1, y)] = M[(0, y)] ^ (1 if trial == 5 else 0)
                if trial in (4, 5):
                    for x in A:
                        M[(x, a + 1)] = M[(x, a)]
            g = L.giraud(n, A, M)
            fm = five_map(g, n)
            dd_quads = set()
            for x, y in combinations(A, 2):
                if not row_deg(M, x, y, B):
                    continue
                for u, w in combinations(B, 2):
                    if col_deg(M, u, w, A):
                        dd_quads.add((x, y, u, w))
            mins = {}
            char_ok = True
            bound_ok = True
            for T in L.subsets(n, 4):
                adm, d = flip_effect(g, n, T, fm)
                ka = sum(1 for v in T if v < a)
                e = L.bit(g, T)
                kind = ("remove" if e else "add") + f"_{ka}{4 - ka}"
                if e and ka == 2:
                    # removal of a cross edge: admissible iff doubly degenerate
                    char_ok &= adm == (T in dd_quads)
                    if adm:
                        char_ok &= d == -Fraction(n - 4, 5)
                if adm:
                    mins[kind] = min(mins.get(kind, d), d)
                    if not (e and ka == 2):
                        bound_ok &= d >= Fraction(n - 4, 10)
                if e and ka in (0, 4):
                    char_ok &= not adm          # internal edge removal never admissible (other part nonempty)
            allok &= char_ok and bound_ok
            res.append({"n": n, "trial": trial, "doubly_degenerate_quads": len(dd_quads),
                        "min_delta_by_kind": {k_: str(v) for k_, v in sorted(mins.items())},
                        "cross_removal_characterisation_ok": char_ok,
                        "all_other_admissible_flips_ge_(n-4)/10": bound_ok})
    out["single_flips"] = res
    out["single_flips_all_ok"] = allok
    return out


def partD():
    t = 1 - sqrt(31) / 8
    return {"1-sqrt(31)/8": t, "floor_1e-5": int(t * 1e5) / 1e5, "2d-d^2_at_5/16": str(2 * Fraction(5, 16) -
                                                                                     Fraction(5, 16) ** 2),
            "135/256_float": 135 / 256, "33/64": 33 / 64}


def frames_count(g, n):
    cnt = 0
    for t in permutations(range(n), 6):
        f, a1, a2, a3, a4, b0 = t
        if L.bit(g, (f, a1, a2, b0)):
            continue
        S = (f, a1, a2, a3, a4)
        if all(L.bit(g, q) for q in combinations(S, 4)):
            cnt += 1
    return cnt


def partE(rng):
    out = []
    ok = True
    for a in range(5, 11):
        for trial in range(2):
            A, B = list(range(a)), list(range(a, 10))
            perm = list(range(10))
            rng.shuffle(perm)
            M = rand_M(A, B, rng)
            g = L.giraud(10, A, M)
            g = L.relabel(g, 10, perm)
            c = frames_count(g, 10)
            pred = (comb(a, 5) * 120 * (10 - a) + comb(10 - a, 5) * 120 * a)
            ok &= c == pred
            out.append({"a": a, "frames": c, "predicted": pred})
    return {"instances": out, "all_equal": ok, "min_nonzero_over_a": min(x["predicted"] for x in out if x["a"] < 10),
            "(10)_6": 151200}


def partF(rng):
    """Claim R1 sanity: n = 13, perturbed balanced Giraud system; random frames of H; every 4-set T disjoint from
    F with H[F u T] Giraud must satisfy H(T) = Gamma_F(T)."""
    n, a = 13, 6
    A, B = list(range(a)), list(range(a, n))
    M = rand_M(A, B, rng)
    g = L.giraud(n, A, M)
    h = g
    for _ in range(6):
        h ^= 1 << rng.randrange(comb(n, 4))
    frames = []
    while len(frames) < 25:
        t = rng.sample(range(n), 6)
        f, a1, a2, a3, a4, b0 = t
        if L.bit(h, (f, a1, a2, b0)) or not all(L.bit(h, q) for q in combinations((f, a1, a2, a3, a4), 4)):
            continue
        frames.append(t)
    tested = viol = giraud_windows = 0
    for F in frames:
        f, a1, a2, a3, a4, b0 = F
        AF = {f, a1, a2, a3, a4} | {x for x in range(n) if x not in F and L.bit(h, (f, a1, a2, x))}
        BF = [x for x in range(n) if x not in AF]
        MF = {}
        for x in AF:
            for y in BF:
                MF[(x, y)] = 0 if (x == f or y == b0) else 1 - L.bit(h, (f, x, b0, y))
        gF = L.giraud(n, sorted(AF), MF)
        rest = [x for x in range(n) if x not in F]
        for T in combinations(rest, 4):
            S = sorted(set(F) | set(T))
            sub = L.induced(h, n, S)
            tested += 1
            if L.is_giraud_struct(sub, 10):
                giraud_windows += 1
                if L.bit(h, T) != L.bit(gF, T):
                    viol += 1
    return {"n": n, "perturbations": 6, "frames": len(frames), "pairs_tested": tested,
            "giraud_windows": giraud_windows, "violations": viol}


def main():
    t0 = time.time()
    rng = random.Random(53_20261006)
    out = {}
    out["A_formula"] = partA(rng)
    print(out["A_formula"], flush=True)
    out["B_numerics"] = partB()
    print(out["B_numerics"], flush=True)
    out["C_flips"] = partC(rng)
    print(json.dumps(out["C_flips"])[:3000], flush=True)
    out["D_cauchy_schwarz"] = partD()
    print(out["D_cauchy_schwarz"], flush=True)
    out["E_frames"] = partE(rng)
    print(out["E_frames"], flush=True)
    out["F_claimR1_sanity"] = partF(rng)
    print(out["F_claimR1_sanity"], flush=True)
    out["time_s"] = round(time.time() - t0, 1)
    json.dump(out, open(os.path.join(OUT, "results_cor.json"), "w"), indent=1)
    print("DONE", out["time_s"])


if __name__ == "__main__":
    main()
