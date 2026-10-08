#!/usr/bin/env python3
"""Independent checker: exact verifier of an obstruction witness (Remark j4limit).  Reviewer code; imports only
jw_lib (reviewer) + stdlib/numpy/scipy.  Nothing from the producer.

Usage: python -B jw_verify.py WITNESS.json [--name NAME] [--nullity] [--fullnull] [--json OUT]
Checks (exact unless marked mod p / float):
  V1 phi: own law of the construction (own recursion over set partitions), own canonical classes, mass 1,
     E_phi f (own codegree formula, cross-checked by the 4-set formula), every class admissible (own tests).
  V2 witness classes (colex -> lex -> own canonical) are classes of phi, mu >= 0, sum 1, no duplicates.
  V3 for every (block, labelled type) occurring in phi: O = flags with positive diagonal of M(phi); every nonzero
     entry of every M(H), H in supp(phi), lies in O x O; M_OO(phi) PSD (fraction-free symmetric elimination);
     kernel of M_OO(phi) by exact RREF, dimension = |O| - rank (rank from the elimination), M_OO(phi) K = 0.
  V4 K^T M_OO(mu) K = 0 exactly for every such key.
  V5 delta = E_phi f - E_mu f exactly; compared with the file.
  --nullity: N1 rank mod p of the full linear condition system mu -> (K^T M(mu) K)_{a<=b} (93 unknowns) -> upper
     bound on the null space; N2 exact null space of the selected (mod-p pivot) rows (exact rows); N3 exact LP
     min E_mu f over {mu in N2, mu >= 0, sum mu = 1}: float HiGHS for the active set, then exact vertex + exact KKT
     multipliers; compared with the witness.  --fullnull additionally checks every basis vector of N2 against ALL
     conditions exactly (then N2 = the true null space).
Exit 0: all checks pass and delta > 0; 2: all pass, delta = 0; 1: a check fails.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import pickle
import sys
import time
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import jw_lib as L  # noqa: E402

LOCAL = os.path.join(os.getcwd(), "local")
PHI_REF = None
P_MOD = 67108859            # prime < 2^26 (int64 matmul of n <= 2048 terms stays exact)


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def get_phi(name):
    cache = os.path.join(LOCAL, f"phi_{name}.pkl")
    if os.path.exists(cache):
        with open(cache, "rb") as fh:
            return pickle.load(fh)
    law7, prob, desc = L.construction_law7(name)
    phi = L.classes_of(law7)
    res = dict(n_labelled=len(law7), phi=phi, prob=prob, desc=desc)
    os.makedirs(LOCAL, exist_ok=True)
    with open(cache, "wb") as fh:
        pickle.dump(res, fh)
    return res


def get_moments(name, classes, mutate=None, phi=None):
    cache = os.path.join(LOCAL, f"mom_{name}{'_' + mutate if mutate else ''}.pkl")
    if os.path.exists(cache):
        with open(cache, "rb") as fh:
            mom = pickle.load(fh)
        if set(mom) == set(classes):
            return mom
    if mutate == "edgeflip":     # MUTATION: moments of the most likely class computed from a graph with 1 triple toggled
        Hs = max(classes, key=lambda H: phi[H])
        mom = {H: L.moment_counts(H ^ 1 if H == Hs else H) for H in classes}
    else:
        mom = {H: L.moment_counts(H, mutate=mutate) for H in classes}
    with open(cache, "wb") as fh:
        pickle.dump(mom, fh)
    return mom


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("witness")
    ap.add_argument("--name", default=None)
    ap.add_argument("--nullity", action="store_true")
    ap.add_argument("--fullnull", action="store_true")
    ap.add_argument("--exactrows", action="store_true", help="count exactly nonzero condition rows (slow for J_4)")
    ap.add_argument("--json", default=None)
    ap.add_argument("--mutate", default=None, choices=["overlap", "typeblind", "edgeflip"],
                    help="SENSITIVITY TEST ONLY: deliberately wrong moment definition")
    a = ap.parse_args()
    t0 = time.time()
    W = json.load(open(a.witness))
    name = a.name or W["construction"]
    rep = {"witness": os.path.abspath(a.witness), "construction": name}
    fails = []

    def check(key, ok, info=""):
        rep[key] = bool(ok)
        print(f"  [{'ok' if ok else 'FAIL'}] {key} {info}")
        if not ok:
            fails.append(key)

    # ---------------- V1
    P = get_phi(name)
    phi = P["phi"]
    prob = P["prob"]
    classes = sorted(phi)
    print(f"# {name} ({P['desc']}), problem {prob}: {P['n_labelled']} labelled 7-vertex masks, {len(classes)} classes"
          f"  ({time.time() - t0:.1f}s)")
    rep.update(n_labelled=P["n_labelled"], n_classes=len(classes), problem=prob)
    check("V1_mass_1", sum(phi.values()) == 1)
    check("V1_all_positive", all(v > 0 for v in phi.values()))
    fv = {H: L.fval(H) for H in classes}
    t = sum(phi[H] * fv[H] for H in classes)
    rep["t_E_phi_f"] = str(t)
    check("V1_t_equals_file_target", t == Fraction(W["target"]), f"E_phi f = {t}, file target {W['target']}")
    check("V1_all_admissible", all(L.admissible(H, prob) for H in classes), f"(own {prob}-freeness test)")
    check("V1_problem_matches_file", prob == W["problem"], f"{prob} vs {W['problem']}")

    # ---------------- V2
    wcl = [L.colex_to_lex(int(c)) for c in W["classes"]]
    wcan = L.canon_many(wcl)
    mu_list = [Fraction(x) for x in W["mu"]]
    check("V2_lengths", len(wcan) == len(mu_list))
    check("V2_no_duplicate_classes", len(set(wcan)) == len(wcan))
    mu = {}
    for c, x in zip(wcan, mu_list):
        mu[c] = mu.get(c, Fraction(0)) + x
    supp_mu = [H for H in mu if mu[H] != 0]
    check("V2_supp_mu_in_supp_phi", all(H in phi for H in supp_mu),
          f"|supp mu| = {len(supp_mu)}, outside phi: {sum(H not in phi for H in supp_mu)}")
    check("V2_mu_nonneg", all(x >= 0 for x in mu.values()))
    check("V2_mu_sum_1", sum(mu.values()) == 1)
    check("V2_file_classes_cover_phi", set(wcan) >= set(classes), "(file lists every class of phi)")
    rep["supp_mu"] = len(supp_mu)
    # report the zero classes by the FILE's representative (colex), to compare with the README
    rep["zero_classes_file_masks"] = sorted(int(c) for c, cc, x in zip(W["classes"], wcan, mu_list) if x == 0)
    rep["mu_max"] = str(max(mu.values()))
    if fails and any(k.startswith("V2_supp") for k in fails):
        return finish(rep, fails, a, t0, None)

    # ---------------- moments
    mom = get_moments(name, classes, a.mutate, phi)
    rep["mutate"] = a.mutate
    print(f"# MUTATED moment definition: {a.mutate}" if a.mutate else "# moments: correct definition")
    print(f"# moments of {len(classes)} classes ({time.time() - t0:.1f}s)")
    tot_ok = all(sum(c for (s, m, ty), cnt in mom[H].items() if (s, m) == bl for c in cnt.values())
                 == L.n_configs(*bl) for H in classes for bl in L.BLOCKS)
    check("V3_config_totals", tot_ok or a.mutate is not None, "(each block: #configurations = 7!/(7-s)! C(7-s,m-s) C(7-m,m-s))")

    Lphi = L.lcm_den(phi.values())
    wphi = {H: int(phi[H] * Lphi) for H in classes}
    Lmu = L.lcm_den(mu.values())
    wmu = {H: int(mu.get(H, 0) * Lmu) for H in classes}
    keys = sorted({k for H in classes for k in mom[H]})
    rep["n_keys"] = len(keys)
    stats = dict(psd_fail=0, outside_O=0, kernel_fail=0, cond_fail=0, nker=0, sumdim=0, sumpairs=0, maxO=0,
                 rank_mismatch=0)
    KER = {}
    per_block = {}
    for key in keys:
        Mp = {}
        for H in classes:
            cnt = mom[H].get(key)
            if cnt:
                for fp, c in cnt.items():
                    Mp[fp] = Mp.get(fp, 0) + wphi[H] * c
        O = sorted({F for (F, G) in Mp if F == G and Mp[(F, G)] > 0})
        Oset = set(O)
        for H in classes:
            cnt = mom[H].get(key, {})
            if any((F not in Oset or G not in Oset) for (F, G), c in cnt.items() if c):
                stats["outside_O"] += 1
        if any((F not in Oset or G not in Oset) for (F, G), x in Mp.items() if x):
            # M(phi) has a zero diagonal entry with a nonzero row: not PSD
            stats["psd_fail"] += 1
            print("   PSD FAIL", key, "zero diagonal with nonzero off-diagonal entry in M(phi)")
            continue
        n = len(O)
        stats["maxO"] = max(stats["maxO"], n)
        ix = {F: i for i, F in enumerate(O)}
        A = [[0] * n for _ in range(n)]
        for (F, G), v in Mp.items():
            A[ix[F]][ix[G]] = v
        psd, rk, piv, zero, why = L.psd_bareiss(A)
        if not psd:
            stats["psd_fail"] += 1
            print("   PSD FAIL", key, why)
            continue
        bl = key[:2]
        pb = per_block.setdefault(bl, [0, 0, 0])
        pb[0] += 1
        pb[1] += n
        pb[2] += rk
        if rk == n:
            continue
        K, rk2 = L.rref_nullspace(A)
        if rk2 != rk or len(K) != n - rk:
            stats["rank_mismatch"] += 1
        if any(sum(A[i][j] * v[j] for j in range(n)) != 0 for v in K for i in range(n)):
            stats["kernel_fail"] += 1
        d = len(K)
        stats["nker"] += 1
        stats["sumdim"] += d
        stats["sumpairs"] += d * (d + 1) // 2
        KER[key] = (O, K)
        # V4: K^T M_OO(mu) K = 0
        Mm = {}
        for H in classes:
            if wmu[H]:
                for fp, c in mom[H].get(key, {}).items():
                    Mm[fp] = Mm.get(fp, 0) + wmu[H] * c
        MK = [[0] * d for _ in range(n)]
        for (F, G), v in Mm.items():
            i, j = ix[F], ix[G]
            row = MK[i]
            for b in range(d):
                kb = K[b][j]
                if kb:
                    row[b] += v * kb
        bad = 0
        for a_ in range(d):
            Ka = K[a_]
            for b in range(d):
                if sum(Ka[i] * MK[i][b] for i in range(n) if Ka[i]) != 0:
                    bad += 1
        if bad:
            stats["cond_fail"] += 1
            stats.setdefault("cond_fail_entries", 0)
            stats["cond_fail_entries"] += bad
    rep.update({f"V3_{k}": v for k, v in stats.items()})
    rep["per_block_keys_dimO_rank"] = {f"{s},{m}": v for (s, m), v in sorted(per_block.items())}
    print(f"# keys {len(keys)}, with kernel {stats['nker']}, sum kernel dims {stats['sumdim']}, "
          f"sum d(d+1)/2 {stats['sumpairs']}, max |O| {stats['maxO']}  ({time.time() - t0:.1f}s)")
    check("V3_psd", stats["psd_fail"] == 0, f"({stats['psd_fail']} failures)")
    check("V3_entries_in_OxO", stats["outside_O"] == 0)
    check("V3_kernel", stats["kernel_fail"] == 0 and stats["rank_mismatch"] == 0)
    check("V4_conditions", stats["cond_fail"] == 0, f"({stats['cond_fail']} keys fail)")

    # ---------------- V5
    Emu = sum(mu.get(H, 0) * fv[H] for H in classes)
    delta = t - Emu
    rep["E_mu_f"] = str(Emu)
    rep["delta"] = str(delta)
    rep["delta_float"] = float(delta)
    check("V5_delta_equals_file", delta == Fraction(W["delta"]), f"delta = {delta} = {float(delta):.10e}")
    check("V5_delta_nonneg", delta >= 0)
    if delta > 0:
        dd = delta.denominator
        fac, x, pr = [], dd, 2
        while pr * pr <= x:
            while x % pr == 0:
                fac.append(pr)
                x //= pr
            pr += 1
        if x > 1:
            fac.append(x)
        rep["delta_den_factorisation"] = fac
        rep["delta_num_prime"] = is_prime(delta.numerator)

    if a.nullity:
        global PHI_REF
        PHI_REF = phi
        nullity_analysis(rep, check, classes, mom, KER, fv, mu, t, a.fullnull)
    if a.exactrows:
        nz = 0
        tot = 0
        for key, (O, K) in KER.items():
            ix = {F: i for i, F in enumerate(O)}
            for ia in range(len(K)):
                for ib in range(ia, len(K)):
                    tot += 1
                    Ka, Kb = K[ia], K[ib]
                    if any(sum(c * Ka[ix[F]] * Kb[ix[G2]] for (F, G2), c in mom[H].get(key, {}).items()) != 0
                           for H in classes):
                        nz += 1
        rep["exact_rows_total"] = tot
        rep["exact_rows_nonzero"] = nz
        print(f"# exact condition rows: {tot}, nonzero {nz}")
    return finish(rep, fails, a, t0, delta)


def nullity_analysis(rep, check, classes, mom, KER, fv, mu, t, fullnull):
    t0 = time.time()
    p = P_MOD
    nc = len(classes)
    rows_p = []
    rowid = []
    for key, (O, K) in KER.items():
        ix = {F: i for i, F in enumerate(O)}
        n, d = len(O), len(K)
        Kp = np.array([[x % p for x in v] for v in K], dtype=np.int64).T        # (n, d)
        G = np.zeros((nc, d, d), dtype=np.int64)
        for h, H in enumerate(classes):
            cnt = mom[H].get(key)
            if not cnt:
                continue
            C = np.zeros((n, n), dtype=np.int64)
            for (F, G2), c in cnt.items():
                C[ix[F], ix[G2]] = c % p
            CK = (C @ Kp) % p
            G[h] = (Kp.T @ CK) % p
        iu = np.triu_indices(d)
        R = G[:, iu[0], iu[1]].T          # (pairs, nc)
        for r_, (aa, bb) in zip(R, zip(*iu)):
            rows_p.append(r_)
            rowid.append((key, int(aa), int(bb)))
    rows_p = np.array(rows_p, dtype=np.int64)
    rk_p, piv = L.rank_mod_p(rows_p, p)
    rep["N1_condition_rows"] = int(rows_p.shape[0])
    rep["N1_rank_mod_p"] = rk_p
    rep["N1_nullity_upper_bound"] = nc - rk_p
    rep["N1_rows_nonzero_mod_p"] = int((rows_p != 0).any(axis=1).sum())
    print(f"# N1: {rows_p.shape[0]} condition rows, rank mod p {rk_p}, nullity <= {nc - rk_p} ({time.time() - t0:.1f}s)")
    # N2: exact rows of the pivot set
    ex_rows = []
    for ri in piv:
        key, aa, bb = rowid[ri]
        O, K = KER[key]
        ix = {F: i for i, F in enumerate(O)}
        Ka, Kb = K[aa], K[bb]
        row = []
        for H in classes:
            cnt = mom[H].get(key, {})
            row.append(sum(c * Ka[ix[F]] * Kb[ix[G2]] for (F, G2), c in cnt.items()))
        ex_rows.append(row)
    NB, rk_ex = L.rref_nullspace(ex_rows) if ex_rows else ([[int(i == j) for i in range(nc)] for j in range(nc)], 0)
    rep["N2_exact_rank_selected"] = rk_ex
    rep["N2_exact_nullity_selected"] = len(NB)
    check("N2_rank_consistent", rk_ex == rk_p, f"(exact rank of the {len(piv)} selected rows = {rk_ex})")
    Lphi = L.lcm_den(PHI_REF.values())
    check("N2_phi_in_null_space", all(sum(r_[h] * int(PHI_REF[H] * Lphi) for h, H in enumerate(classes)) == 0
                                      for r_ in ex_rows))
    # phi in null space?  (it must be: M(phi) K = 0)
    if fullnull:
        allok = True
        for v in NB:
            vv = {H: v[i] for i, H in enumerate(classes)}
            for key, (O, K) in KER.items():
                ix = {F: i for i, F in enumerate(O)}
                M = {}
                for H in classes:
                    if vv[H]:
                        for fp, c in mom[H].get(key, {}).items():
                            M[fp] = M.get(fp, 0) + vv[H] * c
                for Ka in K:
                    for Kb in K:
                        if sum(c * Ka[ix[F]] * Kb[ix[G2]] for (F, G2), c in M.items()) != 0:
                            allok = False
        check("N2_full_null_space_exact", allok, f"(all {len(NB)} basis vectors satisfy every condition exactly)")
    # N3 exact LP
    from scipy.optimize import linprog
    d0 = len(NB)
    Nm = [[NB[j][i] for j in range(d0)] for i in range(nc)]          # nc x d0 ints
    scale = [max(abs(Nm[i][j]) for i in range(nc)) or 1 for j in range(d0)]
    Nf = np.array([[Nm[i][j] / scale[j] for j in range(d0)] for i in range(nc)], dtype=float)
    fvec = np.array([float(fv[H]) for H in classes])
    res = linprog(Nf.T @ fvec, A_ub=-Nf, b_ub=np.zeros(nc), A_eq=(Nf.sum(axis=0))[None, :], b_eq=[1.0],
                  bounds=[(None, None)] * d0, method="highs")
    rep["N3_float_status"] = res.status
    if res.status != 0:
        check("N3_lp_float", False, res.message)
        return
    mu_f = Nf @ res.x
    order = np.argsort(mu_f)
    Zall = [int(i) for i in range(nc) if mu_f[i] <= 1e-12]
    rep["N3_float_zero_count"] = len(Zall)
    rep["N3_float_smallest_positive"] = float(min(x for x in mu_f if x > 1e-12))
    u = [sum(Nm[i][j] for i in range(nc)) for j in range(d0)]
    g = [sum(Nm[i][j] * fv[H] for i, H in enumerate(classes)) for j in range(d0)]
    import itertools as _it
    found = None
    nbases = 0
    for Z in _it.combinations(Zall, d0 - 1):
        Aeq = [[Fraction(Nm[i][j]) for j in range(d0)] + [Fraction(0)] for i in Z] + \
              [[Fraction(x) for x in u] + [Fraction(-1)]]
        ker, rk = L.rref_nullspace(Aeq)
        if not (len(ker) == 1 and ker[0][-1] != 0):
            continue                      # dependent active set
        nbases += 1
        v = ker[0]
        c = [Fraction(v[j], v[-1]) for j in range(d0)]
        mu_star = [sum(Nm[i][j] * c[j] for j in range(d0)) for i in range(nc)]
        if not (all(x >= 0 for x in mu_star) and sum(mu_star) == 1):
            continue
        val = sum(mu_star[i] * fv[H] for i, H in enumerate(classes))
        cols = [u] + [Nm[i] for i in Z]
        Mk = [[Fraction(cols[q][j]) for q in range(len(cols))] + [Fraction(-g[j])] for j in range(d0)]
        kk, _ = L.rref_nullspace(Mk)
        if len(kk) == 1 and kk[0][-1] != 0:
            sol = [Fraction(kk[0][q], kk[0][-1]) for q in range(len(cols))]
            lam, y = sol[0], sol[1:]
            if all(x >= 0 for x in y) and lam == val:
                found = (Z, mu_star, val, lam, y)
                break
    rep["N3_bases_tried"] = nbases
    ok_kkt = found is not None
    check("N3_vertex_feasible", ok_kkt or nbases > 0)
    if not ok_kkt:
        check("N3_vertex_optimal_exact_KKT", False)
        return
    Z, mu_star, val, lam, y = found
    rep["N3_kkt_lambda"] = str(lam)
    rep["N3_kkt_min_y"] = float(min(y)) if y else None
    check("N3_vertex_optimal_exact_KKT", ok_kkt)
    Emu = sum(mu.get(H, 0) * fv[H] for H in classes)
    rep["N3_lp_optimum"] = str(val)
    rep["N3_lp_delta"] = str(t - val)
    check("N3_lp_optimum_equals_witness", val == Emu, f"(LP delta = {t - val} = {float(t - val):.10e})")
    rep["N3_lp_vertex_equals_witness_mu (info)"] = all(mu_star[i] == mu.get(H, 0) for i, H in enumerate(classes))
    rep["N3_time_s"] = round(time.time() - t0, 1)


def finish(rep, fails, a, t0, delta):
    rep["fails"] = fails
    rep["PASS"] = not fails
    rep["time_s"] = round(time.time() - t0, 1)
    code = 1 if fails else (0 if (delta is not None and delta > 0) else 2)
    rep["exit"] = code
    print(f"VERDICT: {'PASS' if not fails else 'FAIL'} exit {code}; delta = {rep.get('delta')}; fails = {fails}"
          f"  ({rep['time_s']}s)")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(rep, fh, indent=1, default=str)
    return code


if __name__ == "__main__":
    sys.exit(main())
