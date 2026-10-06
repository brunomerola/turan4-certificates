"""Independent checker, CLAIM 1 (the dual point of Theorem n7opt): exact re-verification of the N = 7 dual point (own code; producer verifier NOT used).

Reads only the JSON dual point (colex masks, integer numerators, denominator) and, for the consistency check, the
certificate files lpcg_p5_conv_full.cert.{json,npz} (keys and integer factor rows).
Checks:
  (b) den = 2^110, sum num = den, num > 0, masks distinct and pairwise non-isomorphic, admissible (every 5-set has an
      edge); V = sum num e(H) / (35 den) exactly vs the claimed fraction;
  (c) for EVERY labelled type tau on [s] of EVERY block (s,m) in (1,4),(2,4),(3,4),(3,5),(4,5),(5,6), the integer
      matrix Y_tau = sum_H num_H * count_H,tau (count over ordered theta with H[theta] = tau exactly and ordered disjoint
      (U1,U2); proportional to the moment matrix), indexed by all admissible tau-flags (own enumeration): symmetric,
      zero-diagonal rows entirely zero, remaining principal submatrix PSD by exact elimination over Fractions;
  (d) oddness of all support graphs;
  (e) decimals;
  plus: weak-duality consistency b_cert <= sum_H y_H (d(H) - c_Q(H)) <= V with the certified Q; negative control
  (point mass on the certificate's argmin graph must FAIL).
Usage: python rv1_dual.py DUAL.json CERT_PREFIX OUT.json
"""
import json
import sys
import time
from decimal import Decimal, getcontext, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction

import numpy as np

import rv1_lib as L

BLOCKS = [(1, 4), (2, 4), (3, 4), (3, 5), (4, 5), (5, 6)]


def build_Y(support, blocks):
    """support: list of (own mask, integer weight). Returns {(s,m): {tau: {(f1,f2): int}}}."""
    Y = {b: {} for b in blocks}
    for h, w in support:
        for (s, m) in blocks:
            mc = L.moment_counts(h, s, m)
            for tau, d in mc.items():
                D = Y[(s, m)].setdefault(tau, {})
                for key, c in d.items():
                    D[key] = D.get(key, 0) + w * c
    return Y


def check_Y(Y, universes, log):
    allpass = True
    res = []
    for (s, m), bytype in Y.items():
        for tau in L.all_labelled_types(s):
            U = universes[(s, m, tau)]
            pos = {f: i for i, f in enumerate(U)}
            D = bytype.get(tau, {})
            occ = sorted({f for key in D for f in key})
            for f in occ:
                assert f in pos, ("flag outside the admissible universe", s, m, tau, f)
            # zero-diagonal rows must be zero
            diag = {f: D.get((f, f), 0) for f in occ}
            bad0 = [f for f in occ if diag[f] == 0 and any(D.get((f, g), 0) or D.get((g, f), 0) for g in occ)]
            live = [f for f in occ if diag[f] != 0]
            rows = [[D.get((a, b), 0) for b in live] for a in live]
            if bad0:
                ok, rank, why, worst = False, None, f"{len(bad0)} zero-diagonal rows with nonzero entries", None
            elif live:
                ok, rank, why, worst = L.psd_exact(rows)
            else:
                ok, rank, why, worst = True, 0, "zero matrix", None
            # float info: min eigenvalue of the diagonally normalised matrix
            lam = None
            if live:
                A = np.array([[float(Fraction(x)) for x in r] for r in rows])
                dd = 1 / np.sqrt(np.diag(A))
                lam = float(np.linalg.eigvalsh(A * dd[:, None] * dd[None, :]).min())
            ne = L.popcount(tau)
            r = {"s": s, "m": m, "tau_own": tau, "type_edges": ne, "universe": len(U), "occurring": len(occ),
                 "live": len(live), "PSD": ok, "rank": rank, "why": why, "lam_min_normalised_float": lam,
                 "min_pivot_over_diag": (float(worst) if worst is not None else None)}
            res.append(r)
            allpass &= ok
            if live or not ok:
                log(f"  ({s},{m}) tau={tau:>2} e(tau)={ne} universe={len(U):>4} occurring={len(occ):>3} "
                    f"live={len(live):>3} PSD={ok} rank={rank} {why} lam_min(norm,float)={lam} "
                    f"min pivot/diag={float(worst) if worst is not None else None}")
    return allpass, res


def load_cert(prefix):
    cert = json.load(open(prefix + ".cert.json"))
    Z = np.load(prefix + ".cert.npz")
    keys = []
    for k, key in enumerate(cert["keys"]):
        s, m = key["s"], key["m"]
        A = Z[f"A{k}"]
        n = len(key["flags"])
        assert A.shape[1] == n
        if A.shape[0] == 0:
            keys.append(None)
            continue
        # exact Gram matrix: max diagonal < 2^62 => no int64 overflow in any partial sum (Cauchy-Schwarz)
        Ao = [[int(v) for v in row] for row in A.tolist()]
        dg = [sum(Ao[r][i] ** 2 for r in range(len(Ao))) for i in range(n)]
        assert max(dg) < 2 ** 62
        Q = A.astype(np.int64).T @ A.astype(np.int64)
        assert all(int(Q[i, i]) == dg[i] for i in range(n))
        # own type code (lex on [s]) of the certificate's sigma, own canonical codes of its flags
        tau = L.code(L.colex_to_own(key["sigma"], s), tuple(range(s)), s)
        idx = {}
        for i, f in enumerate(key["flags"]):
            hf = L.colex_to_own(f, m)
            assert L.code(hf, tuple(range(s)), s) == tau, "flag root part != sigma"
            c = L.canon(hf, tuple(range(s)), tuple(range(s, m)), m)
            assert c not in idx, "two listed flags isomorphic"
            idx[c] = i
        keys.append({"s": s, "m": m, "tau": tau, "Q": Q, "idx": idx})
    return cert, keys


def eval_N(h, keys, Mfac, mc_cache=None):
    """Exact numerator 144 e M^2 - 5040 M^2 sum c_Q(H); value d - c = N / (5040 M^2)."""
    Z = 0
    bycache = {}
    for kk in keys:
        if kk is None:
            continue
        sm = (kk["s"], kk["m"])
        if mc_cache is not None and sm in mc_cache:
            mc = mc_cache[sm]
        else:
            mc = bycache.get(sm) or L.moment_counts(h, *sm)
            bycache[sm] = mc
        D = mc.get(kk["tau"], {})
        Q, idx = kk["Q"], kk["idx"]
        zb = 0
        for (a, b), c in D.items():
            if a in idx and b in idx:
                zb += c * int(Q[idx[a], idx[b]])
        Z += (5040 // L.conf(*sm)) * zb
    return 144 * L.popcount(h) * Mfac * Mfac - Z


def dec(x, nd, mode):
    getcontext().prec = 60
    q = Decimal(x.numerator) / Decimal(x.denominator)
    return q.quantize(Decimal(1).scaleb(-nd), rounding=mode)


def main():
    t0 = time.time()
    dual_fn, cert_prefix, out_fn = sys.argv[1:4]
    logf = open(out_fn.replace(".json", ".log"), "w")

    def log(*a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    D = json.load(open(dual_fn))
    den = int(D["den"])
    sup = [(int(x["mask"]), int(x["num"])) for x in D["support"]]
    out = {}
    log("# CLAIM 1 (b): basic facts")
    log("support size", len(sup), "den == 2^110:", den == 2 ** 110, "sum num == den:", sum(n for _, n in sup) == den,
        "all num > 0:", all(n > 0 for _, n in sup), "masks < 2^35:", all(0 <= m < 2 ** 35 for m, _ in sup))
    own = [(L.colex_to_own(m), n) for m, n in sup]
    assert all(L.own_to_colex(h) == m for (h, _), (m, _) in zip(own, sup))
    log("distinct masks:", len({h for h, _ in own}) == len(own))
    pm = L.perm_edge_maps()
    canon = [L.canonical_graph(h, pm) for h, _ in own]
    log("pairwise non-isomorphic (canonical forms over all 5040 permutations):", len(set(canon)) == len(canon))
    adm = [L.admissible5(h) for h, _ in own]
    odd = [L.is_odd(h) for h, _ in own]
    log("admissible (every 5-set has an edge):", all(adm), " odd (every 5-set spans 1,3,5 edges):", all(odd))
    prof = {}
    for h, n in own:
        for c in L.five_profile(h):
            prof[c] = prof.get(c, 0) + Fraction(n, den * 21)
    log("5-set edge-count profile under y:", {k: float(v) for k, v in sorted(prof.items())})
    V = Fraction(sum(n * L.popcount(h) for h, n in own), 35 * den)
    claimed = Fraction(D["value"])
    claimed2 = Fraction(1336237682928914994292138923, 7 * 2 ** 89)
    log("V =", V, "=", float(V))
    log("V == JSON value:", V == claimed, " V == 1336237682928914994292138923/(7*2^89):", V == claimed2)
    ecount = {}
    for h, n in own:
        ecount[L.popcount(h)] = ecount.get(L.popcount(h), 0) + Fraction(n, den)
    log("edge-count distribution of y:", {k: round(float(v), 6) for k, v in sorted(ecount.items())})
    out["b"] = {"support": len(sup), "sum_is_one": sum(n for _, n in sup) == den, "distinct_classes": len(set(canon)),
                "admissible": all(adm), "odd": all(odd), "V": str(V), "V_float": float(V),
                "V_equals_claim": V == claimed == claimed2}

    log("# CLAIM 1 (c): flag universes (own enumeration, all labelled types)")
    universes = {}
    for (s, m) in BLOCKS:
        sizes = {}
        for tau in L.all_labelled_types(s):
            universes[(s, m, tau)] = L.flag_universe(s, m, tau)
            sizes.setdefault(L.popcount(tau), set()).add(len(universes[(s, m, tau)]))
        log(f"  ({s},{m}) admissible flags per labelled type, by e(type):", {k: sorted(v) for k, v in sizes.items()})
    log("# CLAIM 1 (c): moment matrices of the dual point, exact PSD test")
    t1 = time.time()
    Y = build_Y(own, BLOCKS)
    log(f"  moment counts built in {time.time() - t1:.1f}s")
    ok, res = check_Y(Y, universes, log)
    log("ALL moment matrices PSD (all labelled types, all six blocks):", ok)
    out["c"] = {"all_psd": ok, "per_type": res}

    log("# negative control: point mass on the certificate argmin graph 27303369217 (colex)")
    hneg = L.colex_to_own(27303369217)
    Yn = build_Y([(hneg, 1)], BLOCKS)
    okn, resn = check_Y(Yn, universes, lambda *a: None)
    fails = [(r["s"], r["m"], r["tau_own"], r["why"]) for r in resn if not r["PSD"]]
    log("negative control all PSD:", okn, "(expected False); failing (block,type):", fails[:12], "...", len(fails))
    out["negative_control_fails"] = (not okn)

    log("# weak-duality consistency with the certified Q (lpcg_p5_conv_full)")
    cert, keys = load_cert(cert_prefix)
    Mfac = int(cert["M"])
    b = Fraction(cert["bound"])
    tot = Fraction(0)
    per = []
    for h, n in own:
        Nh = eval_N(h, keys, Mfac)
        val = Fraction(Nh, 5040 * Mfac * Mfac)
        per.append(val)
        tot += Fraction(n, den) * val
    log("sum_H y_H (d - c_Q) =", float(tot), " b_cert =", float(b), " V =", float(V))
    log("b_cert <= sum y (d - c) <= V :", b <= tot <= V, "; sum <Q, Y> = V - sum y (d - c) =", float(V - tot))
    log("min over support of (d - c) - b =", float(min(per) - b), " (>= 0 required:", min(per) >= b, ")")
    out["consistency"] = {"sum_y_slack": str(tot), "in_[b,V]": b <= tot <= V, "sum_QY": str(V - tot),
                          "min_support_value_minus_b": str(min(per) - b)}

    log("# CLAIM 1 (e): decimals")
    up = dec(V, 12, ROUND_CEILING)
    lo_pi = dec(1 - V, 12, ROUND_FLOOR)
    log("V exact decimal:", dec(V, 30, ROUND_FLOOR), "  ceil12(V) =", up, "  floor12(1-V) =", lo_pi)
    log("claimed 0.308401201195 >= V:", Fraction(str(up)) >= V and str(up) == "0.308401201195",
        " claimed 0.691598798805 <= 1-V:", str(lo_pi) == "0.691598798805")
    bdec = dec(b, 12, ROUND_FLOOR)
    log("two-sided width ceil12(V) - floor12(b) =", Fraction(str(up)) - Fraction(str(bdec)), float(Fraction(str(up)) - Fraction(str(bdec))))
    out["e"] = {"ceil12_V": str(up), "floor12_1mV": str(lo_pi)}
    out["time_s"] = time.time() - t0
    json.dump(out, open(out_fn, "w"), indent=1)
    log(f"done in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
