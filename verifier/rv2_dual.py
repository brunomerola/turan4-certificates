"""Independent checker: exact re-verification of the N = 7 dual points of certificates/n7_dual_points (own code; the
producer's verifier and engine are NOT used or imported).

For one dual JSON (colex masks, integer numerators, den):
  (1) den = 2^110, sum num = den, num >= 1, masks distinct, pairwise non-isomorphic (canonical forms over S_7),
      admissible for the target predicate (every p-set spans >= lam edges); the JSON's p / lam / blocks must equal
      the target's (fixed in TARGETS below from the claims, not from the JSON);
  (2) V = sum num e(H) / (35 den) exactly, = JSON value, = the claimed exact fraction;
  (3) for EVERY labelled type of EVERY listed block, Y_tau = sum_H num_H count_H,tau (own counting, rv1_lib),
      flags inside the own admissible universe, symmetric, zero-diagonal rows entirely zero, live part PSD:
        n <= 64      exact symmetric elimination over Fractions (rv1_lib.psd_exact)
        n <= NB      Bareiss / Sylvester (exact PD)
        all n        own float hint + exact residual (rv2_lib.pd_residual)
      (each criterion alone is a proof; all that apply are run and must agree);
  (4) information: odd graphs, p-set profile, mass concentration;
  (5) decimals ceil12(V), floor12(1 - V), and the derived numbers (V - b, widths);
  (6) b <= V; and, with the certificate (cert.json/npz in certificates/, hashes recorded), the weak-duality
      sandwich b <= sum_H y_H (d(H) - c_Q(H)) <= V with per-block <Q_k, Y_k> >= 0 (own evaluator).
  Controls: target p5no56_*: the (5,6) block is built as well and must FAIL (V < certified b with (5,6));
  perturbation controls of pd_residual on the largest matrix (A - 2 lam_min x x^T must fail, A - lam_min/2 x x^T pass).
Usage: python rv2_dual.py TAG [JSON]      (TAG in TARGETS, or neg:<path> with p,lam,blocks from the target of TAG)
"""
import hashlib
import json
import os
import sys
import time
from decimal import Decimal, getcontext, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction

import numpy as np

import rv1_lib as L
import rv2_lib as R

PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certificates")   # package certificates/
RES = os.path.join(PKG, "n7_dual_points")
ALL = [(1, 4), (2, 4), (3, 4), (3, 5), (4, 5), (5, 6)]
NO56 = ALL[:5]
NB = 160
CERT56 = PKG + "/K5_4/lpcg_p5_conv_full"      # t(5,4) certificate WITH (5,6), b = 35604499940047/115448720916480

TARGETS = {
    "p5no56_W": dict(p=5, lam=1, blocks=NO56, V=Fraction(390498160148386371476706381487, 35 * 2 ** 95),
                     ceil12="0.281644555065", pifloor="0.718355444935", b=None, cert=None),
    "p5no56_Wodd": dict(p=5, lam=1, blocks=NO56, V=Fraction(24406135009720244409729340413, 35 * 2 ** 91),
                        ceil12="0.281644555070", pifloor=None, b=None, cert=None),
    "p6_f3": dict(p=6, lam=1, blocks=ALL, V=Fraction(37867960095539974857512229, 7 * 2 ** 85),
                  ceil12="0.139837689411", pifloor="0.860162310589", b=Fraction(9224492822869, 15 * 2 ** 42),
                  picert="0.860172877518", width="1.056693e-5", cert=PKG + "/K6_4/c6l1"),
    "p7_f1": dict(p=7, lam=1, blocks=ALL, V=Fraction(6413961863433425283788748011, 5 * 2 ** 94),
                  ceil12="0.064764464150", pifloor="0.935235535850", b=Fraction(7476698908057, 115448720916480),
                  picert="0.935237923395", width="2.387545e-6", cert=PKG + "/K7_4/lpcg_p7_full"),
    "h44_nt3": dict(p=5, lam=2, blocks=ALL, V=Fraction(384021679403583125317285924719, 35 * 2 ** 94),
                    ceil12="0.553946861054", pifloor="0.446053138946", b=Fraction(48641756816487, 5 * 2 ** 44),
                    picert="0.447007248631", width="9.541097e-4", cert=PKG + "/K5_4minus/c5l2_f"),
    "k6m_W": dict(p=6, lam=2, blocks=ALL, V=Fraction(153033054540662097928475789999, 35 * 2 ** 94),
                  ceil12="0.220748423194", pifloor="0.779251576806", b=Fraction(202684381843379, 105 * 2 ** 43),
                  picert="0.780547176883", width="1.295600e-3", cert=PKG + "/K6_4minus/c6l2_f"),
}


def dec(x, nd, mode):
    getcontext().prec = 80
    q = Decimal(x.numerator) / Decimal(x.denominator)
    return q.quantize(Decimal(1).scaleb(-nd), rounding=mode)


def sha(fn):
    return hashlib.sha256(open(fn, "rb").read()).hexdigest()


def load_cert(prefix):
    """Certificate keys (as rv1_dual.load_cert): own type code, own canonical flag codes, exact Gram Q."""
    cert = json.load(open(prefix + ".cert.json"))
    Z = np.load(prefix + ".cert.npz")
    keys = []
    for k, key in enumerate(cert["keys"]):
        s, m = key["s"], key["m"]
        A = Z[f"A{k}"]
        n = len(key["flags"])
        assert A.shape[1] == n
        if A.shape[0] == 0:
            continue
        Ao = [[int(v) for v in row] for row in A.tolist()]
        dg = [sum(Ao[r][i] ** 2 for r in range(len(Ao))) for i in range(n)]
        assert max(dg) < 2 ** 62
        Q = A.astype(np.int64).T @ A.astype(np.int64)
        assert all(int(Q[i, i]) == dg[i] for i in range(n))
        tau = L.code(L.colex_to_own(key["sigma"], s), tuple(range(s)), s)
        idx = {}
        for i, f in enumerate(key["flags"]):
            hf = L.colex_to_own(f, m)
            assert L.code(hf, tuple(range(s)), s) == tau, "flag root part != sigma"
            c = L.canon(hf, tuple(range(s)), tuple(range(s, m)), m)
            assert c not in idx, "two listed flags isomorphic"
            idx[c] = i
        keys.append({"k": k, "s": s, "m": m, "tau": tau, "Q": Q, "idx": idx})
    return cert, keys


def key_values(mcs, keys):
    """Per key: zb = sum_{(a,b)} count * Q[a,b] (integer); c_Q(H) = sum_k zb_k / (conf_k M^2)."""
    out = []
    for kk in keys:
        D = mcs[(kk["s"], kk["m"])].get(kk["tau"], {})
        Q, idx = kk["Q"], kk["idx"]
        zb = 0
        for (a, b), c in D.items():
            if a in idx and b in idx:
                zb += c * int(Q[idx[a], idx[b]])
        out.append(zb)
    return out


def check_matrix(acc, universe, log, label, want_controls=False):
    flags, A = acc.matrix()
    n = len(flags)
    pos = set(universe)
    outside = [f for f in flags if f not in pos]
    r = {"label": label, "universe": len(universe), "occurring": n, "outside_universe": len(outside)}
    if outside:
        r.update(PSD=False, why="flag outside the admissible universe")
        log(f"  {label} univ={len(universe):>4} occ={n:>4} PSD=False {r['why']}")
        return r, None
    sym = all(A[i, j] == A[j, i] for i in range(n) for j in range(i + 1, n))
    if not sym:
        r.update(PSD=False, why="not symmetric")
        log(f"  {label} univ={len(universe):>4} occ={n:>4} PSD=False {r['why']}")
        return r, None
    zero = [i for i in range(n) if A[i, i] == 0]
    bad0 = [i for i in zero if any(A[i, j] != 0 for j in range(n))]
    live = [i for i in range(n) if A[i, i] != 0]
    r["live"] = len(live)
    if bad0:
        r.update(PSD=False, why=f"{len(bad0)} zero-diagonal rows with nonzero entries (exact proof of non-PSD)")
        log(f"  {label} univ={len(universe):>4} occ={n:>4} live={len(live):>4} PSD=False {r['why']}")
        return r, None
    if not live:
        r.update(PSD=True, why="zero matrix")
        log(f"  {label} univ={len(universe):>4} occ={n:>4} PSD=True zero matrix")
        return r, None
    Al = A[np.ix_(live, live)]
    proofs, fails = [], []
    if len(live) <= 64:
        ok, rank, why, worst = L.psd_exact([[Al[i, j] for j in range(len(live))] for i in range(len(live))])
        r["ldl_fraction"] = {"PSD": ok, "rank": rank, "why": why}
        (proofs if ok else fails).append("LDL-Fraction" + ("" if ok else ":" + why))
    if len(live) <= NB:
        ok, at = R.pd_bareiss(Al)
        r["bareiss"] = {"PD": ok, "first_nonpositive_minor": at}
        (proofs if ok else fails).append("Bareiss" + ("" if ok else f":minor{at}"))
    res = R.pd_residual(Al)
    r["residual"] = res
    (proofs if res["PASS"] else fails).append("residual" + ("" if res["PASS"] else ":" + res.get("why", "not dd")))
    r["PSD"] = bool(proofs)
    r["proofs"], r["fails"] = proofs, fails
    if not proofs:
        found, lam, q = R.neg_witness(Al)
        r["neg_witness"] = {"found": found, "lam_min": lam, "yAy": str(q)}
        r["why"] = "exact negative witness y^T A y < 0" if found else "UNDECIDED"
    log(f"  {label} univ={len(universe):>4} occ={n:>4} live={len(live):>4} PSD={r['PSD']} proofs={proofs} "
        f"fails={fails} lam_min={res.get('lam_min')} worstT={res.get('worst_ratio')}")
    return r, (Al if want_controls else None)


def main():
    t0 = time.time()
    tag = sys.argv[1]
    neg = tag.startswith("neg:")
    if neg:
        _, ttag, fn = tag.split(":", 2)
        tg = TARGETS[ttag]
        out_fn = f"results/{fn.split('/')[-1].replace('.json', '')}.rv2.json"
    else:
        tg = TARGETS[tag]
        fn = sys.argv[2] if len(sys.argv) > 2 else f"{RES}/dual_{tag}.json"
        out_fn = f"results/{tag}.rv2.json"
    logf = open(out_fn.replace(".json", ".log"), "w")

    def log(*a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    out = {"tag": tag, "json": fn, "json_sha256": sha(fn)}
    log("# rv2_dual", tag, fn, "sha256", out["json_sha256"])
    D = json.load(open(fn))
    p, lam, blocks = tg["p"], tg["lam"], tg["blocks"]
    jb = [tuple(b) for b in D["blocks"]]
    log("# target: every", p, "-set spans >=", lam, "edges; blocks", blocks)
    log("JSON p/lam/blocks equal the target:", (D["p"], D["lam"], jb) == (p, lam, blocks))
    out["json_matches_target"] = (D["p"], D["lam"], jb) == (p, lam, blocks)
    den = int(D["den"])
    sup = [(int(x["mask"]), int(x["num"])) for x in D["support"]]
    log("# (1) basic facts")
    f1 = {"support": len(sup), "den_is_2^110": den == 2 ** 110, "sum_is_den": sum(n for _, n in sup) == den,
          "all_num_ge_1": all(n >= 1 for _, n in sup), "masks_in_range": all(0 <= m < 2 ** 35 for m, _ in sup)}
    own = [(L.colex_to_own(m), n) for m, n in sup]
    assert all(L.own_to_colex(h) == m for (h, _), (m, _) in zip(own, sup))
    f1["distinct_masks"] = len({h for h, _ in own}) == len(own)
    cf = R.canonical_forms([h for h, _ in own])
    f1["pairwise_non_isomorphic"] = len(set(cf)) == len(cf)
    adm = [R.admissible(h, p, lam) for h, _ in own]
    f1["admissible"] = all(adm)
    log(f1)
    out["basic"] = f1
    if not all(adm):
        bad = [m for (m, _), a in zip(sup, adm) if not a]
        log("NOT ADMISSIBLE support graphs (colex):", bad[:10], "count", len(bad), "-> FAIL, stop")
        out["PASS"] = False
        json.dump(out, open(out_fn, "w"), indent=1)
        return
    # (4) information
    odd = [L.is_odd(h) for h, _ in own]
    prof = {}
    for h, n in own:
        for c in R.pset_counts(h, p):
            prof[c] = prof.get(c, 0) + n
    ncp = len(R.pset_counts(0, p))
    w = sorted((Fraction(n, den) for _, n in sup), reverse=True)
    acc, k = Fraction(0), 0
    while 1 - acc > Fraction(1, 10 ** 6):
        acc += w[k]
        k += 1
    oddmass = sum(Fraction(n, den) for (h, n), o in zip(own, odd) if o)
    info = {"odd_graphs": sum(odd), "odd_mass": float(oddmass), "carry_all_but_1e-6": k,
            "pset_profile": {c: float(Fraction(v, den * ncp)) for c, v in sorted(prof.items())}}
    log("# (4) info:", info)
    out["info"] = info
    # (2) V
    V = Fraction(sum(n * L.popcount(h) for h, n in own), 35 * den)
    f2 = {"V": str(V), "V_float": float(V), "V_eq_json": ("value" in D) and V == Fraction(D["value"]), "V_eq_claim": V == tg["V"]}
    log("# (2)", f2)
    out["V"] = f2

    # (3) moment matrices
    build_blocks = list(blocks)
    control56 = tag.startswith("p5no56") and not neg
    if control56:
        build_blocks = ALL
    certp = None if neg else tg.get("cert")
    if control56 and tag == "p5no56_W":
        certp = CERT56
    keys = []
    if certp:
        cert, keys = load_cert(certp)
        out["cert"] = {"prefix": certp, "json_sha256": sha(certp + ".cert.json"), "npz_sha256": sha(certp + ".cert.npz"),
                       "bound": cert["bound"], "M": cert["M"]}
        log("# certificate", out["cert"])
    log("# (3) universes")
    universes = {}
    for (s, m) in build_blocks:
        sizes = {}
        for tau in L.all_labelled_types(s):
            universes[(s, m, tau)] = R.flag_universe_pl(s, m, tau, p, lam)
            sizes.setdefault(L.popcount(tau), set()).add(len(universes[(s, m, tau)]))
        log(f"  ({s},{m}) admissible flags per labelled type, by e(type):", {k_: sorted(v) for k_, v in sizes.items()})
    log("# (3) moment counts")
    t1 = time.time()
    Y = {b: {} for b in build_blocks}
    kv_tot = [0] * len(keys)
    per_val = []
    Mf = int(cert["M"]) if certp else None
    for i, (h, wgt) in enumerate(own):
        mcs = {sm: L.moment_counts(h, *sm) for sm in build_blocks}
        for sm in build_blocks:
            for tau, d in mcs[sm].items():
                a = Y[sm].get(tau)
                if a is None:
                    a = Y[sm][tau] = R.Acc()
                a.add(d, wgt)
        if keys:
            zs = key_values(mcs, keys)
            cval = sum(Fraction(z, L.conf(kk["s"], kk["m"]) * Mf * Mf) for z, kk in zip(zs, keys))
            per_val.append(Fraction(L.popcount(h), 35) - cval)
            for j, z in enumerate(zs):
                kv_tot[j] += wgt * z
        if (i + 1) % 250 == 0:
            log(f"  {i + 1} graphs, {time.time() - t1:.0f}s")
    log(f"  moment counts done in {time.time() - t1:.0f}s")
    out["blocks"] = {}
    allpass = True
    biggest = None
    for sm in build_blocks:
        rs = []
        for tau in L.all_labelled_types(sm[0]):
            a = Y[sm].get(tau)
            label = f"({sm[0]},{sm[1]}) tau={tau:>2} e={L.popcount(tau)}"
            if a is None:
                rs.append({"label": label, "occurring": 0, "PSD": True, "why": "type never occurs (Y = 0)"})
                continue
            r, Al = check_matrix(a, universes[(sm[0], sm[1], tau)], log, label, want_controls=True)
            r["tau"], r["e_tau"] = tau, L.popcount(tau)
            rs.append(r)
            if sm in blocks:
                allpass &= r["PSD"]
                if Al is not None and r["PSD"] and (biggest is None or Al.shape[0] > biggest[0].shape[0]):
                    biggest = (Al, label)
            else:
                log("    ^ CONTROL block (not in the claim): expected to FAIL")
        out["blocks"][str(sm)] = rs
    out["all_listed_blocks_psd"] = allpass
    log("ALL listed blocks PSD (every labelled type):", allpass)
    # summary by e(type)
    for sm in build_blocks:
        summ = {}
        for r in out["blocks"][str(sm)]:
            if r.get("occurring"):
                e = r.get("e_tau")
                summ.setdefault(e, []).append((r["occurring"], r.get("live"), r["PSD"]))
        log(f"  summary {sm}: e(type) -> sorted set of (occurring, live, PSD):",
            {e: sorted(set(v)) for e, v in sorted(summ.items())})
    if control56:
        c = [r for r in out["blocks"][str((5, 6))] if r.get("occurring")]
        nf = sum(1 for r in c if not r["PSD"])
        log(f"CONTROL (5,6) for the no-(5,6) point: {nf} of {len(c)} occurring labelled types FAIL (expected >= 1)")
        out["control56_failing_types"] = nf

    # controls of pd_residual
    if biggest is not None and not neg:
        Al, label = biggest
        Am, Ah, lmin = R.perturbed_controls(Al)
        fm, lm, qm = R.neg_witness(Am)
        rm = R.pd_residual(Am)
        rh = R.pd_residual(Ah)
        log(f"# residual-method controls on {label} (n={Al.shape[0]}, lam_min={lmin:.3e}): "
            f"A-2lam xx^T exactly indefinite={fm}, residual PASS={rm['PASS']} (expected False); "
            f"A-lam/2 xx^T residual PASS={rh['PASS']} (expected True)")
        out["residual_controls"] = {"label": label, "minus_indefinite_exact": fm, "minus_residual_pass": rm["PASS"],
                                    "half_residual_pass": rh["PASS"]}

    # (6) certificate sandwich
    if keys:
        b = Fraction(cert["bound"])
        S = sum(Fraction(n, den) * v for (_, n), v in zip(own, per_val))
        QY = {f"key{kk['k']} ({kk['s']},{kk['m']}) tau={kk['tau']}":
              Fraction(z, den * L.conf(kk["s"], kk["m"]) * Mf * Mf) for z, kk in zip(kv_tot, keys)}
        mn = min(per_val)
        f6 = {"b": str(b), "sum_y_(d-c)": str(S), "sum_y_(d-c)_float": float(S), "b<=sum": b <= S, "sum<=V": S <= V,
              "V-sum=sum<Q,Y>": float(V - S), "min_support_(d-c)-b": float(mn - b), "min_support_ge_b": mn >= b,
              "per_key_QY_float": {k_: float(v) for k_, v in QY.items()},
              "all_key_QY_ge_0": all(v >= 0 for v in QY.values()),
              "neg_keys": [k_ for k_, v in QY.items() if v < 0]}
        assert sum(QY.values()) == V - S
        log("# (6) certificate sandwich:", json.dumps(f6))
        out["sandwich"] = f6
    # (5) decimals
    up = dec(V, 12, ROUND_CEILING)
    lo = dec(1 - V, 12, ROUND_FLOOR)
    f5 = {"V_30": str(dec(V, 30, ROUND_FLOOR)), "ceil12_V": str(up), "floor12_1-V": str(lo),
          "ceil12_claim_ok": str(up) == tg["ceil12"],
          "pifloor_claim_ok": (tg["pifloor"] is None) or str(lo) == tg["pifloor"]}
    if tg.get("b") is not None:
        b = tg["b"]
        fb = dec(b, 12, ROUND_FLOOR)
        f5.update({"b": str(b), "b_float": float(b), "b<=V": b <= V, "V-b": float(V - b), "floor12_b": str(fb),
                   "ceil12_1-b": str(dec(1 - b, 12, ROUND_CEILING)),
                   "picert_claim_ok": str(dec(1 - b, 12, ROUND_CEILING)) == tg["picert"],
                   "width_ceil12V-floor12b": float(Fraction(str(up)) - Fraction(str(fb))),
                   "width_claim": tg["width"],
                   "b_eq_cert_bound": (Fraction(out["cert"]["bound"]) == b) if "cert" in out else None})
    log("# (5) decimals:", json.dumps(f5))
    out["decimals"] = f5
    out["PASS"] = bool(allpass and all(f1.values()) and f2["V_eq_json"] and (neg or f2["V_eq_claim"])
                       and out["json_matches_target"])
    out["time_s"] = time.time() - t0
    json.dump(out, open(out_fn, "w"), indent=1, default=str)
    log(f"OVERALL {tag}: {'PASS' if out['PASS'] else 'FAIL'}  ({out['time_s']:.0f}s)")


if __name__ == "__main__":
    main()
