"""N = 6 control for t(p, 4), p = 5, 6, 7.  Producer code; imports clarabel only (no highspy,
no ortools).

For each p:
 1. plain flag-algebra SDP at N = 6 with ALL blocks (s, m), 2m - s <= 6 (engine ../flags.py, ../sdp.py);
 2. lower bound: Flagmatic-style sharp rounding to the target 1/C(p-1, 3) (de Caen), exact rational certificate
    written to results/cert_turan_r4_p{p}_N6_all_sharp.json (checked independently by ../verify_cert.py);
 3. UPPER bound on the N = 6 plain flag-algebra optimum (the point of the control): a rational "pseudo-density"
    vector y >= 0 on the admissible 6-vertex graphs with sum y = 1 and every moment matrix
        Y_sigma = sum_H y_H M_sigma(H)
    positive semidefinite (exact check).  Weak duality: for every feasible (b, Q),
        b <= sum_H y_H (d(H) - <Q, M(H)>) = sum_H y_H d(H) - sum_sigma <Q_sigma, Y_sigma> <= sum_H y_H d(H),
    so NO plain N = 6 certificate (any blocks with 2m - s <= 6) can prove more than sum_H y_H d(H).
    y is reconstructed from the SDP dual by continued fractions (limit_denominator) and then checked exactly.

Usage: python n6_control.py [--p 5 6 7]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from fractions import Fraction
from math import comb

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
sys.path.insert(0, PARENT)

from exact import certify_eig, certify_sharp, is_psd_exact  # noqa: E402
from flags import build_problem  # noqa: E402
from hypergraphs import EveryPSetHasEdge, PSetEdgeCounts, mask_to_edges  # noqa: E402
from sdp import solve_flag_sdp  # noqa: E402
from turan import write_cert  # noqa: E402

RES = os.path.join(HERE, "results")


def moment_matrices_exact(prob, y):
    """Y_sigma = sum_H y_H M_sigma(H) as Fraction matrices."""
    out = []
    for blk in prob.blocks:
        n = blk.size
        Y = [[Fraction(0)] * n for _ in range(n)]
        for h, a, b, c in zip(blk.h.tolist(), blk.a.tolist(), blk.b.tolist(), blk.cnt.tolist()):
            if y[h]:
                Y[a][b] += y[h] * Fraction(c, blk.denom)
        out.append(Y)
    return out


def rational_dual(prob, yf, dens, maxden):
    y = [Fraction(float(max(v, 0.0))).limit_denominator(maxden) for v in yf]
    y = [v if v > Fraction(1, 10 * maxden) else Fraction(0) for v in y]
    s = sum(y)
    if s == 0:
        return None
    y = [v / s for v in y]
    Ys = moment_matrices_exact(prob, y)
    psd = [is_psd_exact(Y) for Y in Ys]
    val = sum(v * dens[i] for i, v in enumerate(y))
    return {"y": y, "value": val, "psd_all": all(ok for ok, _, _ in psd), "ranks": [rk for _, rk, _ in psd],
            "maxden": maxden}


def run_p(p, N=6, r=4, lam=1, target=None):
    t0 = time.time()
    pred = EveryPSetHasEdge(p) if lam == 1 else PSetEdgeCounts(p, range(lam, comb(p, r) + 1))
    prob = build_problem(N, r, pred, type_mode="all")
    num, den = prob.density_num_den()
    dens = [Fraction(int(x), den) for x in num]
    target = Fraction(1, comb(p - 1, r - 1)) if target is None else Fraction(target)
    out = {"p": p, "N": N, "r": r, "lam": lam, "n_admissible": prob.nH, "target_deCaen": str(target),
           "blocks": [(b.fam.s, b.fam.m, b.size) for b in prob.blocks],
           "min_H_density_noSOS": str(min(dens))}
    sol = solve_flag_sdp(prob, num / den)
    out["sdp"] = {k: sol[k] for k in ("status", "b", "obj_dual", "iterations", "r_prim", "r_dual")}
    ce = certify_eig(prob, dens, sol["Q"], bits=40)
    out["eig_b_exact_float"] = float(ce["b_exact"])
    slf = np.array([float(x) for x in ce["slacks"]]) - float(ce["b_exact"])
    cs = certify_sharp(prob, dens, sol["Q"], slf, target, verbose=False)
    out["sharp"] = {"success": cs.get("success"), "b_exact": str(cs.get("b_exact")), "psd_ok": cs.get("psd_ok"),
                    "reason": cs.get("reason")}
    tag = "" if lam == 1 else f"_l{lam}"
    cert_path = os.path.join(RES, f"cert_turan_r4_p{p}{tag}_N6_all_sharp.json")
    if cs.get("success"):
        write_cert(cert_path, prob, cs["Qexact"], cs["b_exact"], r, p, N,
                   None if lam == 1 else list(range(lam, comb(p, r) + 1)))
        out["sharp"]["cert"] = os.path.relpath(cert_path, HERE)
    # dual (upper bound on the N = 6 optimum)
    best = None
    for maxden in (64, 256, 1024, 4096, 1 << 16):
        rd = rational_dual(prob, sol["y"], dens, maxden)
        if rd and rd["psd_all"]:
            best = rd
            break
    if best is None:
        out["dual"] = {"success": False}
    else:
        supp = [(i, v) for i, v in enumerate(best["y"]) if v]
        out["dual"] = {"success": True, "value": str(best["value"]), "value_float": float(best["value"]),
                       "equals_target": best["value"] == target, "maxden": best["maxden"],
                       "moment_ranks": best["ranks"],
                       "support": [{"y": str(v), "edges": int(num[i]),
                                    "graph_edges": mask_to_edges(int(prob.H[i]), N, r)} for i, v in supp]}
    out["time_s"] = time.time() - t0
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", type=int, nargs="+", default=[5, 6, 7])
    ap.add_argument("--lam", type=int, default=1, help="every p-set spans >= lam edges (lam = 2, p = 5: K_5^{4-})")
    ap.add_argument("--target", default=None, help="sharp target (default de Caen 1/C(p-1,3) for lam = 1)")
    a = ap.parse_args()
    os.makedirs(RES, exist_ok=True)
    allres = {}
    for p in a.p:
        res = run_p(p, lam=a.lam, target=a.target)
        allres[p] = res
        print(json.dumps(res, default=str), flush=True)
        tag = "" if a.lam == 1 else f"_l{a.lam}"
        with open(os.path.join(RES, f"n6_control_p{p}{tag}.json"), "w") as f:
            json.dump(res, f, indent=1, default=str)


if __name__ == "__main__":
    main()
