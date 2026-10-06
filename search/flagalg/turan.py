"""Flag-algebra lower bounds on the Turan density in complement form,
    t(p, r) = lim min{ e(G)/C(n, r) : every p-set of V(G) contains an edge } = 1 - pi(K_p^{(r)}).

Usage:  python turan.py --r 3 --p 4 --N 6 [--types standard|all] [--sharp] [--out results/x.json]
"""
from __future__ import annotations

import argparse
import json
import resource
import time
from fractions import Fraction
from math import comb

import numpy as np

from exact import certify_eig, certify_sharp
from flags import build_problem
from hypergraphs import EveryPSetHasEdge, PSetEdgeCounts, mask_to_edges
from sdp import solve_flag_sdp


def peak_rss_mb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def run(r: int, p: int, N: int, types: str = "standard", sharp: bool = False, target=None, bits: int = 40,
        tol: float = 1e-10, verbose: bool = True, cert_out: str | None = None, allowed=None,
        settings: dict | None = None) -> dict:
    t0 = time.time()
    pred = EveryPSetHasEdge(p) if allowed is None else PSetEdgeCounts(p, allowed)
    prob = build_problem(N, r, pred, type_mode=types, verbose=verbose)
    t_build = time.time() - t0
    num, den = prob.density_num_den()
    d = num / den
    t1 = time.time()
    sol = solve_flag_sdp(prob, d, tol=tol, settings=settings)
    t_sdp = time.time() - t1
    if verbose:
        print(f"  SDP: {sol['status']} b = {sol['b']:.12f} ({sol['iterations']} it, {t_sdp:.2f}s)", flush=True)
    obj_exact = [Fraction(int(x), den) for x in num]
    t2 = time.time()
    ce = certify_eig(prob, obj_exact, sol["Q"], bits=bits)
    t_eig = time.time() - t2
    if verbose:
        print(f"  eig rounding: b_exact = {float(ce['b_exact']):.12f}  PSD checks {ce['psd']}", flush=True)
    out = {
        "problem": f"minimise edge density of r-graphs with: {pred.name}",
        "r": r, "p": p, "N": N, "types": types,
        "classes_per_n": {n: int(v.size) for n, v in prob.classes.items()},
        "n_admissible": prob.nH,
        "blocks": [{"s": b.fam.s, "m": b.fam.m, "sigma_edges": mask_to_edges(b.fam.sigma, b.fam.s, r),
                    "n_flags": b.size, "denominator": b.denom, "nnz": int(b.cnt.size)} for b in prob.blocks],
        "sdp": {k: sol[k] for k in ("status", "b", "iterations", "solve_time", "r_prim", "r_dual", "obj_dual",
                                    "n_vars", "n_rows", "nnz")},
        "eig_rounding": {"bits": bits, "b_exact": str(ce["b_exact"]), "b_exact_float": float(ce["b_exact"]),
                         "psd_checks": [{"psd": ok, "rank": rk} for ok, rk in ce["psd"]],
                         "argmin_graph_edges": mask_to_edges(int(prob.H[ce["argmin"]]), N, r),
                         "gap_to_float": float(sol["b"]) - float(ce["b_exact"])},
        "time_s": {"build": t_build, "sdp": t_sdp, "eig_rounding": t_eig},
        "min_H_density": str(Fraction(int(num.min()), den)),
    }
    y = sol["y"]
    top = np.argsort(-y)[:8]
    out["dual_support"] = [{"y": float(y[i]), "edges": int(num[i]), "density": float(d[i])} for i in top
                           if y[i] > 1e-6]
    if sharp:
        tgt = Fraction(target) if target is not None else Fraction(sol["b"]).limit_denominator(1000)
        slf = np.array([float(x) for x in ce["slacks"]]) - float(ce["b_exact"])
        # float slacks relative to the float bound (d(H) - <Q, M(H)> - b)
        t3 = time.time()
        cs = certify_sharp(prob, obj_exact, sol["Q"], slf, tgt, verbose=verbose)
        out["sharp_rounding"] = {"target": str(tgt), "success": cs.get("success"), "n_sharp": cs.get("n_sharp"),
                                 "b_exact": str(cs.get("b_exact")), "psd_ok": cs.get("psd_ok"),
                                 "X_rank": cs.get("X_rank"), "reason": cs.get("reason"),
                                 "warnings": cs.get("warnings"), "time_s": time.time() - t3}
        if verbose:
            print(f"  sharp rounding: {out['sharp_rounding']}", flush=True)
        if cert_out and cs.get("success"):
            write_cert(cert_out, prob, cs["Qexact"], cs["b_exact"], r, p, N, allowed)
            cert_out = None
    if cert_out:
        write_cert(cert_out, prob, ce["Qexact"], ce["b_exact"], r, p, N, allowed)
    out["peak_rss_mb"] = peak_rss_mb()
    out["time_s"]["total"] = time.time() - t0
    return out


def write_cert(path, prob, Qexact, b, r, p, N, allowed=None):
    cert = {"r": r, "p": p, "N": N, "bound": str(b), "allowed": list(allowed) if allowed else None,
            "blocks": [{"s": blk.fam.s, "m": blk.fam.m, "sigma": int(blk.fam.sigma),
                        "flags": [int(x) for x in blk.fam.flags], "den": str(den),
                        "Qnum": [[str(int(x)) for x in row] for row in Qnum.tolist()]}
                       for blk, (Qnum, den) in zip(prob.blocks, Qexact)]}
    with open(path, "w") as f:
        json.dump(cert, f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r", type=int, required=True)
    ap.add_argument("--p", type=int, required=True)
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--types", default="standard", choices=["standard", "all", "none"])
    ap.add_argument("--allowed", default=None, help="comma list: allowed edge counts in every p-set (induced problems)")
    ap.add_argument("--static-reg", type=float, default=None, help="clarabel static_regularization_constant")
    ap.add_argument("--sharp", action="store_true")
    ap.add_argument("--target", default=None)
    ap.add_argument("--bits", type=int, default=40)
    ap.add_argument("--tol", type=float, default=1e-10)
    ap.add_argument("--out", default=None)
    ap.add_argument("--cert", default=None)
    a = ap.parse_args()
    print(f"# turan r={a.r} p={a.p} N={a.N} types={a.types}", flush=True)
    allowed = [int(x) for x in a.allowed.split(",")] if a.allowed else None
    st = {"static_regularization_constant": a.static_reg} if a.static_reg else None
    res = run(a.r, a.p, a.N, a.types, a.sharp, a.target, a.bits, a.tol, cert_out=a.cert, allowed=allowed,
              settings=st)
    js = json.dumps(res, indent=1, default=str)
    if a.out:
        with open(a.out, "w") as f:
            f.write(js)
    print(js)


if __name__ == "__main__":
    main()
