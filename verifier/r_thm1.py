"""Independent checker (Theorem 1): the explicit O(1/n) error term c(n) of the robust consequence, from the certificate
file only (own code; no producer module imported).

Q_k = G_k / D with G_k = q W_k^T W_k + 2^(2 kb) B_k^T Xq_k B_k, D = 2^(2 kb) q (the formula of certificates/sigma_K5_4).
For each key k: lambda_k = max_i (Q_k)_ii exactly (diag(W^T W)_i = sum_r W_ri^2, diag(B^T X B)_i = b_i^T X b_i).
Also: exact max |(Q_k)_ij| on the small keys (all entries), and a PSD-implied bound |Q_ij| <= max diag elsewhere.
Error term per key: lambda_k * p_k/(1-p_k), p_k = P[two independent uniform (m-s)-subsets of an (n-s)-set meet].
Checks: (a) the producer's per-key max_diag fractions equal ours (only where its results file is available;
not in the package); (b) the per-block
sums 61.2022134, 0.9434973, 1.4007420, 0.0002868; (c) c(n) <= 68/(n-7) exactly for 10 <= n <= 20000 and the
analytic per-block bounds; (d) the sharper constant obtained with the max (not the sum) over the types of one block;
(e) mu, 1/mu, K = 21894/mu, the Theorem 3 constants.
Output: results_thm1.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from math import comb

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.getcwd()   # outputs (results_thm1.json) go to the working directory
PKG = os.path.join(os.path.dirname(HERE), "certificates", "sigma_K5_4")
CERT = os.path.join(PKG, "sharp_klp5.cert.npz")
KEYS = os.path.join(PKG, "sharp_klp5.cert.json")
PROD = None   # the producer's results file is not part of the package; comparison (a) is skipped
MU = Fraction(176165518826891630107, 6280747422216628134215680)   # the margin; minimum checked by the raw scan (rs_eval.c)


def main():
    z = np.load(CERT)
    kb = int(z["kb"])
    q = int(str(z["q"]))
    D = int(str(z["D"]))
    assert D == (1 << (2 * kb)) * q
    corr = {c["key"]: c for c in json.loads(str(z["corr"]))}
    keys = json.load(open(KEYS))["keys"]
    res = {"kb": kb, "q": str(q), "D": str(D), "keys": []}
    for k, kd in enumerate(keys):
        nf = kd["n_flags"]
        W = z[f"W{k}"]
        Wi = [[int(str(x)) for x in row] for row in W] if W.size else []
        assert all(len(r) == nf for r in Wi)
        diag = [q * sum(r[i] * r[i] for r in Wi) for i in range(nf)]
        if k in corr:
            B = z[f"B{k}"].astype(object)
            X = [[int(x) for x in row] for row in corr[k]["X_times_q"]]
            r_ = len(X)
            assert B.shape == (r_, nf)
            Bi = [[int(B[a, i]) for i in range(nf)] for a in range(r_)]
            sym = all(X[a][b] == X[b][a] for a in range(r_) for b in range(r_))
            assert sym
            for i in range(nf):
                col = [Bi[a][i] for a in range(r_)]
                diag[i] += (1 << (2 * kb)) * sum(col[a] * X[a][b] * col[b] for a in range(r_) for b in range(r_))
        md = max(diag) if diag else 0
        ent = {"key": k, "s": kd["s"], "m": kd["m"], "sigma": kd["sigma"], "n_flags": nf,
               "lambda": str(Fraction(md, D)), "lambda_float": float(Fraction(md, D)),
               "min_diag_nonneg": min(diag) >= 0 if diag else True}
        # exact max |entry| for small keys (full G), else PSD bound
        if nf <= 40:
            G = [[0] * nf for _ in range(nf)]
            for r in Wi:
                for i in range(nf):
                    if r[i]:
                        for j in range(nf):
                            G[i][j] += q * r[i] * r[j]
            if k in corr:
                for i in range(nf):
                    for j in range(nf):
                        G[i][j] += (1 << (2 * kb)) * sum(Bi[a][i] * X[a][b] * Bi[b][j]
                                                         for a in range(r_) for b in range(r_))
            mx = max(abs(G[i][j]) for i in range(nf) for j in range(nf))
            ent["max_abs_entry_equals_lambda"] = mx == md
        res["keys"].append(ent)
    if PROD is not None:
        prod = json.load(open(PROD))
        pk = {e["key"]: Fraction(e["max_diag"]) for e in prod["keys"]}
        res["producer_max_diag_equal"] = all(Fraction(e["lambda"]) == pk[e["key"]] for e in res["keys"])
    blocks = {}
    for e in res["keys"]:
        blocks.setdefault((e["s"], e["m"]), []).append(Fraction(e["lambda"]))
    res["block_sum"] = {f"{s},{m}": float(sum(v)) for (s, m), v in blocks.items()}
    res["block_max"] = {f"{s},{m}": float(max(v)) for (s, m), v in blocks.items()}

    def pratio(n, s, m):
        kk = m - s
        disj = Fraction(comb(n - m, kk), comb(n - s, kk))   # P[U1, U2 disjoint]
        return (1 - disj) / disj

    def c(n, agg):
        return sum((agg(v) * pratio(n, s, m) for (s, m), v in blocks.items()), Fraction(0))

    # analytic per-block bounds  p/(1-p) <= w_b/(n-7), checked exactly
    wts = {(5, 6): 1, (4, 5): 1, (3, 5): 4, (1, 4): 12, (2, 4): 0, (3, 4): 0}
    an_ok = True
    for n in range(10, 20001):
        for (s, m), w in wts.items():
            if any(x != 0 for x in blocks[(s, m)]):
                an_ok &= pratio(n, s, m) <= Fraction(w, n - 7)
    # formulas for p/(1-p)
    form_ok = all(pratio(n, 5, 6) == Fraction(1, n - 6) and pratio(n, 4, 5) == Fraction(1, n - 5)
                  and pratio(n, 3, 5) == Fraction(4 * n - 18, (n - 5) * (n - 6))
                  and pratio(n, 1, 4) == Fraction(9 * n * n - 63 * n + 114, (n - 4) * (n - 5) * (n - 6))
                  for n in range(10, 2001))
    const = sum((sum(v) * wts[b] for b, v in blocks.items()), Fraction(0))
    const_max = sum((max(v) * wts[b] for b, v in blocks.items()), Fraction(0))
    cn_ok = all(c(n, sum) <= Fraction(68, n - 7) for n in range(10, 20001))
    res["c_n"] = {"p_over_1mp_formulas_ok_10..2000": form_ok, "analytic_block_bounds_ok_10..20000": an_ok,
                  "constant_sum_form": float(const), "constant_sum_form_le_68": const <= 68,
                  "c(n)<=68/(n-7)_exact_10..20000": cn_ok, "c(10)": float(c(10, sum)),
                  "n*c(n)_n=20000": float(20000 * c(20000, sum)),
                  "constant_with_max_over_types": float(const_max),
                  "c_max(n)<=27/(n-7)_10..5000": all(c(n, max) <= Fraction(27, n - 7) for n in range(10, 5001))}
    # mu, K and Theorem 3 numbers
    K = Fraction(21894) / MU
    rho_num = 1 - Fraction(5, 3) * (Fraction(33, 64) + Fraction(1, 100) - Fraction(2, 5)) - Fraction(1, 10)
    res["thm3"] = {"mu": str(MU), "mu_float": float(MU), "inv_mu": float(1 / MU), "inv_mu_le_35653": 1 / MU <= 35653,
                   "K=21894/mu": float(K), "K_le_7.81e8": K <= Fraction(781, 10 ** 6) * 10 ** 12,
                   "rho_factor": str(rho_num), "120*126/rho_factor": float(Fraction(120 * 126) / rho_num),
                   "120*126/rho_factor_le_21894": Fraction(120 * 126) / rho_num <= 21894,
                   "mu/1200_implies_f_le_33/64+1/100": MU / 1200 <= Fraction(1, 100),
                   "1200/mu_le_K": Fraction(1200) / MU <= K,
                   "frames_min_over_a": min(comb(a, 5) * 120 * (10 - a) + comb(10 - a, 5) * 120 * a
                                            for a in range(5, 10)),
                   "(10)_6": 151200,
                   "n_where_69K/(n-7)<1": float(69 * K + 7),
                   "n_where_68/(n-7)<mu": float(68 / MU + 7)}
    json.dump(res, open(os.path.join(OUT, "results_thm1.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
