"""Primal-dual interior-point method (HKM direction, Mehrotra predictor-corrector) for the dual flag-algebra SDP.
Float64; produces a strictly feasible interior pseudo-density y (and an approximate certificate (b, Q)); nothing is
certified here -- the exact check is verify_dual.py.

    (D)  min d.y   s.t.  Y_k(y) = sum_i y_i M_i^k  PSD (live flags),  y >= 0,  1.y = 1
    (P)  max b     s.t.  z = d - b 1 - A^T(Q) >= 0,  Q_k PSD,      A^T(Q)_i = sum_k <M_i^k, Q_k>.
Central path Y_k Q_k = mu I, y_i z_i = mu.  The y-iterates stay strictly feasible (Y_k > 0, y > 0, sum y = 1
exactly up to rounding); (b, Q, z) may be infeasible at the start (residual r_d -> 0).
Schur matrix S_ij = sum_k tr(M_i Z_k M_j Q_k) + delta_ij z_i / y_i,  computed with the Aut-invariance trick of
barrier.py: column i = |Aut| * S_k @ vec(Z[:, R_i] C_i Q[R_i, :]).

Usage: python pdip.py results/support_X.npz --start results/X_phase1.npz --out results/pdip_X.npz [--gap 1e-9]
"""
import argparse
import json
import time

import numpy as np
import scipy.linalg as sla

from barrier import Problem, phase1


def schur_key(K, Z, Q, m, chunk=64):
    S = np.zeros((m, m))
    n = K.n
    for c0 in range(0, len(K.blocks), chunk):
        blk = K.blocks[c0:c0 + chunk]
        B = np.empty((n * n, len(blk)))
        for q, (Rj, Cj) in enumerate(blk):
            B[:, q] = (Z[:, Rj] @ Cj @ Q[Rj, :]).ravel()
        S[:, K.graphs[c0:c0 + len(blk)]] = K.G * (K.S @ B)
    return 0.5 * (S + S.T)


def sym(X):
    return 0.5 * (X + X.T)


def group_avg(K, X):
    out = np.zeros_like(X)
    for g in range(K.G):
        p = K.P[g]
        out += X[np.ix_(p, p)]
    return out / K.G


def max_step(X, dX):
    """largest alpha with X + alpha dX PSD (X PD)."""
    L = np.linalg.cholesky(X)
    W = sla.solve_triangular(L, sla.solve_triangular(L, dX, lower=True).T, lower=True)
    lam = float(np.linalg.eigvalsh(sym(W))[0])
    return np.inf if lam >= 0 else -1.0 / lam


def max_step_vec(x, dx):
    neg = dx < 0
    return np.inf if not neg.any() else float(np.min(-x[neg] / dx[neg]))


def AT(pb, Xs):
    out = np.zeros(pb.m)
    for K, X in zip(pb.keys, Xs):
        out += K.grad(X)               # S_k @ vec(X) = <M_i, X>
    return out


def run(pb, y, log, gap=1e-9, max_it=200, tau=0.95, save=None, mu0=None, save_best=None):
    m = pb.m
    nu = pb.nu_psd + m
    Ys = [K.Y(y) for K in pb.keys]
    Zs = [np.linalg.inv(Y) for Y in Ys]
    mu = mu0 if mu0 is not None else 1e-2
    Qs = [group_avg(K, sym(mu * Z)) for K, Z in zip(pb.keys, Zs)]
    z = mu / y
    b = float(np.min(pb.d - AT(pb, Qs) - z)) - 1e-3
    one = np.ones(m)
    hist = []
    best_obj = np.inf
    for it in range(1, max_it + 1):
        tt = time.time()
        Ys = [K.Y(y) for K in pb.keys]
        Zs = [np.linalg.inv(Y) for Y in Ys]
        Zs = [sym(Z) for Z in Zs]
        rd = pb.d - b * one - AT(pb, Qs) - z
        rp = 1.0 - y.sum()
        comp = sum(float(np.sum(Y * Q)) for Y, Q in zip(Ys, Qs)) + float(y @ z)
        mu = comp / nu
        row = {"it": it, "obj_y": float(pb.d @ y), "b": b, "mu": mu, "nu_mu": nu * mu,
               "rd": float(np.abs(rd).max()), "rp": rp}
        if nu * mu < gap and row["rd"] < 1e-12:
            log(row)
            break
        S = np.diag(z / y)
        for K, Z, Q in zip(pb.keys, Zs, Qs):
            S += schur_key(K, Z, Q, m)
        try:
            cf = sla.cho_factor(S, lower=True, check_finite=False)
            solve = lambda r: sla.cho_solve(cf, r, check_finite=False)  # noqa: E731
        except np.linalg.LinAlgError:
            lu = sla.lu_factor(S)
            solve = lambda r: sla.lu_solve(lu, r)  # noqa: E731
        S1 = solve(one)

        def direction(Cs, c):
            # Cs: per key target C_k (sigma mu I - Y Q - corr); c: sigma mu - y z - corr
            h = c / y + AT(pb, [sym(Z @ C) for Z, C in zip(Zs, Cs)]) - rd
            Sh = solve(h)
            db = (rp - one @ Sh) / (one @ S1)
            dy = Sh + db * S1
            dYs = [K.Y(dy) for K in pb.keys]
            dQs = [sym(Z @ C) - sym(Z @ dY @ Q) for Z, C, dY, Q in zip(Zs, Cs, dYs, Qs)]
            dz = (c - z * dy) / y
            return dy, db, dYs, dQs, dz

        def steps(dy, dYs, dQs, dz):
            ap = min([max_step_vec(y, dy)] + [max_step(Y, dY) for Y, dY in zip(Ys, dYs)])
            ad = min([max_step_vec(z, dz)] + [max_step(Q, dQ) for Q, dQ in zip(Qs, dQs)])
            return ap, ad

        # predictor
        Cs = [-(Y @ Q) for Y, Q in zip(Ys, Qs)]
        dy, db, dYs, dQs, dz = direction(Cs, -y * z)
        ap, ad = steps(dy, dYs, dQs, dz)
        ap, ad = min(1.0, ap), min(1.0, ad)
        comp_aff = sum(float(np.sum((Y + ap * dY) * (Q + ad * dQ))) for Y, dY, Q, dQ in zip(Ys, dYs, Qs, dQs)) + \
            float((y + ap * dy) @ (z + ad * dz))
        sigma = min(1.0, (comp_aff / comp) ** 3)
        # corrector (Mehrotra second-order term)
        Cs = [sigma * mu * np.eye(Y.shape[0]) - Y @ Q - dY @ dQ for Y, Q, dY, dQ in zip(Ys, Qs, dYs, dQs)]
        c = sigma * mu - y * z - dy * dz
        dy, db, dYs, dQs, dz = direction(Cs, c)
        ap, ad = steps(dy, dYs, dQs, dz)
        ap, ad = min(1.0, tau * ap), min(1.0, tau * ad)
        y = y + ap * dy
        b = b + ad * db
        Qs = [group_avg(K, sym(Q + ad * dQ)) for K, Q, dQ in zip(pb.keys, Qs, dQs)]
        z = z + ad * dz
        row.update({"sigma": sigma, "ap": ap, "ad": ad, "sec": round(time.time() - tt, 2)})
        log(row)
        hist.append(row)
        if save is not None:
            save(y, b, Qs, z)
        obj = float(pb.d @ y)
        if save_best is not None and obj < best_obj:
            best_obj = obj
            save_best(y, it)
    return y, b, Qs, z, hist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("support")
    ap.add_argument("--start", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--gap", type=float, default=1e-9)
    ap.add_argument("--max-it", type=int, default=200)
    ap.add_argument("--tau", type=float, default=0.95)
    a = ap.parse_args()
    t0 = time.time()
    pb = Problem(a.support)
    print(json.dumps({"m": pb.m, "dims": [K.n for K in pb.keys], "setup_s": round(time.time() - t0, 1)}),
          flush=True)

    def log(r):
        print(json.dumps(r), flush=True)

    if a.start:
        st = np.load(a.start)
        assert np.array_equal(st["masks"], pb.masks)
        y = st["y"]
    else:
        y, s = phase1(pb, log)
        assert s < 0, "no strictly feasible point"
        np.savez(a.out.replace(".npz", "_phase1.npz"), y=y, s=s, masks=pb.masks)
    y = y / y.sum()

    def save(y, b, Qs, z):
        np.savez(a.out, y=y, b=b, z=z, masks=pb.masks, d=pb.d, **{f"Q{i}": Q for i, Q in enumerate(Qs)})

    def save_best(y, it):
        np.savez(a.out.replace(".npz", "_best.npz"), y=y, it=it, masks=pb.masks, d=pb.d)

    y, b, Qs, z, hist = run(pb, y, log, gap=a.gap, max_it=a.max_it, tau=a.tau, save=save, save_best=save_best)
    save(y, b, Qs, z)
    log({"final_obj_y": float(pb.d @ y), "b": b, "min_eigs": pb.min_eigs(y), "rel_min_eigs": pb.rel_min_eigs(y),
         "min_y": float(y.min()), "time_s": round(time.time() - t0, 1)})


if __name__ == "__main__":
    main()
