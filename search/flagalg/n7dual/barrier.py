"""Primal barrier method for the dual flag-algebra SDP (float64; produces an interior point, nothing certified here).

    minimise  d . y   over  y >= 0, sum y = 1 (y on the alive support S),
    s.t.      Y_k(y) = sum_H y_H M_k(H)  PSD  for every key k (restricted to the live flags; dead rows are exactly 0).

Phase I:  minimise s  s.t. Y_k(y) + s E_k > 0 (E_k = Diag Y_k(uniform) > 0 on live flags); stops as soon as s < 0.
Phase II: barrier  F_t(y) = t d.y - sum_k log det Y_k(y) - sum_H log y_H,  t increased geometrically.
Newton steps with the equality constraint sum y = 1 (KKT by elimination), backtracking line search keeping every
Y_k positive definite (Cholesky) and y > 0.

Hessian: H_k[i, j] = tr(Z M_i Z M_j), Z = Y_k^{-1}.  Z is Aut(sigma)-invariant, so
    tr(Z M_i Z M_j) = |Aut| tr(Z Mt_i Z M_j),   Z Mt_j Z = Z[:, R_j] C_j Z[R_j, :]
with R_j the (few) flags of the canonical terms of graph j; then H_k[:, j] = |Aut| * S_k @ vec(Z Mt_j Z), S_k the
sparse (graphs x n^2) matrix of the symmetrised M_i.

Usage: python barrier.py results/support_X.npz --out results/barrier_X.npz [--gap 1e-9]
"""
import argparse
import json
import time

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp


class Key:
    def __init__(self, data, k, m):
        self.n = n = int(data[f"live{k}"].size)
        self.P = data[f"P{k}"].astype(np.int64)
        self.G = self.P.shape[0]
        ti, ta, tb = data[f"ti{k}"], data[f"ta{k}"], data[f"tb{k}"]
        c = data[f"tw{k}"].astype(np.float64) / 10080.0
        G = self.G
        R = np.tile(ti, G)
        X = self.P[:, ta].ravel()
        Yb = self.P[:, tb].ravel()
        C = np.tile(c, G)
        rows = np.concatenate([R, R])
        cols = np.concatenate([X * n + Yb, Yb * n + X])
        vals = np.concatenate([C, C])
        self.S = sp.csr_matrix((vals, (rows, cols)), shape=(m, n * n))
        self.ST = self.S.T.tocsr()
        # canonical per-graph small blocks
        order = np.argsort(ti, kind="stable")
        ti, ta, tb, c = ti[order], ta[order], tb[order], c[order]
        bounds = np.searchsorted(ti, np.arange(m + 1))
        self.blocks = []
        self.graphs = []
        for j in range(m):
            lo, hi = bounds[j], bounds[j + 1]
            if lo == hi:
                continue
            a, b, cc = ta[lo:hi], tb[lo:hi], c[lo:hi]
            Rj, inv = np.unique(np.concatenate([a, b]), return_inverse=True)
            ia, ib = inv[:a.size], inv[a.size:]
            Cj = np.zeros((Rj.size, Rj.size))
            np.add.at(Cj, (ia, ib), cc)
            np.add.at(Cj, (ib, ia), cc)
            self.blocks.append((Rj, Cj))
            self.graphs.append(j)
        self.graphs = np.array(self.graphs, dtype=np.int64)

    def Y(self, y):
        return (self.ST @ y).reshape(self.n, self.n)

    def grad(self, Z):
        return self.S @ Z.ravel()

    def hess(self, Z, m, chunk=64):
        """H (m x m) = [tr(Z M_i Z M_j)]."""
        Hk = np.zeros((m, m))
        n = self.n
        for c0 in range(0, len(self.blocks), chunk):
            blk = self.blocks[c0:c0 + chunk]
            B = np.empty((n * n, len(blk)))
            for q, (Rj, Cj) in enumerate(blk):
                ZR = Z[:, Rj]
                B[:, q] = (ZR @ Cj @ ZR.T).ravel()
            Hk[:, self.graphs[c0:c0 + len(blk)]] = self.G * (self.S @ B)
        return 0.5 * (Hk + Hk.T)


class Problem:
    def __init__(self, path):
        data = np.load(path)
        self.masks = data["masks"]
        self.m = m = self.masks.size
        self.d = data["e"].astype(np.float64) / 35.0
        self.keys = [Key(data, k, m) for k in range(int(data["nkeys"]))]
        self.keys = [K for K in self.keys if K.n > 0]
        self.nu_psd = sum(K.n for K in self.keys)

    def chol_all(self, y, s=0.0, E=None):
        out = []
        for i, K in enumerate(self.keys):
            Y = K.Y(y)
            if E is not None:
                Y = Y + s * E[i]
            try:
                L = np.linalg.cholesky(Y)
            except np.linalg.LinAlgError:
                return None
            out.append(L)
        return out

    @staticmethod
    def logdet(L):
        return 2.0 * np.log(np.diag(L)).sum()

    @staticmethod
    def inv_from_chol(L):
        Li = sla.solve_triangular(L, np.eye(L.shape[0]), lower=True)
        return Li.T @ Li

    def min_eigs(self, y):
        return [float(np.linalg.eigvalsh(K.Y(y))[0]) for K in self.keys]

    def rel_min_eigs(self, y):
        """min eigenvalue of D^-1/2 Y D^-1/2 (D = diag Y) per key: scale-free margin."""
        out = []
        for K in self.keys:
            Y = K.Y(y)
            dg = np.sqrt(np.diag(Y))
            out.append(float(np.linalg.eigvalsh(Y / np.outer(dg, dg))[0]))
        return out


def newton_dir(H, g, a):
    """min g.dx + dx H dx / 2 s.t. a.dx = 0."""
    try:
        c = sla.cho_factor(H, lower=True, check_finite=False)
        Hg = sla.cho_solve(c, g, check_finite=False)
        Ha = sla.cho_solve(c, a, check_finite=False)
    except np.linalg.LinAlgError:
        w, V = np.linalg.eigh(H)
        w = np.maximum(w, 1e-14 * w.max())
        Hg = V @ ((V.T @ g) / w)
        Ha = V @ ((V.T @ a) / w)
    wq = -(a @ Hg) / (a @ Ha)
    dx = -(Hg + wq * Ha)
    return dx, float(-(g @ dx))


def phase1(pb, log, max_newton=400):
    m = pb.m
    y = np.full(m, 1.0 / m)
    E = [np.diag(np.diag(K.Y(y))) for K in pb.keys]
    lam = min(float(sla.eigh(K.Y(y), np.diag(np.diag(K.Y(y))), eigvals_only=True)[0]) for K in pb.keys)
    s = max(0.0, -lam) + 0.5
    SLB = 1.0
    t = 2.0 * pb.nu_psd / s
    a = np.concatenate([np.ones(m), [0.0]])
    it = 0

    def F(y, s, Ls):
        return t * s - sum(pb.logdet(L) for L in Ls) - np.log(y).sum() - np.log(s + SLB)

    Ls = pb.chol_all(y, s, E)
    log({"phase": 1, "start_s": s, "lam_min_rel": lam})
    while it < max_newton:
        it += 1
        Zs = [pb.inv_from_chol(L) for L in Ls]
        g = np.zeros(m + 1)
        H = np.zeros((m + 1, m + 1))
        for K, Z, Ek in zip(pb.keys, Zs, E):
            g[:m] -= K.grad(Z)
            g[m] -= float(np.sum(Z * Ek))
            H[:m, :m] += K.hess(Z, m)
            ZEZ = Z @ Ek @ Z
            hs = K.grad(ZEZ)
            H[:m, m] += hs
            H[m, :m] += hs
            H[m, m] += float(np.sum(ZEZ * Ek))
        g[:m] -= 1.0 / y
        g[m] += t - 1.0 / (s + SLB)
        H[np.arange(m), np.arange(m)] += 1.0 / y ** 2
        H[m, m] += 1.0 / (s + SLB) ** 2
        dx, dec = newton_dir(H, g, a)
        f0 = F(y, s, Ls)
        al = 1.0
        neg = dx[:m] < 0
        if neg.any():
            al = min(al, 0.99 * float(np.min(-y[neg] / dx[:m][neg])))
        if dx[m] < 0:
            al = min(al, 0.99 * (s + SLB) / -dx[m])
        while True:
            yn, sn = y + al * dx[:m], s + al * dx[m]
            Ln = pb.chol_all(yn, sn, E)
            if Ln is not None and F(yn, sn, Ln) <= f0 - 0.25 * al * dec:
                break
            al *= 0.5
            if al < 1e-12:
                break
        if Ln is None:
            log({"phase": 1, "it": it, "error": "line search failed"})
            break
        y, s, Ls = yn / yn.sum(), sn, Ln
        Ls = pb.chol_all(y, s, E)
        log({"phase": 1, "it": it, "t": t, "s": s, "dec": dec, "step": al, "obj": float(pb.d @ y)})
        if s < 0:
            return y, s
        if dec < 1e-2:
            t *= 4.0
    return y, s


def phase2(pb, y, log, gap=1e-9, mu=4.0, t0=None, max_newton=2000, save=None):
    m = pb.m
    nu = pb.nu_psd + m
    a = np.ones(m)
    t = t0 if t0 is not None else nu / max(float(pb.d @ y) - 0.25, 1e-3)
    it = 0
    hist = []

    def F(y, Ls):
        return t * float(pb.d @ y) - sum(pb.logdet(L) for L in Ls) - np.log(y).sum()

    Ls = pb.chol_all(y)
    assert Ls is not None, "phase 2 needs a strictly feasible start"
    while it < max_newton:
        it += 1
        tt = time.time()
        Zs = [pb.inv_from_chol(L) for L in Ls]
        g = t * pb.d - 1.0 / y
        H = np.diag(1.0 / y ** 2)
        for K, Z in zip(pb.keys, Zs):
            g -= K.grad(Z)
            H += K.hess(Z, m)
        dx, dec = newton_dir(H, g, a)
        f0 = F(y, Ls)
        al = 1.0
        neg = dx < 0
        if neg.any():
            al = min(al, 0.99 * float(np.min(-y[neg] / dx[neg])))
        Ln = None
        while al > 1e-14:
            yn = y + al * dx
            Ln = pb.chol_all(yn)
            if Ln is not None and F(yn, Ln) <= f0 - 0.25 * al * dec:
                break
            Ln = None
            al *= 0.5
        if Ln is None:
            log({"phase": 2, "it": it, "error": "line search failed", "dec": dec})
            break
        y = yn
        Ls = Ln
        row = {"phase": 2, "it": it, "t": t, "obj": float(pb.d @ y), "gap_est": nu / t, "dec": dec, "step": al,
               "sum_y_minus_1": float(y.sum() - 1.0), "sec": round(time.time() - tt, 2)}
        log(row)
        hist.append(row)
        if dec < 0.5:
            if nu / t < gap:
                break
            t *= mu
            if save is not None:
                save(y, t)
    return y, t, hist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("support")
    ap.add_argument("--out", required=True)
    ap.add_argument("--gap", type=float, default=1e-9)
    ap.add_argument("--mu", type=float, default=4.0)
    ap.add_argument("--start", default=None, help="npz with a strictly feasible y (skips phase I)")
    a = ap.parse_args()
    t0 = time.time()
    pb = Problem(a.support)
    print(json.dumps({"m": pb.m, "dims": [K.n for K in pb.keys], "aut": [K.G for K in pb.keys],
                      "nnz": [K.S.nnz for K in pb.keys], "setup_s": round(time.time() - t0, 1)}), flush=True)

    def log(r):
        print(json.dumps(r), flush=True)

    def save(y, t):
        np.savez(a.out, y=y, t=t, masks=pb.masks, d=pb.d)

    if a.start:
        y = np.load(a.start)["y"]
    else:
        y, s = phase1(pb, log)
        if s >= 0:
            log({"phase1": "FAILED (no strictly feasible point found on this face)", "s": s})
            np.savez(a.out.replace(".npz", "_phase1.npz"), y=y, s=s, masks=pb.masks)
            return
        np.savez(a.out.replace(".npz", "_phase1.npz"), y=y, s=s, masks=pb.masks)
    log({"start_obj": float(pb.d @ y), "min_eigs": pb.min_eigs(y)})
    y, t, hist = phase2(pb, y, log, gap=a.gap, mu=a.mu, save=save)
    save(y, t)
    log({"final_obj": float(pb.d @ y), "t": t, "gap_est": (pb.nu_psd + pb.m) / t, "min_eigs": pb.min_eigs(y),
         "rel_min_eigs": pb.rel_min_eigs(y), "min_y": float(y.min()), "time_s": round(time.time() - t0, 1)})


if __name__ == "__main__":
    main()
