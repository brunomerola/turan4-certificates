"""Consistency tests: (1) Problem.Y equals the lpcg7-style moment matrices (restricted to live flags) for a random y;
(2) gradient and Hessian of -log det Y_k against finite differences."""
import numpy as np

import paths
from analyse_support import moment
from barrier import Problem
from common import build_model
from lpcg7 import ALL_BLOCKS
from prep_support import graph_terms


def main():
    pb = Problem(paths.RES + "/support_W.npz")
    data = np.load(paths.RES + "/support_W.npz")
    model = build_model(5, ALL_BLOCKS)
    T, _ = graph_terms(model, pb.masks)
    rng = np.random.default_rng(0)
    y = rng.random(pb.m)
    y /= y.sum()
    kk = [k for k in range(int(data["nkeys"])) if data[f"live{k}"].size]
    for K, k in zip(pb.keys, kk):
        bi, t = model.keys[k]
        Yf = moment(T[k], model.blocks[bi].aut[t], y, int(model.qdim[bi, t]))
        lv = data[f"live{k}"]
        dead = np.setdiff1d(np.arange(Yf.shape[0]), lv)
        err = np.abs(K.Y(y) - Yf[np.ix_(lv, lv)]).max()
        dz = np.abs(Yf[dead]).max() if dead.size else 0.0
        print(f"key {k}: |Y - Y_lpcg| = {err:.2e}, max |dead rows| = {dz:.2e}")
        assert err < 1e-15 and dz == 0.0
    # derivative checks on a strictly feasible-ish point: use keys 0..5 only (small, PD at random y?)
    for K, k in zip(pb.keys, kk):
        Y = K.Y(y)
        w = np.linalg.eigvalsh(Y)
        if w[0] <= 0:
            continue
        Z = np.linalg.inv(Y)
        g = -K.grad(Z)
        H = K.hess(Z, pb.m)
        dirn = rng.standard_normal(pb.m) * 1e-3 * y
        f = lambda yy: -np.linalg.slogdet(K.Y(yy))[1]
        eps = 1e-4
        fd = (f(y + eps * dirn) - f(y - eps * dirn)) / (2 * eps)
        fd2 = (f(y + dirn) - 2 * f(y) + f(y - dirn))
        print(f"key {k}: dir. derivative {g @ dirn:.8e} vs fd {fd:.8e}; 2nd {dirn @ H @ dirn:.8e} vs fd {fd2:.8e}")
    # Hessian formula with an Aut-invariant Z (Z = Y_k(y') for random y') against dense traces
    for K, k in zip(pb.keys, kk):
        Z = K.Y(rng.random(pb.m))
        H = K.hess(Z, pb.m)
        idx = rng.choice(K.graphs, size=min(6, K.graphs.size), replace=False)
        Ms = {i: K.S[i].toarray().reshape(K.n, K.n) for i in idx}
        err = max(abs(H[i, j] - np.trace(Z @ Ms[i] @ Z @ Ms[j])) / max(abs(H[i, j]), 1e-300) for i in idx for j in idx)
        print(f"key {k}: Hessian formula rel. error {err:.2e} (|Aut| = {K.G})")
        assert err < 1e-10
    print("tests done")


if __name__ == "__main__":
    main()
