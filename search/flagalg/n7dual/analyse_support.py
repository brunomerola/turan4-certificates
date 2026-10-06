"""Exploration: facial reduction of W and of supp(y_LP); spectra of the LP moment matrices on the live flags."""
import json

import numpy as np

import paths
from common import build_model
from facial import live_support
from lpcg7 import ALL_BLOCKS


def moment(T, P, y, n):
    rows, aa, bb, ww = T
    R = np.zeros((n, n))
    v = y[rows] * ww / 5040.0 * 0.5
    np.add.at(R, (aa, bb), v)
    np.add.at(R, (bb, aa), v)
    Y = np.zeros((n, n))
    for g in range(P.shape[0]):
        Y += R[np.ix_(P[g], P[g])]
    return Y


def main():
    model = build_model(5, ALL_BLOCKS)
    d = np.load(paths.RES + "/lp_y.npz")
    W, y = d["W"], d["y"]
    nk = len(model.keys)
    T = [[d[f"T{k}_{j}"] for j in range(4)] for k in range(nk)]
    P = [model.blocks[bi].aut[t] for bi, t in model.keys]
    dims = [int(model.qdim[bi, t]) for bi, t in model.keys]
    print("# facial reduction on W", flush=True)
    aliveW, liveW = live_support(T, P, dims, W.size)
    print(f"# y_LP mass on graphs killed in W: {float(y[~aliveW].sum()):.3e}, y_LP support killed: "
          f"{int(((y > 0) & ~aliveW).sum())}", flush=True)
    print("# facial reduction on supp(y_LP)", flush=True)
    aliveY, liveY = live_support(T, P, dims, W.size, alive=y > 0)
    out = {"W": int(W.size), "alive_W": int(aliveW.sum()), "live_W": [int(l.sum()) for l in liveW],
           "supp_y": int((y > 0).sum()), "alive_supp_y": int(aliveY.sum()), "live_supp_y": [int(l.sum()) for l in liveY]}
    # spectra of Y_LP on the live flags of W
    for k in range(nk):
        Y = moment(T[k], P[k], y, dims[k])
        lv = liveW[k]
        off = np.abs(Y[~lv][:, :]).max() if (~lv).any() else 0.0
        w = np.linalg.eigvalsh(Y[np.ix_(lv, lv)])
        print(f"key {k} dim {dims[k]} live {int(lv.sum())}: max |Y| on dead rows {off:.2e}; eig live: min "
              f"{w[0]:.3e}, #<-1e-12 {(w < -1e-12).sum()}, #|.|<=1e-12 {(np.abs(w) <= 1e-12).sum()}, "
              f"#>1e-9 {(w > 1e-9).sum()}, max {w[-1]:.3e}", flush=True)
    with open(paths.RES + "/analyse_support.json", "w") as f:
        json.dump(out, f, indent=1)
    np.savez(paths.RES + "/facial_W.npz", alive=aliveW, **{f"live{k}": liveW[k] for k in range(nk)})


if __name__ == "__main__":
    main()
